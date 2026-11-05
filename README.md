library(dplyr)
library(mclust)
library(tableone)
library(ggplot2)
library(survminer)
library(survival)
library(tidyr)
# 1. 数据清洗与预处理
data(NHANES4)
#定义 CKD-EPI 2021 计算函数
calc_egfr_2021 <- function(scr, age, sex) {
  egfr <- numeric(length(scr))
  for(i in seq_along(scr)) {
    s <- scr[i]; a <- age[i]; sx <- sex[i]
    if(is.na(s) || is.na(a) || is.na(sx)) { egfr[i] <- NA; next }
    if(sx == 2) {
      egfr[i] <- 142 * (min(s/0.7, 1)^-0.241) * (max(s/0.7, 1)^-1.200) * (0.9938^a) * 1.012
    } else {
      egfr[i] <- 142 * (min(s/0.9, 1)^-0.302) * (max(s/0.9, 1)^-1.200) * (0.9938^a) * 1.000
    }
  }
  return(egfr)
}

nhanes_work <- NHANES4 %>%
  mutate(
    egfr = calc_egfr_2021(scr = creat, age = age, sex = gender),
    new_hyperten = if_else(sbp >= 140 | dbp >= 90 | hyperten == 1, 1, 0, missing = 0),
    new_diabetes = if_else(hba1c > 6.5 | glucose_fasting >= 7.0 | diabetes == 1, 1, 0, missing = 0)
  )
#规范排纳流程
#0. 记录初始总样本量
N0 <- nrow(nhanes_work)
cat("=========================================\n")
cat("N0 (初始原始样本量):", N0, "人\n")
#步骤 1：年龄筛选 (纳入 >= 45 岁中老年人群)
step1_data <- nhanes_work %>% filter(!is.na(age) & age >= 45)
N1 <- nrow(step1_data)
cat("N1 (筛选年龄 >= 45 岁): 剩余", N1, "人 | 剔除", N0 - N1, "人\n")
#步骤 2：排除孕妇 (pregnant == 1)
step2_data <- step1_data %>% filter(pregnant != 1 | is.na(pregnant))
N2 <- nrow(step2_data)
cat("N2 (排除孕妇): 剩余", N2, "人 | 剔除", N1 - N2, "人\n")
#步骤 3：剔除核心聚类变量及生存随访数据缺失者
#注：严格纳入 phenoage0、egfr、bmi、sbp 以及生存随访和基础协变量
step3_data <- step2_data %>% filter(
  !is.na(gender),
  !is.na(creat),
  !is.na(egfr),
  !is.na(phenoage0),
  !is.na(bmi),          # 补全聚类核心变量 bmi 缺失过滤
  !is.na(sbp),          # 补全聚类核心变量 sbp 缺失过滤
  !is.na(status),       # 生存状态随访数据
  !is.na(permth_int)    # 随访时间数据
)
N3 <- nrow(step3_data)
cat("N3 (剔除核心变量及随访缺失): 剩余", N3, "人 | 剔除", N2 - N3, "人\n")
#步骤 4：剔除极端离群值 (eGFR < 15)
nhanes_analytic <- step3_data %>% filter(egfr >= 15)
N4 <- nrow(nhanes_analytic)
cat("N4 (剔除 eGFR < 15 极端值，最终分析样本): 剩余", N4, "人 | 剔除", N3 - N4, "人\n")
cat("=========================================\n")

# 2. 基于 PhenoAge 的 GMM 聚类建模
#提取聚类核心变量（仅保留 phenoage0 及代谢/肾功能指标）
cluster_vars <- c("phenoage0", "egfr", "bmi", "sbp")
cluster_matrix <- scale(nhanes_analytic[, cluster_vars])
#计算不同群数 (G) 和协方差模型下的 BIC 值
bic_values <- mclustBIC(cluster_matrix)
gmm_fit$BIC
#最终锁定 G = 4 聚类
set.seed(123) 
optimal_G <- 4 
gmm_fit_final <- Mclust(cluster_matrix, G = optimal_G)
cat("=== 最终拟合模型摘要 ===\n")
print(summary(gmm_fit_final))
nhanes_analytic_clustered <- nhanes_analytic
nhanes_analytic_clustered$cluster_gmm <- as.factor(gmm_fit_final$classification)
cat("=== 4 组亚型人数最终分布 ===\n")
print(table(nhanes_analytic_clustered$cluster_gmm, useNA = "ifany"))

# 3. Table 1 人群基线表生成
vars_to_compare <- c(
  "age", "gender", "race", "education", "poverty_ratio", 
  "bmi", "sbp", "dbp", "new_hyperten", "new_diabetes", 
  "phenoage0", "egfr"
)
nhanes_analytic_clustered <- nhanes_analytic_clustered %>%
  mutate(
    gender = as.factor(gender),
    race = as.factor(race),
    education = as.factor(education),
    new_hyperten = as.factor(new_hyperten),
    new_diabetes = as.factor(new_diabetes),
    cluster_gmm = as.factor(cluster_gmm)
  )
categorical_vars <- c("gender", "race", "education", "new_hyperten", "new_diabetes")
table1_obj <- CreateTableOne(
  vars = vars_to_compare,
  strata = "cluster_gmm",
  data = nhanes_analytic_clustered,
  factorVars = categorical_vars,
  test = TRUE
)
print(table1_obj, showAllLevels = TRUE, smd = TRUE)
table1_df <- print(table1_obj, printToggle = FALSE, showAllLevels = TRUE)
write.csv(table1_df, file = "Baseline_Table1_PhenoAge_Clusters.csv", row.names = TRUE)

# 4. K-M 生存曲线分析
library(survival)
library(survminer)
library(ggplot2)
#确保因子水平严格为 1, 2, 3, 4 自然顺序
nhanes_analytic_clustered$cluster_gmm <- factor(
  as.character(nhanes_analytic_clustered$cluster_gmm),
  levels = c("1", "2", "3", "4")
)
#拟合生存模型
km_fit <- survfit(Surv(permth_int, status) ~ cluster_gmm, data = nhanes_analytic_clustered)
#顶刊经典色系 (Lancet/Nature：海蓝、琥珀黄、深橙红、祖母绿)
sci_palette <- c("#0072B2", "#E69F00", "#D55E00", "#009E73")
# 绘制 K-M 曲线（图例与统计值精准置于左下角纯白区）
km_plot_topjournal <- ggsurvplot(
  km_fit,
  data = nhanes_analytic_clustered,
    # --- 曲线与置信区间 ---
  size = 1.0,               # 典雅适中的主曲线粗细
  censor = FALSE,           # 关键：大样本不画密集十字毛刺
  conf.int = TRUE,          # 保留置信区间
  conf.int.alpha = 0.12,    # 薄雾通透阴影，不显脏乱
    # --- 坐标轴 ---
  ylim = c(0.75, 1.00),     
  break.x.by = 50,          # 0, 50, 100, 150, 200 个月标准刻度
  xlab = "Time (Months)", 
  ylab = "Survival Probability",
    # --- 统计学显著性标注 ---
  pval = "Log-rank p < 0.001",
  pval.coord = c(10, 0.86),    # 放置在左下角上方空白处 (Time=10, Y=0.86)
  pval.size = 4.2,
    # --- 图例（关键修改：移动到左下角完全空白区） ---
  palette = sci_palette,
  legend = c(0.18, 0.23),      # 相对坐标 (X: 18%, Y: 23%)，精准落入左下角无曲线区域
  legend.title = "",           # 去除冗余的 Subtypes 标题
  legend.labs = c("Subtype 1", "Subtype 2", "Subtype 3", "Subtype 4"),
    # --- 风险表 (Risk Table) 样式 ---
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.22,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  fontsize = 3.6,
    # --- 经典 L 型顶刊极简主题 ---
  ggtheme = theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      axis.text = element_text(size = 11, face = "bold", color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      legend.text = element_text(size = 10, face = "bold"),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA), # 半透明纯白底
      plot.margin = margin(t = 10, r = 15, b = 5, l = 15)
    )
)
#微调风险表去除多余边框
km_plot_topjournal$table <- km_plot_topjournal$table + 
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 11, face = "bold", hjust = 0)
  )
#打印查看
print(km_plot_topjournal)
#导出发表级 PDF 矢量图与高分 TIFF
pdf("KM_Survival_TopJournal_Clean.pdf", width = 7.5, height = 6.5, onefile = FALSE)
print(km_plot_topjournal)
dev.off()
tiff("KM_Survival_TopJournal_Clean.tiff", width = 2250, height = 1950, res = 300, compression = "lzw")
print(km_plot_topjournal)
dev.off()

# 5. 单因素 Cox 回归与 PH 假定检验
library(survival)
library(dplyr)
#【核心设定】：显式将 Subtype 3（最健康/预后最好组）设定为参照基准 (ref = "3", HR = 1.00)
nhanes_analytic_clustered$cluster_gmm <- relevel(factor(nhanes_analytic_clustered$cluster_gmm), ref = "3")
#拟合单因素 Cox 比例风险回归模型
uni_cox_model <- coxph(
  Surv(permth_int, status) ~ cluster_gmm, 
  data = nhanes_analytic_clustered
)
#提取详细结果并构建标准 SCI 论文格式的表格
uni_summary <- summary(uni_cox_model)
uni_hr <- exp(uni_summary$coefficients[, "coef"])
uni_ci <- exp(confint(uni_cox_model))
uni_pval <- uni_summary$coefficients[, "Pr(>|z|)"]
uni_results_table <- data.frame(
  Comparison = rownames(uni_summary$coefficients),
  HR = round(uni_hr, 3),
  Lower_95CI = round(uni_ci[, 1], 3),
  Upper_95CI = round(uni_ci[, 2], 3),
  P_value = signif(uni_pval, 3)
)
rownames(uni_results_table) <- NULL
cat("=== 1. 单因素 Cox 回归结果汇总表 (Ref = Subtype 3) ===\n")
print(uni_results_table)
#导出为表格文件
write.csv(uni_results_table, file = "Univariate_Cox_PhenoAge_Results.csv", row.names = FALSE)
#运行比例风险假定检验 (cox.zph)
ph_test <- cox.zph(uni_cox_model)
cat("\n=== 2. 等比例风险假定 (cox.zph) 统计检验结果 ===\n")
print(ph_test)
#导出 Schoenfeld 残差诊断图
pdf("Schoenfeld_Residuals_Plot.pdf", width = 7, height = 6)
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
plot(
  ph_test, 
  col = "black", 
  lwd = 1.5, 
  lty = c(1, 2, 2), # 实线为平滑拟合曲线，虚线为 95% 置信区间
  xlab = "Follow-up Time (Months)", 
  ylab = "Beta(t) for PhenoAge Subtypes"
)
abline(h = 0, col = "red", lty = 2, lwd = 1)
title(main = "Schoenfeld Residuals Test (PH Assumption)")
dev.off()
cat("\n=== 提示：Schoenfeld 诊断图已成功保存！===\n\n")

# 核心指标均值画像与顶刊热图绘制
library(dplyr)
library(tidyr)
library(ggplot2)
#提取各亚型在 4 个核心指标上的均值
cluster_vars_profile <- c("phenoage0", "egfr", "bmi", "sbp")
subtype_profile <- nhanes_analytic_clustered %>%
  group_by(cluster_gmm) %>%
  summarise(across(all_of(cluster_vars_profile), ~ mean(.x, na.rm = TRUE)))
#进行标准化处理，并对保护性指标 eGFR 实施方向反转（取负）
subtype_profile_aligned <- subtype_profile %>%
  mutate(
    z_PhenoAge = scale(phenoage0)[, 1],
    z_eGFR_impairment = -scale(egfr)[, 1],  # 核心：取负值，代表肾功能减退/损伤程度
    z_BMI = scale(bmi)[, 1],
    z_SBP = scale(sbp)[, 1]
  ) %>%
  select(cluster_gmm, z_PhenoAge, z_eGFR_impairment, z_BMI, z_SBP)
#转换为长数据格式
heatmap_long <- subtype_profile_aligned %>%
  pivot_longer(
    cols = starts_with("z_"),
    names_to = "Variable",
    values_to = "Z_Score"
  )
#规范化指标标签与列排列顺序
heatmap_long$Variable <- factor(
  heatmap_long$Variable,
  levels = c("z_PhenoAge", "z_eGFR_impairment", "z_BMI", "z_SBP"),
  labels = c("PhenoAge", "Renal Impairment\n(-eGFR)", "BMI", "Systolic BP")
)
#设定纵轴阶梯顺序：由下至上由健康走向高危
#Subtype 3 (底/最健康) -> Subtype 4 -> Subtype 1 -> Subtype 2 (顶/最危重)
heatmap_long$Subtype <- paste0("Subtype ", heatmap_long$cluster_gmm)
heatmap_long$Subtype <- factor(
  heatmap_long$Subtype, 
  levels = c("Subtype 3", "Subtype 4", "Subtype 1", "Subtype 2")
)
#绘制临床语义统一对齐的顶刊热图
p_heatmap_aligned <- ggplot(heatmap_long, aes(x = Variable, y = Subtype, fill = Z_Score)) +
  geom_tile(color = "white", linewidth = 1.2, na.rm = FALSE) +
  geom_text(aes(label = sprintf("%.2f", Z_Score)), color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",      # 优质深蓝（健康/生理机能优异）
    mid = "#FFFFBF",      # 柔和米黄（全人群平均水准）
    high = "#B2182B",     # 优质深红（病理加重/高度危险）
    midpoint = 0,
    limit = c(-1.5, 1.5), 
    name = "Risk\nZ-Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, lineheight = 0.9),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 11),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Clinical Risk & Phenotypic Aging Profile Across GMM Subtypes"
  )
#控制台预览与高清导出
print(p_heatmap_aligned)
ggsave("Publication_Quality_Heatmap_Aligned.pdf", plot = p_heatmap_aligned, width = 8.5, height = 5, dpi = 300)
cat("=== 极性对齐版顶刊热图已生成并导出为 PDF 文件！===\n")
cat("\n=== 提示：单因素回归表格与规范热图已全部生成完毕！===\n")
