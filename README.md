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
nhanes_analytic_clustered$cluster_gmm <- relevel(nhanes_analytic_clustered$cluster_gmm, ref = "4")
km_fit <- survfit(Surv(permth_int, status) ~ cluster_gmm, data = nhanes_analytic_clustered)
km_plot_optimized <- ggsurvplot(
  km_fit,
  data = nhanes_analytic_clustered,
  pval = TRUE,             
  pval.coord = c(5, 0.82), 
  conf.int = TRUE,         
  ylim = c(0.75, 1.00),    
  censv = TRUE,            
  surv.median.line = "none", 
  risk.table = TRUE,       
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  palette = c("#2E9FDF", "#E7B800", "#FC4E07", "#2ca02c"),
  xlab = "Follow-up Time (Months)", 
  ylab = "Survival Probability",
  legend.title = "Subtypes",
  legend.labs = c("Subtype 1", "Subtype 2", "Subtype 3", "Subtype 4"),
  legend = c(0.85, 0.25), 
  ggtheme = theme_bw() + theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.8))
  )
)
print(km_plot_optimized)
ggsave("KM_Survival_Curve_PhenoAge.pdf", plot = print(km_plot_optimized), width = 8, height = 7)

# 5. 单因素 Cox 回归与 PH 假定检验
# ==========================================
# 1. 单因素 Cox 回归分析与比例风险假定 (PH) 检验
# ==========================================
library(survival)
library(dplyr)

# 【核心设定】：显式将 Subtype 4（最健康/预后最好组）设定为参照基准 (ref = "4", HR = 1.00)
nhanes_analytic_clustered$cluster_gmm <- relevel(factor(nhanes_analytic_clustered$cluster_gmm), ref = "4")

# 拟合单因素 Cox 比例风险回归模型
uni_cox_model <- coxph(
  Surv(permth_int, status) ~ cluster_gmm, 
  data = nhanes_analytic_clustered
)

# 提取详细结果并构建标准 SCI 论文格式的表格
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

cat("=== 1. 单因素 Cox 回归结果汇总表 (Ref = Subtype 4) ===\n")
print(uni_results_table)

# 导出为表格文件
write.csv(uni_results_table, file = "Univariate_Cox_PhenoAge_Results.csv", row.names = FALSE)

# 运行比例风险假定检验 (cox.zph)
ph_test <- cox.zph(uni_cox_model)
cat("\n=== 2. 等比例风险假定 (cox.zph) 统计检验结果 ===\n")
print(ph_test)

# 导出 Schoenfeld 残差诊断图 (PDF 矢量图)
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


# ==========================================
# 2. 核心指标均值画像与顶刊热图绘制
# ==========================================
library(tidyr)
library(ggplot2)

# 定义核心指标变量
cluster_vars_profile <- c("phenoage0", "egfr", "bmi", "sbp")

# 计算各亚型在核心指标上的均值
subtype_profile <- nhanes_analytic_clustered %>%
  group_by(cluster_gmm) %>%
  summarise(across(all_of(cluster_vars_profile), ~ mean(.x, na.rm = TRUE)))

# 进行 Z-score 标准化（均值为0，方差为1）
subtype_profile_scaled <- subtype_profile %>%
  mutate(across(all_of(cluster_vars_profile), ~ scale(.x)[, 1]))

# 调整绘图因子顺序：将最健康的 Subtype 4 排在最下方，高危组在上方
subtype_profile_scaled$cluster_gmm <- factor(
  subtype_profile_scaled$cluster_gmm, 
  levels = c("2", "3", "1", "4") # 4在最下，2在最上
)

# 转换为长数据格式
heatmap_long <- subtype_profile_scaled %>%
  pivot_longer(
    cols = all_of(cluster_vars_profile),
    names_to = "Variable",
    values_to = "Z_Score"
  )

# 规范化变量显示名称
heatmap_long$Variable <- factor(
  heatmap_long$Variable,
  levels = c("phenoage0", "egfr", "bmi", "sbp"),
  labels = c("PhenoAge", "eGFR", "BMI", "Systolic BP")
)

heatmap_long$Subtype <- paste0("Subtype ", heatmap_long$cluster_gmm)
heatmap_long$Subtype <- factor(
  heatmap_long$Subtype, 
  levels = c("Subtype 2", "Subtype 3", "Subtype 1", "Subtype 4")
)

# 绘制高颜值、顶刊风格的精细化热图
p_heatmap_pro <- ggplot(heatmap_long, aes(x = Variable, y = Subtype, fill = Z_Score)) +
  geom_tile(color = "white", linewidth = 1.2, na.rm = FALSE) +
  geom_text(aes(label = sprintf("%.2f", Z_Score)), color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",      # 优质深蓝（代表低风险/负向偏离）
    mid = "#FFFFBF",      # 柔和过渡白黄色
    high = "#B2182B",     # 优质深红（代表高风险/正向偏离）
    midpoint = 0,
    limit = c(-1.5, 1.5), 
    name = "Z-Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 11),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Clinical and PhenoAge Profile Across GMM Subtypes"
  )

# 在控制台打印预览
print(p_heatmap_pro)

# 导出为高分辨率出版级 PDF 矢量图
ggsave("Publication_Quality_Heatmap_PhenoAge.pdf", plot = p_heatmap_pro, width = 8, height = 5, dpi = 300)
cat("\n=== 提示：单因素回归表格与规范热图已全部生成完毕！===\n")
