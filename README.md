library(dplyr)
# 1. 载入原始数据集
data(NHANES4)
# 2. 定义 CKD-EPI 2021（无种族修正版）计算函数
calc_egfr_2021 <- function(scr, age, sex) {
  egfr <- numeric(length(scr))
  for(i in seq_along(scr)) {
    s <- scr[i]
    a <- age[i]
    sx <- sex[i]
    if(is.na(s) || is.na(a) || is.na(sx)) {
      egfr[i] <- NA
      next
    }
    if(sx == 2) { # 女性
      kappa <- 0.7; alpha <- -0.241; mult <- 1.012
      egfr[i] <- 142 * (min(s/kappa, 1)^alpha) * (max(s/kappa, 1)^-1.200) * (0.9938^a) * mult
    } else { # 男性
      kappa <- 0.9; alpha <- -0.302; mult <- 1.000
      egfr[i] <- 142 * (min(s/kappa, 1)^alpha) * (max(s/kappa, 1)^-1.200) * (0.9938^a) * mult
    }
  }
  return(egfr)
}
# 3. 初始化记录初始样本量
N0 <- nrow(NHANES4)
cat("=== 队列排纳过程统计 (含生存随访过滤) ===\n")
cat("N0 (原始总样本量):", N0, "\n")
# 新增 egfr 变量
nhanes_work <- NHANES4 %>%
  mutate(egfr = calc_egfr_2021(scr = creat, age = age, sex = gender))
# 步骤 1：年龄筛选 (纳入 >= 45 岁中老年人群)
step1_data <- nhanes_work %>% filter(!is.na(age) & age >= 45)
N1 <- nrow(step1_data)
cat("N1 (筛选年龄 >= 45 岁): 剩余", N1, "人 | 剔除", N0 - N1, "人\n")
# 步骤 2：排除孕妇 (pregnant == 1)
step2_data <- step1_data %>% filter(pregnant != 1 | is.na(pregnant))
N2 <- nrow(step2_data)
cat("N2 (排除孕妇): 剩余", N2, "人 | 剔除", N1 - N2, "人\n")
# 步骤 3：剔除核心原生变量及生存随访数据缺失者 
# 注：这里保留了 phenoage0，同时去掉kdm0，并加入了 status 和 time 随访变量过滤
step3_data <- step2_data %>% filter(
  !is.na(gender),
  !is.na(creat),
  !is.na(egfr),
  !is.na(phenoage0),
  !is.na(status), # 生存状态随访数据
  !is.na(permth_int)    # 随访时间数据
)
N3 <- nrow(step3_data)
cat("N3 (剔除核心变量及随访缺失): 剩余", N3, "人 | 剔除", N2 - N3, "人\n")
# 步骤 4：剔除极端离群值 (eGFR < 15)
nhanes_analytic <- step3_data %>% filter(egfr >= 15)
N4 <- nrow(nhanes_analytic)
cat("N4 (剔除 eGFR < 15 极端值，最终分析样本): 剩余", N4, "人 | 剔除", N3 - N4, "人\n")
cat("=========================================\n")
