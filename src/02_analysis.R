# 概率论与数理统计X，2026秋，大作业：两个实验训练
# 运行：Rscript src/02_analysis.R [项目根目录]
# 仅使用 R 自带的 base/stats/graphics/grDevices，无需安装第三方 R 包。
# 与 src/01_extract_bin.py 提取的实际 BIN 日志数据配合。
# 本脚本读取固定 CSV 输入，重新计算结果及绘图；用户无须自己运行。

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
in_data <- function(file) file.path(root, "data", file)
out_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================= 实验训练1 =======================
# 统计对象：约25m相对Home高度、接近目标航向的两段近似直线飞行。
# 具体筛选标准和原始日志解码在 01_extract_bin.py 中，R读取已审计的结果。
samples <- read.csv(in_data("airspeed_samples.csv"), stringsAsFactors = FALSE,
                    check.names = FALSE, fileEncoding = "UTF-8")
stopifnot(nrow(samples) >= 30, all(is.finite(samples$airspeed_mps)))
x <- samples$airspeed_mps
n <- length(x)
mean_x <- mean(x)
ordered_x <- sort(x, decreasing = FALSE)
median_x <- median(x)
var_x <- var(x)  # 分母 n-1，无偏样本方差
sd_x <- sd(x)
q1_x <- unname(quantile(x, 0.25))
q3_x <- unname(quantile(x, 0.75))

stat_names <- c("sample_n", "mean_mps", "median_mps", "sample_variance_mps2",
                "sample_sd_mps", "min_mps", "q1_mps", "q3_mps", "max_mps",
                "coefficient_of_variation", "iqr_mps")
stat_values <- c(n, mean_x, median_x, var_x, sd_x,
                 min(x), q1_x, q3_x, max(x), sd_x / mean_x, IQR(x))
write.csv(data.frame(metric = stat_names, value = stat_values),
          file.path(out_dir, "experiment1_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(rank = seq_along(ordered_x), airspeed_mps = ordered_x),
          file.path(out_dir, "experiment1_order_statistics.csv"), row.names = FALSE)

# 分段统计用于描述；切勿把这些连续高频数据当作相互独立的重复试验。
segments <- split(x, samples$segment_id)
segment_summary <- do.call(rbind, lapply(names(segments), function(nm) {
  values <- segments[[nm]]
  data.frame(segment_id = nm, n = length(values), mean_mps = mean(values),
             median_mps = median(values), sample_sd_mps = sd(values),
             min_mps = min(values), max_mps = max(values))
}))
write.csv(segment_summary, file.path(out_dir, "experiment1_by_segment.csv"), row.names = FALSE)

# 频率密度直方图：每个柱的面积等于该组的相对频率。
# 为便于复核，组距固定为 0.5 m/s，而不依赖 hist() 的自动分组。
breaks_x <- seq(floor(min(x) * 2) / 2,
                ceiling(max(x) * 2) / 2 + 0.5, by = 0.5)
png(file.path(fig_dir, "02_airspeed_histogram_R.png"), width = 2100, height = 1350, res = 240)
par(mar = c(4.5, 5.1, 3.1, 1.5))
hist(x, breaks = breaks_x, freq = FALSE, col = "#bfd8e7", border = "white",
     xlab = "Measured airspeed (m/s)", ylab = "Relative-frequency density (s/m)",
     main = sprintf("Near-level straight-leg airspeed (n = %d)", n), cex.lab = 1.15)
grid(nx = NA, ny = NULL, col = "#e5ebee", lty = 3)
curve(dnorm(x, mean = mean_x, sd = sd_x), col = "#be6b4e", lwd = 2.4, add = TRUE)
abline(v = mean_x, col = "#1c6585", lty = 2, lwd = 2)
legend("topright", c("Histogram (0.5 m/s bin)", "Matched normal reference", "Sample mean"),
       col = c("#bfd8e7", "#be6b4e", "#1c6585"),
       lwd = c(9, 2.4, 2), lty = c(1, 1, 2), bty = "n")
dev.off()

# 经验分布函数（阶梯图），正态分布函数仅作参考拟合，非正态性检验。
png(file.path(fig_dir, "03_empirical_cdf_R.png"), width = 2100, height = 1350, res = 240)
par(mar = c(4.5, 5.1, 3.1, 1.5))
plot(ecdf(x), do.points = FALSE, verticals = TRUE, pch = NA,
     xlab = "Measured airspeed (m/s)", ylab = "Empirical cumulative probability",
     main = "Empirical distribution function", col = "#1c6585", lwd = 2.2,
     xlim = range(x) + c(-0.35,0.35))
grid(col = "#e5ebee", lty = 3)
xx <- seq(min(x)-0.35, max(x)+0.35, length.out = 500)
lines(xx, pnorm(xx, mean = mean_x, sd = sd_x), col = "#be6b4e", lwd = 2.3, lty = 2)
legend("topleft", c("Empirical CDF", "Matched normal reference"),
       col = c("#1c6585", "#be6b4e"), lty = c(1,2), lwd = c(2.2,2.3), bty = "n")
dev.off()

# 定高直航段的观测时序，不将两段数据跨时间连起来。
png(file.path(fig_dir, "04_selected_timeseries_R.png"), width = 2100, height = 1350, res = 240)
par(mar = c(4.5, 5.1, 3.1, 1.5))
plot(samples$time_s, samples$airspeed_mps, type = "n",
     xlab = "Flight controller elapsed time (s)", ylab = "Measured airspeed (m/s)",
     main = "Selected airspeed records by flight leg", ylim = range(x) + c(-0.5,0.5))
grid(col = "#e5ebee", lty = 3)
for (seg in unique(samples$segment_id)) {
  part <- samples[samples$segment_id == seg, , drop = FALSE]
  color <- if (seg == "S") "#1c6585" else "#c17c4d"
  lines(part$time_s, part$airspeed_mps, col = color, lwd = 1.8, type = "b", pch = 16, cex = 0.65)
}
abline(h = mean_x, col = "#74828a", lty = 3, lwd = 1.8)
legend("topright", c("Southbound: WP3-WP4", "Northbound: WP7-WP9", "Pooled mean"),
       col = c("#1c6585", "#c17c4d", "#74828a"),
       lty = c(1,1,3), lwd = 2, bty = "n")
dev.off()

# 二段对比，仅描述各段分布，不作推断性显著性检验。
png(file.path(fig_dir, "05_segment_boxplot_R.png"), width = 1800, height = 1250, res = 240)
par(mar = c(4.5, 5.1, 3.1, 1.5))
boxplot(airspeed_mps ~ segment_id, data = samples, names = c("Northbound", "Southbound"),
        col = c("#f1dac6", "#bcd6e6"), border = "#6a7580",
        main = "Airspeed by selected straight leg", ylab = "Measured airspeed (m/s)",
        outline = TRUE)
grid(nx = NA, col = "#e5ebee", lty = 3)
stripchart(airspeed_mps ~ segment_id, data = samples, vertical = TRUE, method = "jitter",
           add = TRUE, jitter = 0.1, pch = 16, cex = 0.65, col = "#354d5b")
dev.off()

# ================ 实验训练2：正态总体推断 ================
A <- c(88, 87, 92, 90, 91)
B <- c(89, 89, 90, 84, 88)
alpha <- 0.05
confidence_level <- 1 - alpha

# 1. 方差未知时对甲总体均值的 95% t 区间
mean_CI <- t.test(A, conf.level = confidence_level)$conf.int
# 2. 两个独立正态总体的方差比 95% F 区间，方向为甲/乙
ratio_test <- var.test(A, B, ratio = 1, conf.level = confidence_level)
ratio_CI <- ratio_test$conf.int
# 3. 检验 H0: mu_A = 90；H1: mu_A != 90 （双侧）
t_result <- t.test(A, mu = 90, alternative = "two.sided", conf.level = confidence_level)

rows <- data.frame(
  metric = c("n_A", "n_B", "mean_A", "mean_B", "sample_var_A", "sample_var_B",
             "mean_A_ci95_lower", "mean_A_ci95_upper",
             "variance_ratio_A_over_B", "variance_ratio_ci95_lower", "variance_ratio_ci95_upper",
             "t_statistic_muA_eq_90", "t_df", "t_pvalue_two_sided", "alpha"),
  value = c(length(A), length(B), mean(A), mean(B), var(A), var(B),
            mean_CI[1], mean_CI[2], unname(ratio_test$estimate), ratio_CI[1], ratio_CI[2],
            unname(t_result$statistic), unname(t_result$parameter), t_result$p.value, alpha)
)
write.csv(rows, file.path(out_dir, "experiment2_results.csv"), row.names = FALSE)

cat("Exercise 1: n=", n, ", mean=", mean_x, ", median=", median_x,
    ", var=", var_x, ", sd=", sd_x, "\n", sep = "")
cat("Exercise 1: ordered sample =", paste(sprintf("%.6f", ordered_x), collapse = ", "), "\n")
cat("Exercise 2: mean 95% CI =", paste(signif(mean_CI, 9), collapse = " to "), "\n")
cat("Exercise 2: variance-ratio 95% CI =", paste(signif(ratio_CI, 9), collapse = " to "), "\n")
print(t_result)
if (t_result$p.value < alpha) {
  cat("Decision: reject H0 at alpha=0.05.\n")
} else {
  cat("Decision: fail to reject H0 at alpha=0.05.\n")
}
