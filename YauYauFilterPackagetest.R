# ---------------------------------------------------------------------------- #
#                                    清除环境
# ---------------------------------------------------------------------------- #

clean_environment <- function() {
  # 清除全局环境中的所有对象（包括变量和函数）
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  
  # 卸载所有非基础包
  base_packages <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
  loaded_packages <- search()[grepl("^package:", search())]
  to_detach <- setdiff(loaded_packages, paste0("package:", base_packages))
  lapply(to_detach, detach, character.only = TRUE, unload = TRUE)
  
  # 清除图形设备
  while (dev.cur() > 1) dev.off()
  
  # 清空 Console
  cat("\014")
}
clean_environment()

# ---------------------------------------------------------------------------- #
#                               初始化参数
# ---------------------------------------------------------------------------- #

library(YauYauFilter)
ls("package:YauYauFilter")

required_packages <- c(
  "usethis", 
  "devtools", 
  "roxygen2", 
  "roxygen2md", 
  "pkgdown", 
  "testthat", 
  "Rcpp", 
  "RcppArmadillo", 
  "RcppEigen", 
  "Matrix", 
  "Deriv", 
  "reshape2", 
  "gridExtra", 
  "viridis"
)
lapply(required_packages, library, character.only = TRUE)


Dim <- 2
T <- 20
Dt <- 0.001
Dtau <- 5 * Dt
Ds <- 0.5
Nt <- as.integer(Dtau / Dt)
Ntau <- as.integer(T / Dtau)
NtNtau <- as.integer(T / Dt)
# 
f <- function(x) {return(c(cos(x[1]), cos(x[2])))}
h <- function(x) {return(c(x[1]^3, x[2]^3))}

# f <- function(x) {return(c(cos(x[1]), cos(x[2]), cos(x[3])))}
# h <- function(x) {return(c(x[1]^3-x[1]^2, x[2]^3-x[2]^2, x[3]^3-x[3]^2))}


df <- generate_derivative(f)


# ---------------------------------------------------------------------------- #
#                             模拟状态和观测
# ---------------------------------------------------------------------------- #

library("RcppArmadillo")
seed_value <- 42
result <- Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed = seed_value)
x <- result$x
y <- result$y

source("~/YauYauFilterRfunction/plot_simulation.R")
plot_combined(x, y, Dim, T, NtNtau)


# ---------------------------------------------------------------------------- #
#                                   离散化
# ---------------------------------------------------------------------------- #


s <- seq(min(x), max(x)+Ds, by = Ds)
Ns <- length(s)
s <- ExpandGrid(Dim, s)
D <- generateD(Dim, Ns, Ds)
B <- computeB(s,D,Dt,Ds,f,df,h)
Lambda <- computeLambda(Dim, Ns, Dt, Ds)


# ---------------------------------------------------------------------------- #
#                                   离散化
# ---------------------------------------------------------------------------- #


Iu <- wrap_outiu_function(s, NtNtau, Ntau, Nt, Dim, y, h, Lambda, B, Ns, NormalizedExp, DST_Solver)


plot_Iu <- function(x, Iu, Dim, T, NtNtau) {
  # 设置绘图布局，增加右侧外边距，为图例预留空间
  par(mfrow = c(ceiling(Dim / 3), min(3, Dim)), oma = c(4, 4, 1, 8), mar = c(3, 3, 2, 1))
  
  for (i in 1:Dim) {
    # 获取当前维度的 x 和 Iu 范围，用于限定网格线
    x_range <- range(1:NtNtau)
    y_range <- range(c(x[, i], Iu[, i]))
    
    # 绘制 x 和 Iu 的折线图
    plot(x[, i], type = "l", col = "blue", lty = 1, lwd = 2, 
         xlab = "", ylab = "", main = paste("Dimension", i), xaxt = "n", yaxt = "n",
         xlim = x_range, ylim = y_range)
    lines(Iu[, i], col = "red", lty = 1, lwd = 2)  # Iu 使用实线并保持红色区分
    
    # 手动添加网格线，仅限于子图范围内
    x_ticks <- pretty(x_range)  # 获取 x 轴网格刻度
    y_ticks <- pretty(y_range)  # 获取 y 轴网格刻度
    
    # 水平网格线
    segments(x0 = min(x_range), x1 = max(x_range), y0 = y_ticks, y1 = y_ticks, col = "gray", lty = "dotted")
    # 垂直网格线
    segments(x0 = x_ticks, x1 = x_ticks, y0 = min(y_range), y1 = max(y_range), col = "gray", lty = "dotted")
    
    # 添加坐标轴
    axis(1, at = seq(0, NtNtau, length.out = 5), labels = round(seq(0, T, length.out = 5), 2))
    axis(2)
  }
  
  # 添加公共横纵轴标签
  mtext("Time", side = 1, outer = TRUE, line = 0.5, cex = 1)
  mtext("Value", side = 2, outer = TRUE, line = 0.5, cex = 1)
  
  # 在右侧绘制统一图例，进一步右移
  par(xpd = NA)  # 允许绘图超出边界
  legend("topright", inset = c(-0.4, 0), legend = c("x", "Iu"), 
         col = c("blue", "red"), lty = c(1, 1), lwd = c(2, 2), bty = "n", cex = 1.2)
}



# 调用绘图函数
plot_Iu(x, Iu, Dim, T, NtNtau)






































































