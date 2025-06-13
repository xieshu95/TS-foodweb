# 安装必要包（如尚未安装）
# install.packages(c("deSolve","ggplot2","gganimate","transformr"))

library(deSolve)
library(ggplot2)
library(gganimate)

# ---- 1. 定义带温度依赖的捕食-猎物模型 ----
predprey <- function(t, state, parms) {
  N <- state[1]
  P <- state[2]
  # 举例：温度影响攻击率 a 和代谢率 m
  T <- parms$T0 + 0.01 * t  # 这里可以改为函数：parms$temp_fun(t)
  a <- parms$a0 * exp(parms$Ea_a * (T - parms$T0))
  m <- parms$m0 * exp(parms$Ea_m * (T - parms$T0))
  dN <- parms$r * N * (1 - N/parms$K) - a * N * P
  dP <- a * parms$e * N * P - m * P
  list(c(dN, dP))
}

# ---- 2. 设置参数与初始值 ----
parms <- list(
  r    = 1,     # 猎物内禀增长率
  K    = 10,    # 猎物环境承载力
  a0   = 0.5,   # 攻击率基线
  Ea_a = 0.05,  # 攻击率温度灵敏度
  m0   = 0.2,   # 代谢率基线
  Ea_m = 0.03,  # 代谢率温度灵敏度
  e    = 0.1,   # 转化效率
  T0   = 20,    # 参考温度
  temp = 25     # 恒定温度示例
)
state0 <- c(N = 5, P = 2)

# ---- 3. 生成数值解 ----
times <- seq(0, 100, by = 0.5)
out <- ode(y = state0, times = times, func = predprey, parms = parms)
df  <- as.data.frame(out)

# ---- 4. 时间序列动画 ----
p1 <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = N, colour = "Prey")) +
  geom_line(aes(y = P, colour = "Predator")) +
  labs(y = "Density", colour = "") +
  theme_minimal() +
  transition_reveal(time)

# 渲染动画（可输出 gif 或 mp4）
anim1 <- animate(p1, nframes = 200, fps = 20, width = 600, height = 400)
# anim_save("time_series.gif", anim1)

# ---- 5. 相平面动画 ----
# 我们保留整个轨迹并随着时间动态高亮当前点
df_phase <- df
df_phase$current <- seq_len(nrow(df_phase))

p2 <- ggplot(df_phase, aes(x = N, y = P)) +
  geom_path(alpha = 0.3) +
  geom_point(aes(frame = current), size = 2, colour = "red") +
  labs(x = "Prey (N)", y = "Predator (P)") +
  theme_minimal() +
  transition_time(time) +
  ease_aes('linear')

anim2 <- animate(p2, nframes = 200, fps = 20, width = 600, height = 400)
# anim_save("phase_plane.gif", anim2)

# ---- 结束 ----
# 以上 anim1, anim2 就是两个动态可视化对象，
# 可以直接在 RStudio Viewer 中播放，或使用 anim_save() 导出。




