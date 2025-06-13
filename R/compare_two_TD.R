# # compare two TD models
library(tidyr)

# 假设已有三个结果数据框 sim1, sim2, sim3
# 下面先给出示例数据（可删掉）
library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)


#' predator-prey ODE with T dependence but without size dependence(type I)
#' N is the Population of prey and P is the Population of predator
#'
predator_prey_T <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp))  # Arrhenius公式
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp))
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp))

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P

    list(c(dN, dP),
         Temp = Temp)
    # dN  = dN,
    # dP = dP,
    # r = r,
    # a = a,
    # m = m)
  })
}
# 定义根函数：当 x<=0 或 y<=0 时触发 root = 0
root_func <- function(t, state, parameters) {
  with(as.list(state), {
    # 返回两个根，任意一个为 0 时触发事件
    c(N, P)
  })
}
event_func <- function(t, state, parameters) {
  # 将状态设置为 NA，或者直接退出。这里设为 NA，deSolve 会终止并返回
  state[] <- NA
  return(state)
}
parameters <- c(T0 =270,k = 8.618*10^(-5),n = 0.01,
                a0 = 0.3, m0 = 0.4, r0 = 0.4,
                e = 0.45, K = 100, Ea = 0.67, Em = 0.67,Er = 0.8)
state <- c(N = 50, P = 20)
times <- seq(0, 10/parameters[3], by = 1)

sim1 = as.data.frame(
  ode(y = state,
      times = times,
      func = predator_prey_T,
      parms = parameters,
      rootfun = root_func,
      events  = list(func = event_func,
                     root = TRUE))
)



pp_results_time = gather(sim1, Species, Population, -time)
pp_results_time$Temp <- parameters[1] + parameters[3]* pp_results_time$time

ggthemr('dust')
ggplot(data = pp_results_time,aes(Temp, Population, colour = Species)) +
  geom_line(size = 1)+
  ylim(0,30)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","orange"))



parameters <- c(T0 =283.15,k = 8.618*10^(-5),n = 0.001,
                a0 = 0.3, m0 = 0.5, r0 = 0.5,
                e = 0.45, K = 50, Ea = 0.67, Em = 0.67,Er = 0.67)
state <- c(N = 10, P = 20)
times <- seq(0, 600, by = 0.01)

sim2 = as.data.frame(
  ode(y = state,
      times = times,
      func = predator_prey_T,
      parms = parameters,
      rootfun = root_func)
)


parameters <- c(T0 =283.15,k = 8.618*10^(-5),n = 0.00001,
                a0 = 0.3, m0 = 0.5, r0 = 0.5,
                e = 0.45, K = 50, Ea = 0.67, Em = 0.67,Er = 0.67)
state <- c(N = 10, P = 20)
times <- seq(0, 600, by = 0.01)

sim3 = as.data.frame(
  ode(y = state,
      times = times,
      func = predator_prey_T,
      parms = parameters,
      rootfun = root_func)
)

# 添加“模拟编号”
sim1$run <- "run1"
sim2$run <- "run2"
sim3$run <- "run3"

# 合并
all_sims <- rbind(sim2, sim3)

# 转成长格式：time、run、species、Population
df_long <- pivot_longer(all_sims,
                        cols = c("P", "N"),
                        names_to  = "species",
                        values_to = "Population")

# 绘图
p <- ggplot(df_long, aes(x = time, y = Population,
                         color   = species,
                         linetype = run)) +
  geom_line(size = 0.8) +
  labs(x = "Time",
       y = "Population",
       color = "Species",
       linetype = "Simulation run") +
  theme_minimal(base_size = 14) +
  ggplot2::scale_colour_manual("Species",values = c("red2","blue2","#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","orange"))+
  theme(legend.position = "bottom")

print(p)

ggplot(df_long, aes(x = time, y = Population,
                    color   = species,
                    alpha   = run)) +
  geom_line(size = 1) +
  scale_alpha_manual(name = "Simulation run",
                     values = c("run1" = 0.3,
                                "run2" = 0.6,
                                "run3" = 1.0)) +
  labs(x = "Time",
       y = "Population",
       color = "Species") +
  theme_minimal(base_size = 14) +
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","orange"))+
  theme(legend.position = "bottom")
