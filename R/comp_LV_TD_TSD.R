# compare PP vs TDpp and TSDpp model (三个模型完全一致的条件)
library(tidyr)
library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)

predator_prey_TI <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # update temperature
    # temperature-dependent rates
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pr)   # Arrhenius公式
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pa1) * (z2 ^ pa2)
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * (z2 ^ pm)

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P

    list(c(dN, dP),
         z1 = z1,
         z2 = z2,
         Temp = Temp,
         dN  = dN,
         dP = dP,
         r = r,
         a = a,
         m = m)
  })
}

root_func <- function(t, state, parms) {
  with(as.list(state), {
    c(N, P)
  })
}
event_func <- function(t, state, parameters) {
  # 将状态设置为 NA，或者直接退出。这里设为 NA，deSolve 会终止并返回
  state[] <- NA
  return(state)
}
parms <- c(T0 =273.15,k = 8.618*10^(-5),n = 0.01,Temp = 283.15,
           a0 = 0.5, m0 = 0.5, r0 = 0.5,
           e = 0.45, K = 50, Ea = 0.2, Em = 0.67,Er = 0.67,
           pr = -0.25, pa1 = -0.25, pa2 = -0.25, pm = -0.25, #pr影响
           mu1 = 1, mu2 = 1, z1 = 1, z2 = 1) #决定N P body size变化幅度
state <- c(N = 20, P = 10)
times <- seq(0, 500, by = 0.1)

sim1 = as.data.frame(
  ode(y = state,
      times = times,
      func = predator_prey_TI,
      parms = parms,
      rootfun = root_func,
      events  = list(func = event_func,
                     root = TRUE))
)



#' predator-prey ODE with T dependence but without size dependence(type I)
#' N is the biomass of prey and P is the biomass of predator
#'
predator_prey_T <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    # a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp))  # Arrhenius公式
    # m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp))
    # r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp))
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pr)
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pa1) * (z2 ^ pa2)  # Arrhenius公式
    # metabolic rate of the predator
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * (z2 ^ pm)

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P

    list(c(dN, dP),
         z1 = z1,
         z2 = z2,
         Temp = Temp,
         dN  = dN,
         dP = dP,
         r = r,
         a = a,
         m = m)
  })
}

root_func <- function(t, state, parms) {
  with(as.list(state), {
    c(N, P)
  })
}
parms <- c(T0 =273.15,k = 8.618*10^(-5),n = 0.02,
           a0 = 0.5, m0 = 0.5, r0 = 0.5,
           e = 0.45, K = 50, Ea = 0.2, Em = 0.67,Er = 0.67,
           pr = -0.25, pa1 = -0.25, pa2 = -0.25, pm = -0.25, #pr影响
           mu1 = 1, mu2 = 1, z1 = 1, z2 = 1) #决定N P body size变化幅度
state <- c(N = 20, P = 10)
times <- seq(0, 500, by = 0.1)

sim2 = as.data.frame(
  ode(y = state,
      times = times,
      func = predator_prey_T,
      parms = parms,
      rootfun = root_func)
)

library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)

predator_prey_TS <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    # grouth rate of the prey
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pr)
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pa1) * (z2 ^ pa2)  # Arrhenius公式
    # metabolic rate of the predator
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * (z2 ^ pm)

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - pa1 * P * a)
    dz2 <- mu2 * (1 / z2) * (e * N * pa2 * a - pm * m)


    list(c(dN, dP, dz1, dz2),
         Temp = Temp,
         r = r,
         a = a,
         m = m,
         dN = dN,
         dP = dP,
         dz1 = dz1,
         dz2 = dz2)

    # list(c(dN, dP, dz1, dz2))
  })
}
parms <- c(T0 =273.15,k = 8.618*10^(-5),n = 0.02,
           a0 = 0.5, m0 = 0.5, r0 = 0.5,
           e = 0.45, K = 50, Ea = 0.2, Em = 0.67,Er = 0.67,
           pr = -0.25, pa1 = -0.25, pa2 = -0.25, pm = -0.25, #pr影响
           mu1 = 3, mu2 = 3) #决定N P body size变化幅度
state <- c(N = 20, P = 10, z1 = 1, z2 = 1)
times <- seq(0, 500, by = 0.1)



sim3 = as.data.frame(
  ode(y = state, times = times, func = predator_prey_TS, parms = parms)
)
sim1 <- sim1[,1:5]
sim2 <- sim2[,1:5]
sim3 <- sim3[,1:5]

# sim2$time = sim2$time * 10
sim1$run <- "TI"
sim2$run <- "TD"
sim3$run <- "TSD"
all_sims <- rbind(sim1, sim2,sim3)
df_long <- pivot_longer(all_sims,
                        cols = c("N", "P","z1","z2"),
                        names_to  = "species",
                        values_to = "biomass")


p <- ggplot(df_long, aes(x = time, y = biomass,
                         color   = species,
                         linetype = run)) +
  # scale_alpha_manual(values = c("PP" = 0.3,
  #                               "TD" = 0.6,
  #                               "TSD" = 1.0)) +
  geom_line(size = 0.8) +
  labs(x = "Time",
       y = "Biomass",
       color = "Species",
       linetype = "Model") +
  theme_minimal(base_size = 14) +
  ggplot2::scale_colour_manual("Species",values = c("orange","#4daf4a","#E90F44","#0B34E7","#FFC839","#984ea3"))+ #"red2","blue2",
  theme(legend.position = "bottom")

print(p)
