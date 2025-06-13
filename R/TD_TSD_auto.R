## plot TD VS TSD Predator-Prey model

library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)
library(gridExtra)

library(manipulate)
manipulate(
  plot(1:x, (1:x)^2, type="b"),
  x = slider(5, 50, initial=20, step=1)
)
predator_prey_TS <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # update temperature
    Temp <- T0 + n * t
    # update temperature-dependent rates
    # grouth rate of the prey
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pr)
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pa1) * (z2 ^ pa2)
    # metabolic rate of the predator
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * (z2 ^ pm)

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - pa1 * P * a)
    dz2 <- mu2 * (1 / z2) * (e * N * pa2 * a - pm * m)
    list(c(dN, dP, dz1, dz2),
         Temp = Temp)
  })
}

predator_prey_T <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pr)
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * (z1 ^ pa1) * (z2 ^ pa2) # Arrhenius公式
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * (z2 ^ pm)

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    list(c(dN, dP),
         # z1 = z1,
         # z2 = z2,
         Temp = Temp)
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
manipulate({
  parms <- list(
    n    = n,
    mu1  = mu1,
    mu2  = mu2,
    z1   = z1,
    z2   = z2,
    k    = 8.618*10^(-5),
    r0   = r0,
    a0   = a0,
    m0   = m0,
    Er   = Er,
    Ea   = Ea,
    Em   = Em,
    pr   = pr,
    pa1  = pa1,
    pa2  = pa2,
    pm   = pm,
    K    = K,
    e    = e,
    T0   = 270
  )
  # parms <- list(
  #   n    = 0.01,
  #   mu1  = 1,
  #   mu2  = 1,
  #   z1   = 2,
  #   z2   = 20,
  #   k    = 8.618*10^(-5),
  #   r0   = 0.4,
  #   a0   = 0.3,
  #   m0   = 0.4,
  #   Er   = 0.67,
  #   Ea   = 0.67,
  #   Em   = 0.67,
  #   pr   = -0.25,
  #   pa1  = -0.25,
  #   pa2  = -0.25,
  #   pm   = -0.25,
  #   K    = 100,
  #   e    = 0.45,
  #   T0   = 270
  # )
  times <- seq(0, 10/parms$n, by = 1)  # 10/parms$n
  out1 <- ode(c(N = 50, P = 20, z1 = parms$z1, z2 = parms$z2), times, predator_prey_TS, parms,
              rootfun = root_func,events  = list(func = event_func,root = TRUE))
  out2 <- ode(c(N = 50, P = 20), times, predator_prey_T, parms,
              rootfun = root_func, events  = list(func = event_func,root = TRUE))
  colnames(out2) <- c("time", "NT","PT","Temp")

  df  <- data.frame(out1,out2[,c(2,3)])

  ggplot(df, aes(x = Temp)) +
    ylim(0,120)+
    geom_line(aes(y = NT, colour = "TD Prey"),size = 1) +
    geom_line(aes(y = PT, colour = "TD Predator"),size = 1) +
    geom_line(aes(y = N, colour = "TSD Prey"),size = 1) +
    geom_line(aes(y = P, colour = "TSD Predator"),size = 1) +
    geom_line(aes(y = (z1), colour = "Size Prey"),size = 1) +
    geom_line(aes(y = (z2), colour = "Size Predator"),size = 1) +
    geom_line(aes(y = (z2/z1), colour = "Ratio"),size = 1) +
    scale_colour_manual(values = c("TSD Prey" = "red",
                                   "TSD Predator" = "blue",
                                   "Size Prey" = "gold",
                                   "Size Predator" = "green",
                                   "TD Prey" = "pink",
                                   "TD Predator" = "lightblue",
                                   "Ratio" = "black")) +
    labs(x = "Temperature",y = "Density", colour = "") +
    theme_minimal()



  # p2 <- ggplot(df, aes(x = N, y = P)) +
  #   geom_path(color = "red", size = 1.2) +
  #   labs(x = "Prey (N)", y = "Predator (P)") +
  #   theme_minimal() +
  #   ggtitle("Phase Plot: N vs P")
  #
  # grid.arrange(p1, p2, nrow = 1)
},

n  = slider(0, 0.1, initial = 0.01, step = 0.001),
mu1  = slider(0, 50, initial = 1, step = 0.01),
mu2  = slider(0, 500, initial = 1, step = 0.1),
z1     = slider(0.001, 100, initial = 1, step = 0.001),
z2     = slider(0.01, 1000, initial = 100, step = 0.01),
r0    = slider(0.01, 1, initial = 0.4, step = 0.01),
a0    = slider(0.01, 1, initial = 0.1, step = 0.01),
m0    = slider(0.01, 1, initial = 0.2, step = 0.01),
Er  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
Ea  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
Em  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
pr   = slider(-1, 0.9, initial = -0.25, step = 0.005),
pa1  = slider(-1, 0.9, initial = -0.25, step = 0.005),
pa2  = slider(-1, 0.9, initial = -0.25, step = 0.005),
pm   = slider(-1, 0.9, initial = -0.25, step = 0.005),
K     = slider(50, 150, initial = 100, step = 1),
e     = slider(0.1, 0.5, initial = 0.45, step = 0.01)
)

