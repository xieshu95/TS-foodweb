library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)

library(manipulate)
manipulate(
  plot(1:x, (1:x)^2, type="b"),
  x = slider(5, 50, initial=20, step=1)
)
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
manipulate({
  parms <- list(
    # K    = 100,
    # k    = 8.618*10^(-5),
    # a0   = 0.3,
    # Ea   = 0.67,
    # m0   = 0.4,
    # Em   = 0.67,
    # r0   = 0.4,
    # Er   = 0.67,
    # e    = 0.45,
    # T0   = 270,
    # n    = 0.01
    n    = n,
    K    = K,
    k    = 8.618*10^(-5),
    r0   = r0,
    a0   = a0,
    m0   = m0,
    Er   = Er,
    Ea   = Ea,
    Em   = Em,
    e    = e,
    T0   = 270
  )
  times <- seq(0, 20/parms$n, by = 1)
  out <- ode(c(N = 50, P = 20), times, predator_prey_T, parms)
  df  <- as.data.frame(out)

  ggplot(df, aes(x = Temp)) +
    geom_line(aes(y = N, colour = "Prey"),size = 1) +
    geom_line(aes(y = P, colour = "Predator"),size = 1) +
    scale_colour_manual(values = c("Prey" = "red", "Predator" = "blue")) +
    labs(y = "Density", colour = "") +
    theme_minimal()
},
n  = slider(0.001, 0.1, initial = 0.01, step = 0.002),
K     = slider(50, 150, initial = 100, step = 1),
r0    = slider(0.05, 1, initial = 0.4, step = 0.05),
a0    = slider(0.1, 1, initial = 0.3, step = 0.05),
m0    = slider(0.05, 1, initial = 0.4, step = 0.05),
Er  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
Ea  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
Em  = slider(-0.5, 0.9, initial = 0.67, step = 0.005),
e     = slider(0.1, 0.5, initial = 0.45, step = 0.01)
)



predator_prey_T <- function(t, state, parms) {
  N <- state[1]; P <- state[2]
  # update temperature
  Temp <- parms$T0 + parms$n * t
  # temperature-dependent rates
  a <- parms$a0 * exp(parms$Ea *(Temp - parms$T0)/ (parms$k * parms$T0 * Temp))  # Arrhenius公式
  m <- parms$m0 * exp(parms$Em *(Temp - parms$T0)/ (parms$k * parms$T0 * Temp))
  r <- parms$r0 * exp(parms$Er *(Temp - parms$T0)/ (parms$k * parms$T0 * Temp))

  dN <-  (1 - N/parms$K)*r*N - a*N*P
  dP <-  parms$e*a*N*P - m*P

  list(c(dN, dP))
}
# parms <- list(
#   K    = 100,
#   k    = 8.618*10^(-5),
#   a0   = 0.3,
#   Ea   = 0.67,
#   m0   = 0.4,
#   Em   = 0.67,
#   r0   = 0.4,
#   Er   = 0.67,
#   e    = 0.45,
#   T0   = 270,
#   n    = 0.01)
#
# state0 <- c(N = 50, P = 20)
#
# times <- seq(0, 100, by = 0.5)
# out <- ode(y = state0, times = times, func = predator_prey_T, parms = parms)
# df  <- as.data.frame(out)
# p1 <- ggplot(df, aes(x = time)) +
#   geom_line(aes(y = N, colour = "Prey")) +
#   geom_line(aes(y = P, colour = "Predator")) +
#   labs(y = "Density", colour = "") +
#   theme_minimal() +
#   transition_reveal(time)
