## auto parameetr setting for TSD Predator-Prey model

library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)

library(manipulate)
manipulate(
  plot(1:x, (1:x)^2, type="b"),
  x = slider(5, 50, initial=20, step=1)
)
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
         Temp = Temp)
         # r = r,
         # a = a,
         # m = m,
         # dN = dN,
         # dP = dP,
         # dz1 = dz1,
         # dz2 = dz2)

    # list(c(dN, dP, dz1, dz2))
  })
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
  times <- seq(0, 20/parms$n, by = 1)
  out <- ode(c(N = 50, P = 20, z1 = parms$z1, z2 = parms$z2), times, predator_prey_TS, parms)
  df  <- as.data.frame(out)

  ggplot(df, aes(x = Temp)) +
    geom_line(aes(y = N, colour = "Prey"),size = 1) +
    geom_line(aes(y = P, colour = "Predator"),size = 1) +
    geom_line(aes(y = log10(z1), colour = "Size Prey"),size = 1) +
    geom_line(aes(y = log10(z2), colour = "Size Predator"),size = 1) +
    scale_colour_manual(values = c("Prey" = "red", "Predator" = "blue", "Size Prey" = "gold", "Size Predator" = "green")) +
    labs(y = "Density", colour = "") +
    theme_minimal()
},

n  = slider(0.001, 0.1, initial = 0.01, step = 0.002),
mu1  = slider(0, 5, initial = 0.5, step = 0.005),
mu2  = slider(0, 5, initial = 0.5, step = 0.005),
z1     = slider(0, 100, initial = 1, step = 1),
z2     = slider(0, 10000, initial = 100, step = 10),
r0    = slider(0.05, 1, initial = 0.4, step = 0.05),
a0    = slider(0.1, 1, initial = 0.3, step = 0.05),
m0    = slider(0.05, 1, initial = 0.4, step = 0.05),
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

