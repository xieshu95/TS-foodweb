library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)
#' predator-prey ODE with T dependence but without size dependence(type I)
#' N is the biomass of prey and P is the biomass of predator
#'
predator_prey_T <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    a <- a0 * exp(Er *(Temp - T0)/ (k * T0 * Temp))  # Arrhenius公式
    m <- m0 * exp(Er *(Temp - T0)/ (k * T0 * Temp))
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp))

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    list(c(dN, dP))
         # Temp = Temp,
         # dN  = dN,
         # dP = dP,
         # r = r,
         # a = a,
         # m = m)
  })
}
parameters <- c(T0 =300,k = 8.617*10^(-5),n = 0.001,
                a0 = 0.5, m0 = 0.3, r0 = 0.6,
                e = 0.45, K = 100, Ea = 0.6, Em = 0.6,Er = 0.6)
state <- c(N = 20, P = 20)
times <- seq(0, 2000, by = 1)

pp_results = as.data.frame(
  ode(y = state, times = times, func = predator_prey_T, parms = parameters)
)
# pp_results$Temp <- pp_results$Temp - 300


pp_results_time = gather(pp_results, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ylim(0,30)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","orange"))

# plot N-P relationships over time
ggthemr('dust')
ggplot2::ggplot(pp_results) +
  geom_path(aes(N, P, colour = time), size = 1)+
  scale_color_gradientn(colours = rainbow(6))

T0 =300
k = 8.617*10^(-5)
n = 0.01
a0 = 10^7
m0 = 10^6
r0 = 0.5*10^7
e = 0.85
K = 20
Ea = 0.6
Em = 0.5
Er = 0.6


