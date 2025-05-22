# pp_TS衍生，test parameters

library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)

predator_prey_TS <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    # grouth rate of the prey
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * z1 ^ pr
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * z1 ^ pa1 * z2 ^ pa2  # Arrhenius公式
    # handling time, z1 and z2 are body sizes of prey and predator
    h <- h0 * exp(Eh *(Temp - T0)/ (k * T0 * Temp)) * z1 ^ ph1 * z2 ^ ph2
    # metabolic rate of the predator
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * z2 ^ pm

    dN <-  (1 - N/K)*r*N - a*N*P/(1+h*a*N)
    dP <-  e*a*N*P/(1+h*a*N) - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - P * a * (pa1-ph1*h*a*N)/((1+h*a*N)^2))
    dz2 <- mu2 * (1 / z2) * (e * N * a * (pa2-ph2*h*a*N)/((1+h*a*N)^2) - pm * m)


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
root_func <- function(t, state, parameters) {
  with(as.list(state), {
    # 返回两个根，任意一个为 0 时触发事件
    c(N, P,z1,z2)
  })
}
# parameters <- c(T0 =283.15,k = 8.618*10^(-5),n = 0.001,
#                 a0 = 0.5, m0 = 0.5, r0 = 0.5,
#                 e = 0.45, K = 50, Ea = 0.67, Em = 0.67,Er = 0.67,
#                 pr = -0.25, pa1 = -0.25, pa2 = 0.75, pm = -0.25, #pr影响
#                 mu1 = 0.5, mu2 = 0.5) #决定N P body size变化幅度
# Binzer 2015
parameters <- c(T0 =283.15,k = 8.618*10^(-5),n = 1,
                a0 = exp(-13.1), m0 = exp(-16.54), r0 = exp(-15.68),h0 = exp(9.66),
                e = 0.85, K = 20, Ea = -0.38, Em = -0.69,Er = -0.84,Eh = 0.26,
                pr = -0.25, pa1 = 0.25, pa2 = -0.8, pm = -0.31, ph1 = -0.45, ph2 = 0.47,
                mu1 = 1, mu2 = 1)
state <- c(N = 10, P = 5, z1 = 10, z2 = 100)
times <- seq(0, 10, by = 0.1)



pp_results = as.data.frame(
  ode(y = state, times = times, func = predator_prey_TS, parms = parameters)
)

pp_results1 <- pp_results[,1:5]
colnames(pp_results1) <- c("time","BN","BP","SN","SP")
pp_results_time = gather(pp_results1, Species, Biomass, -time)

pp_results2 <- pp_results[,c(1,10,11)]
pp_rates_time = gather(pp_results2, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  # ylim(0,30)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))
