#' temperature-size dependent predator-prey model

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
    r <- r0 * exp(Er *(Temp - T0)/ (k * T0 * Temp)) * z1 ^ pr
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(Ea *(Temp - T0)/ (k * T0 * Temp)) * z1 ^ pa1 * z2 ^ pa2  # Arrhenius公式
    # metabolic rate of the predator
    m <- m0 * exp(Em *(Temp - T0)/ (k * T0 * Temp)) * z2 ^ pm

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - pa1 * P * a)
    dz2 <- mu2 * (1 / z2) * (e * N * pa2 * a - pm * m)


    list(c(dN, dP, dz1, dz2),
         Temp = Temp,
         r = r,
         a = a,
         m = m,
         dz1 = dz1,
         dz2 = dz2)

    # list(c(dN, dP, dz1, dz2))
  })
}
parms <- c(T0 =300,k = 8.618*10^(-5),n = 0.001,
                a0 = 0.5, m0 = 0.3, r0 = 0.6,
                e = 0.45, K = 100, Ea = 0.6, Em = 0.6,Er = 0.6,
                pr = -0.25, pm = -0.25, pa1 = -0.25, pa2 = -0.25,
                mu1 = 1, mu2 = 1)
state <- c(N = 20, P = 20, z1 = 2, z2 = 20)
times <- seq(0, 2000, by = 1)



pp_results = as.data.frame(
  ode(y = state, times = times, func = predator_prey_TS, parms = parms)
)

pp_results1 <- pp_results[,1:5]
colnames(pp_results1) <- c("time","BN","BP","SN","SP")
pp_results_time = gather(pp_results1, Species, Biomass, -time)

pp_results2 <- pp_results[,c(1,10,11)]
pp_rates_time = gather(pp_results2, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ylim(0,30)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))

# ggthemr('dust')
ggplot(data = pp_rates_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))

# plot N-P relationships over time
ggthemr('dust')
ggplot2::ggplot(pp_results) +
  geom_path(aes(N, P, colour = time), size = 1)+
  scale_color_gradientn(colours = rainbow(6))


# plot all the other parms
predator_prey_TS <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # update temperature
    Temp <- T0 + n * t
    # temperature-dependent rates
    # grouth rate of the prey
    r <- r0 * exp(-Er / (k * Temp)) * z1 ^ pr
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(-Ea / (k * Temp)) * z1 ^ pa1 * z2 ^ pa1  # Arrhenius公式
    # metabolic rate of the predator
    m <- m0 * exp(-Em / (k * Temp)) * z2 ^ pm

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - pa1 * P * a)
    dz2 <- mu2 * (1 / z2) * (e * N * pa2 * a - pm * m)


    list(c(dN, dP, dz1, dz2),
         Temp = Temp,
         r = r,
         a = a,
         m = m)

    # list(c(dN, dP, dz1, dz2))
  })
}
pp_results1 <- pp_results[,1:5]
pp_results_time = gather(pp_results1, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))

pp_results2 <- pp_results[,c(1,7:9)]
pp_results_time = gather(pp_results2, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))







