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
    r <- r0 * exp(-Er / (k * Temp)) * z1 ^ pr
    # attack or feeding rate, z1 and z2 are body sizes of prey and predator
    a <- a0 * exp(-Ea / (k * Temp)) * z1 ^ pa1 * z2 ^ pa1  # Arrhenius公式
    # metabolic rate of the predator
    m <- m0 * exp(-Em / (k * Temp)) * z2 ^ pm

    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    dz1 <- mu1 * (1 / z1) * (pr * r * (1 - N / K) - pa1 * P * a)
    dz2 <- mu2 * (1 / z2) * (e * N * pa2 * a - pm * m)


    # list(c(dN, dP, dz1, dz2),
    #      Temp = Temp,
    #      r = r,
    #      a = a,
    #      m = m)

    list(c(dN, dP, dz1, dz2))
  })
}
parameters <- c(T0 =300,k = 8.617*10^(-5),n = 0.01,
                a0 = 10^7, m0 = 0.5*10^7, r0 = 10^7,
                  e = 0.5, K = 20, Ea = 0.8, Em = 0.8,Er = 0.8,
                pr = 0.3, pm = 0.3, pa1 = 0.3, pa2 = 0.3,
                mu1 = 0.5, mu2 = 0.5)
state <- c(N = 10, P = 10, z1 = 10, z2 = 20)
times <- seq(0, 20000, by = 1)

pp_results = as.data.frame(
  ode(y = state, times = times, func = predator_prey_TS, parms = parameters)
)


pp_results_time = gather(pp_results, Species, Biomass, -time)

ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ggplot2::scale_colour_manual("Species",values = c("#FFC839","#4daf4a","#E90F44","#63ADEE","#984ea3","black","grey"))

# plot N-P relationships over time
ggthemr('dust')
ggplot2::ggplot(pp_results) +
  geom_path(aes(N, P, colour = time), size = 1)+
  scale_color_gradientn(colours = rainbow(6))


# plot all the other parameters
predator_prey_TS <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
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







