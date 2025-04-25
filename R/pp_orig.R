library(tidyverse)
library(deSolve)
library(ggplot2)
library(ggthemr)
#' original predator-prey ODE without T and size dependence(type I)
#' N is the biomass of prey and P is the biomass of predator
#'
predator_prey_type1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dN <-  (1 - N/K)*r*N - a*N*P
    dP <-  e*a*N*P - m*P
    list(c(dN, dP))
  })
}
parameters <- c(r = 0.3, a = 0.1, m = 0.4, e = 0.85, K = 20)
state <- c(N = 10, P = 10)
times <- seq(0, 300, by = 0.5)

pp_results = as.data.frame(
  ode(y = state, times = times, func = predator_prey_type1, parms = parameters)
)


pp_results_time = gather(pp_results, Species, Biomass, -time)

# plot biomass N and P dynamics over time
ggthemr('dust')
ggplot(data = pp_results_time,aes(time, Biomass, colour = Species)) +
  geom_line(size = 1)+
  ggplot2::scale_colour_manual("Species",values = c("#E90F44","#63ADEE"))

# plot N-P relationships over time
ggthemr('dust')
ggplot2::ggplot(pp_results) +
  geom_path(aes(N, P, colour = time), size = 1)+
  scale_color_gradientn(colours = rainbow(6))

  # scale_color_gradient(low = "yellow", high = "red", na.value = NA)



# lotka_volterra <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#     dN <-  r*N - a*N*P
#     dP <-  e*a*N*P - m*P
#     list(c(dN, dP))
#   })
# }
# parameters <- c(r = 0.3, a = 0.1, m = 0.4, e = 0.5)
# state <- c(N = 2, P = 2)
# times <- seq(0, 100, by = 0.2)
# sol = as.data.frame(
#   ode(y = state, times = times, func = lotka_volterra, parms = parameters)
# )
#
# sol2 = gather(sol, species, individuals, -time)
#
# ggplot(sol2) +
#   geom_line(aes(time, individuals, colour = species))+
#   theme_minimal() +
#   theme(
#     plot.background = element_rect(fill = rgb(.2,.21,.27)),
#     text = element_text(colour = 'grey'),
#     axis.text = element_text(colour = 'grey'),
#     panel.grid = element_line(colour = 'grey')
#   )
#
# ggplot(sol) +
#   geom_path(aes(N, P), colour = 'red', size = 1)+
#   theme_minimal() +
#   theme(
#     plot.background = element_rect(fill = rgb(.2,.21,.27)),
#     text = element_text(colour = 'grey'),
#     axis.text = element_text(colour = 'grey'),
#     panel.grid = element_line(colour = 'grey')
#   )
