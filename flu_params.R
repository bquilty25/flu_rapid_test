
influenza_gen_time <- epidist_db(
  disease = "influenza",
  epi_dist = "generation_time",
  single_epidist = TRUE
)

plot(influenza_gen_time)

influenza_gen_time$prob_dist

#sample from weihbull distribution
rweibull(10000, shape = 2.4, scale = 3.2) %>% hist()


influenza_incubation <- epidist_db(
  disease = "influenza",
  epi_dist = "incubation_period",
  single_epidist = TRUE
)

influenza_incubation$prob_dist
