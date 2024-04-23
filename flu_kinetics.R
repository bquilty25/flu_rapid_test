pacman::p_load(tidyverse,tidytable,epiparameter,data.table,ggdist)

source("flu_functions.R")

# Example usage
flu_params <- list(
  peak_x_mean = 2.5,
  peak_x_sd = 1,
  peak_y_mean = 1000,
  peak_y_sd = 20,
  max_x_mean = 10,
  max_x_sd = 2
)

# Simulate index cases
set.seed(123)

index_cases <- generate_index_cases(num_cases = 100)

infections <- index_cases %>% 
  generate_system_delays() %>% 
  mutate(sec_cases = map(flu_generation_time, ~ simulate_trajectories(flu_generation_time = .x, 
                                                                      params = flu_params,
                                                                      num_sec_cases = 10,
                                                                      num_points = 100)))

res <- infections %>% 
  unnest(sec_cases) %>%
  crossing(symptom_isolation = c(TRUE,FALSE),
           rapid= c(TRUE,FALSE))

res2 <- res[, simulate_isolation(.SD, symptom_isolation = symptom_isolation, rapid = rapid), by = c("index", "symptom_isolation","rapid")]

res2 %>% 
  ggplot(aes(x = rapid, y = 1 - inf_prevented, fill = factor(symptom_isolation))) +
  stat_dotsinterval() +
  facet_grid(~symptom_isolation,labeller = label_both)+
  lims(y = c(0, 1)) +
  labs(x = "Rapid testing and tracing", y = "Transmission potential averted", title = "Effect of Testing, Isolation, and Symptom Onset on Infection Prevention",
       fill = "Symptom Isolation") +
  theme_minimal()

