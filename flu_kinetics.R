pacman::p_load(tidyverse,tidytable,epiparameter,data.table,ggdist,qs,ggh4x)

source("flu_functions.R")
source("titre_dat_processing.R")

# Simulate index cases
set.seed(123)

t_interval = 0.2

index_cases <- generate_index_cases(num_cases = 1000)

infections <- index_cases %>% 
  generate_system_delays() %>% 
  mutate(sec_cases = map(flu_generation_time, 
                         ~ simulate_empirical_traj(flu_generation_time = .x,
                                                   kinetics_dat = h1n1_titres,
                                                   num_sec_cases = 100)))

# Interpolate trajectories
infections_m <- infections %>% 
  unnest(sec_cases) %>%
  select(index,j,m,flu_generation_time) %>% 
  mutate(x = seq(0, 20, by = t_interval)) %>%
  mutate(res = map2(m, x, ~interpolate_trajectories(m = .x, x = .y))) %>% 
  select(index, j, res,flu_generation_time)

parms <- infections %>% 
  unnest(sec_cases) %>%
  select(-m) %>%
  crossing(
       bind_rows(
      crossing(scenario = "pcr",
                n_tests = NA,
               symptom_isolation = c(TRUE,FALSE)),
      crossing(scenario = c("quarantine_fast", "quarantine_slow"),
               n_tests = NA,
               symptom_isolation = c(TRUE,FALSE)),
      crossing(scenario = "rapid",
               n_tests = c(3,5,7),
               symptom_isolation = c(TRUE,FALSE)),
       crossing(scenario = "no_test_iso",
                n_tests = NA,
                symptom_isolation = c(TRUE,FALSE)
       )
       )
  )


res <- infections_m %>% 
  unnest(res) %>% 
  left_join(parms)

 #simulation isolation
 res2 <- res %>% 
      group_split(scenario,index,symptom_isolation,n_tests,.keep = T) %>% 
      map(~simulate_isolation(.x, scenario = .x$scenario, 
                              symptom_isolation = .x$symptom_isolation, 
                              n_tests = .x$n_tests))
    
  qsave(res2, file = "res2.qs")


res2 %>% 
  bind_rows(.id = "index") %>% 
  group_by(index, scenario, symptom_isolation, n_tests) %>% 
  summarise(inf_prevented = (sum(area_all)-sum(area_infectious))/sum(area_all)) %>% 
  group_by(scenario, symptom_isolation, n_tests) %>% 
  summarise(median_inf_prevented = median(inf_prevented),
            lower_inf_prevented = quantile(inf_prevented, 0.025),
            upper_inf_prevented = quantile(inf_prevented, 0.975))

res2 %>% 
  bind_rows(.id = "index") %>% 
  group_by(index, scenario, symptom_isolation, n_tests) %>% 
  summarise(inf_prevented = (sum(area_all)-sum(area_infectious))/sum(area_all)) %>% 
  #mutate(n_tests_lab=ifelse(is.na(n_tests),"PCR testing",paste0("Days of rapid testing:\n",n_tests))) %>% 
  mutate(scenario = fct_recode(scenario, "PCR and formal tracing" = "pcr", 
                    "Rapid testing and rapid tracing followed by daily testing" = "rapid", 
                    "Quarantine of contacts following rapid informal tracing" = "quarantine_fast", 
                    "Quarantine of contacts following formal tracing with delays" = "quarantine_slow")) %>% 
  #filter(symptom_isolation == FALSE) %>%
  ggplot(aes(x = str_wrap(scenario,20),
             y = inf_prevented, group = n_tests, colour=scenario, shape=factor(n_tests))) +
  stat_pointinterval(position = position_dodge(width = 0.5),point_size = 3,)+
  # geom_violin(alpha=0.7) +
  # geom_boxplot(width = 0.1, fill=NA, outliers=F) +
  facet_grid(~symptom_isolation, 
              axes="all", 
              #remove_labels = "y",
              scales = "free",
              labeller = labeller(symptom_isolation = function(x) paste("Symptomatic self-isolation:", tolower(x))))+
  scale_colour_brewer(name = "Scenario", palette = "Set2",guide=F)+
  scale_shape_manual(name = "Number of days of rapid testing", values=c(15:19),na.value = 1)+
  lims(y = c(0, 1)) +
  labs(x = "", y = "Transmission potential averted", 
       #title = str_wrap("Effectiveness of rapid testing and informal tracing compared to a centralised PCR testing and formal tracing system in reducing influenza transmission",60),
 #      subtitle = "\nAssumed delays:\n- Rapid testing: contacts traced informally immediately upon index's symptom onset.\n- PCR testing: mean 1 day for index case to get test after symptom onset + mean 2 days for results + mean 2 days for contact tracing.\n\nAssumed test sensitivity:\n- Rapid testing: 50% sensitive at 2 log titre.\n- PCR testing: 50% sensitive at 1 log titre.\n\nOther assumptions:\n- 66.9% of cases are symptomatic.\n- 100% compliance with testing, tracing, isolation upon test result or symptoms.",
       fill = "Symptom Isolation") +
  plotting_theme+
  coord_flip()

ggsave("flu_rapid_testing_versus_pcr.png", width = 210, height = 150, dpi = 600, units = "mm", bg= "white")
