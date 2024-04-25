pacman::p_load(tidyverse,tidytable,epiparameter,data.table,ggdist,qs,ggh4x)

source("flu_functions.R")
source("titre_dat_processing.R")

# Simulate index cases
set.seed(123)

t_interval = 0.2

index_cases <- generate_index_cases(num_cases = 100)

infections <- index_cases %>% 
  generate_system_delays() %>% 
  mutate(sec_cases = map(flu_generation_time, 
                         ~ simulate_empirical_traj(flu_generation_time = .x,
                                                   kinetics_dat = h1n1_titres,
                                                   num_sec_cases = 100)))

# Interpolate trajectories
infections_m <- infections %>% 
  unnest(sec_cases) %>%
  select(index,j,m) %>% 
  mutate(x = seq(0, 20, by = t_interval)) %>%
  mutate(res = map2(m, x, ~interpolate_trajectories(m = .x, x = .y))) %>% 
  select(index, j, res)

parms <- infections %>% 
  unnest(sec_cases) %>%
  select(-m) %>%
  crossing(
       bind_rows(
       crossing(scenario = "pcr",
               symptom_isolation = c(TRUE,FALSE)),
       crossing(scenario = c("quarantine_fast", "quarantine_slow"),
                                   symptom_isolation = c(TRUE,FALSE)),
       crossing(scenario = "rapid", 
                                   n_tests = c(1,3,5),  
                                   symptom_isolation = c(TRUE,FALSE)))
)

res <- infections_m %>% 
  unnest(res) %>% 
  left_join(parms)

rm(infections_m,parms)
gc()

#safe_simulate_isolation <- safely(simulate_isolation)

#plot sec cases
# infections_m %>% 
#   unnest(res) %>% 
#   filter(j==1) %>% 
#   ggplot(aes(x = x, y = y, group=j)) +
#   geom_line(alpha=0.1) +
#   labs(x = "Days since symptom onset", y = "Viral titre", title = "Simulated Viral Titre Trajectories",
#        color = "Symptom Isolation") +
#   theme_minimal()

#if res2 is not in environment, load from file. if not in file, run simulation
if(!exists("res2")){
  
  if(file.exists("res2.qs")){
    
    res2 <- qread("res2.qs")
    
  } else {
    
    #simulation isolation
    res2 <- res %>% 
      group_split(scenario,index,symptom_isolation,n_tests,.keep = T) %>% 
      map(~simulate_isolation(.x, scenario = .x$scenario, symptom_isolation = .x$symptom_isolation, n_tests = .x$n_tests),
          .progress = list(
            type = "iterator", 
            format = "Calculating {cli::pb_bar} {cli::pb_percent}",
            clear = TRUE))
    
    qsave(res2, file = "res2.qs")
    
  }
}

res2 %>% 
  bind_rows(.id = "index") %>% 
  group_by(index, scenario, symptom_isolation, n_tests) %>% 
  summarise(inf_prevented = sum(area_test_symptom)/sum(area_all)) %>% 
  #mutate(n_tests_lab=ifelse(is.na(n_tests),"PCR testing",paste0("Days of rapid testing:\n",n_tests))) %>% 
  #filter(symptom_isolation == FALSE) %>%
  ggplot(aes(x = fct_recode(scenario, "PCR and formal tracing" = "pcr", 
                            "Rapid testing and rapid tracing followed by daily testing" = "rapid", 
                            "Quarantine of contacts following rapid informal tracing" = "quarantine_fast", 
                            "Quarantine of contacts following formal tracing delays" = "quarantine_slow")
             , y = inf_prevented, group = n_tests, colour=scenario, shape=factor(n_tests))) +
  stat_pointinterval(position = position_dodge(width = 0.5),point_size = 3)+
  # geom_violin(alpha=0.7) +
  # geom_boxplot(width = 0.1, fill=NA, outliers=F) +
  facet_grid(~symptom_isolation, 
              axes="all", 
              #remove_labels = "y",
              scales = "free",
              labeller = labeller(symptom_isolation = function(x) paste("Symptomatic self-isolation:", tolower(x))))+
  scale_colour_brewer(name = "Scenario", palette = "Set2",guide=F)+
  scale_shape_manual(name = "Number of days of testing", values=c(15:19),na.value = 19)+
  lims(y = c(0, 1)) +
  labs(x = "", y = "Transmission potential averted", title = "Effectiveness of rapid testing and informal tracing compared to a centralised PCR testing and formal tracing system in reducing influenza transmission",
       subtitle = "\nAssumed delays:\n- Rapid testing: contacts traced informally immediately upon index's symptom onset.\n- PCR testing: mean 1 day for index case to get test after symptom onset + mean 2 days for results + mean 2 days for contact tracing.\n\nAssumed test sensitivity:\n- Rapid testing: 50% sensitive at 2 log titre.\n- PCR testing: 50% sensitive at 1 log titre.\n\nOther assumptions:\n- 66.9% of cases are symptomatic.\n- 100% compliance with testing, tracing, isolation upon test result or symptoms.",
       fill = "Symptom Isolation") +
  plotting_theme+
  coord_flip()

ggsave("flu_rapid_testing_versus_pcr.png", width = 297, height = 210, dpi = 600, units = "mm", bg= "white")

res2 %>%
  mutate(n_tests_lab=ifelse(is.na(n_tests),"PCR testing",paste0("5 days of rapid testing: ",n_tests))) %>%
  #filter(symptom_isolation == FALSE) %>%
  ggplot(aes(x = n_tests_lab, y = inf_prevented, fill = factor(symptom_isolation))) +
  stat_dotsinterval()+
  facet_grid(~symptom_isolation,labeller = label_both)+
  scale_fill_brewer(palette = "Set2")+
  lims(y = c(0, 1)) +
  labs(x = "Daily rapid tests and instant tracing", y = "Transmission potential averted", title = "Effect of Testing, Isolation, and Symptom Onset on Infection Prevention",
       fill = "Symptom Isolation") +
  theme_minimal()+
  ggpubr::rotate_x_text(angle = 45)

