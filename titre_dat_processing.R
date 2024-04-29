#Carrat F, Vergu E, Ferguson NM, Lemaitre M, Cauchemez S, Leach S, Valleron AJ. Time lines of infection and disease in human influenza: a review of volunteer challenge studies. Am J Epidemiol. 2008 Apr 1;167(7):775-85. doi: 10.1093/aje/kwm375. Epub 2008 Jan 29. PMID: 18230677.

h1n1_titres <- read_csv("h1n1_titres_over_time.csv",col_names = c("day","log_titre"))

h1n1_titres <- h1n1_titres %>% 
  mutate(day = round(day,1)) %>% 
  group_by(day) %>% 
  summarise(mean_log_titre = median(log_titre),
            upper_log_titre = max(log_titre),
            lower_log_titre = min(log_titre)) %>%
  mutate(se=(upper_log_titre-lower_log_titre),
         sd=se/sqrt(n())) 


h1n1_titres %>%
  crossing(i=1:100) %>% 
  mutate(trajectory = rnorm(n(),mean = mean_log_titre,sd = sd)) %>%
  ggplot(aes(x=day,y=mean_log_titre,ymin = lower_log_titre,ymax = upper_log_titre))+
  geom_pointrange()+
  geom_line(aes(y=trajectory,group=i),alpha=0.1)+
  labs(x = "Days since exposure", y = "Log viral titre",
       color = "Symptom Isolation") +
  plotting_theme

ggsave("sec_cases_trajectories.png", width = 210, height = 150, dpi = 600, units = "mm", bg= "white")

h3n2_titres <- read_csv("h3n2_titres_over_time.csv",col_names = c("day","log_titre"))

h3n2_titres <- h3n2_titres %>% 
  mutate(day = round(day,1)) %>% 
  group_by(day) %>% 
  summarise(mean_log_titre = median(log_titre),
            upper_log_titre = max(log_titre),
            lower_log_titre = min(log_titre)) %>%
  mutate(se=(upper_log_titre-lower_log_titre)/2,
         sd=se/sqrt(n())) 


h3n2_titres %>%
  crossing(i=1:100) %>% 
  mutate(trajectory = rnorm(n(),mean = mean_log_titre,sd = sd)) %>%
  ggplot(aes(x=day,y=mean_log_titre,ymin = lower_log_titre,ymax = upper_log_titre))+
  geom_pointrange()+
  geom_line(aes(y=trajectory,group=i),alpha=0.1)
