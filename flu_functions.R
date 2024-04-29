# Define the logistic function
sensitivity_function <- function(x, rapid = F) {
  if(rapid){
    p <- 1 / (1 + exp(-4 * (x - 2)))
  } else {
    p <- 1 / (1 + exp(-4 * (x - 1)))
  }
  return(p)
}

#plot sensitivity function with ggplot2
x <- seq(0, 4, 0.1)
y_rapid <- sensitivity_function(x, rapid = T)
y_pcr <- sensitivity_function(x, rapid = F)

data.frame(x = x, 
           y_rapid = y_rapid,
           y_pcr = y_pcr) %>%
  pivot_longer(-x, names_to = "test_type", values_to = "y") %>% 
  ggplot(aes(x = x, y = y, colour=test_type)) +
  geom_line() +
  scale_color_brewer(palette = "Set1", labels=c("PCR", "Rapid test")) +
  labs(x = "Log viral titre", y = "Probability of positive test",
       color = "Symptom Isolation") +
  plotting_theme

ggsave("sensitivity_curves.png", width = 210, height = 150, dpi = 600, units = "mm", bg= "white")


# Generate index cases
generate_index_cases <- function(num_cases) {
  
  flu_generation_time <- rweibull(n=num_cases, shape = 2.4, scale = 3.2)
  
  index_symptom_onset <- rweibull(n=num_cases, shape = 2.1, scale = 3.8)
  
  data.frame(
    index = 1:num_cases,
    index_symptom_onset = index_symptom_onset,
    flu_generation_time = flu_generation_time
  )
  
}

generate_system_delays <- function(index_cases) {
  
  index_cases %>% 
    mutate(
      onset_to_test_delay = rnorm(n(), mean = 1, sd = 0.5),
      test_to_result_delay = rnorm(n(), mean = 2, sd = 0.5),
      test_to_tracing_delay = rnorm(n(), mean = 2, sd = 0.5)
    ) %>% 
    mutate(
      total_delay = index_symptom_onset + onset_to_test_delay + test_to_result_delay + test_to_tracing_delay
    )
}

# # Function to generate piecewise linear data
# generate_piecewise_linear <- function(num_points, params, flu_generation_time) {
#   
#   peak_x_mean <- params$peak_x_mean
#   peak_x_sd <- params$peak_x_sd
#   peak_y_mean <- params$peak_y_mean
#   peak_y_sd <- params$peak_y_sd
#   max_x_mean <- params$max_x_mean
#   max_x_sd <- params$max_x_sd
#   
#   # Duration
#   duration <- rnorm(1, mean = max_x_mean, sd = max_x_sd)
#   
#   # Generate breakpoints
#   breakpoints <- c(0, rnorm(1, mean = peak_x_mean, sd = peak_x_sd), duration)
#   
#   # Generate values
#   values <- c(0, rnorm(1, mean = peak_y_mean, sd = peak_y_sd), 0)
#   
#   # Interpolate between breakpoints
#   piecewise_linear <- function(x) {
#     ifelse(x <= breakpoints[2],
#            values[1] + (values[2] - values[1]) / (breakpoints[2] - breakpoints[1]) * (x - breakpoints[1]),
#            values[2] + (values[3] - values[2]) / (breakpoints[3] - breakpoints[2]) * (x - breakpoints[2]))
#   }
#   
#   # Generate x values
#   x_values <- seq(0, duration, length.out = num_points)
#   
#   # Generate y values using the piecewise linear function
#   y_values <- piecewise_linear(x_values)
#   
#   # Add generation time to x values
#   x_values <- x_values + flu_generation_time
#   
#   return(data.frame(x = x_values, y = y_values))
# }
# 
# # Function to simulate viral load trajectories 
# simulate_trajectories <- function(num_sec_cases, num_points, params, flu_generation_time) {
#   
#   # Generate multiple trajectories
#   trajectories <- lapply(1:num_sec_cases, function(i) {
#     generate_piecewise_linear(num_points = num_points, params, flu_generation_time)
#   })
#   
#   # Combine all trajectories into one data frame
#   trajectories <- do.call(rbind, lapply(1:length(trajectories), function(i) {
#     cbind(trajectories[[i]], trajectory = i)
#   }))
#   
#   return(trajectories)
# }

# Function to simulate empirical trajectories from the kinetics data
simulate_empirical_traj <- function(kinetics_dat, num_sec_cases, flu_generation_time){
  
  #browser()
  
  # Create individual trajectories
  kinetics_dat %>% 
    crossing(j=1:num_sec_cases) %>%
    mutate(y_values = rnorm(n(), mean = mean_log_titre, sd = sd), 
           # Add generation time to x values
           x_values = day + flu_generation_time) %>% 
    select(j, x = x_values, y = y_values) %>% 
    group_by(j) %>% 
    nest() %>% 
    mutate(m = map(data,~approxfun(x=.x$x, y=.x$y))) %>% 
    select(j,m) %>% 
    ungroup() %>% 
    #chance of being symptomatic: 66.9% Carrat et al.
    mutate(symptomatic = rbinom(n=n(), 1, prob = 0.669),
           sec_symptom_onset = flu_generation_time + rweibull(n=n(), shape = 2.1, scale = 3.8))
  
}

# Function to interpolate viral load trajectories
interpolate_trajectories <- function(m, x) {
  
  #browser()
  
  y = m(x)
  
  data.frame(x=x, y = y) %>% 
    mutate(y = replace_na(y,0))

}

# Function to simulate testing, isolation due to symptoms, and isolation due to positive test
simulate_isolation <- function(df, scenario, symptom_isolation = FALSE, n_tests) {
  
  #browser()
  
  scenario <-  unique(scenario)
  symptom_isolation <- unique(symptom_isolation)
  index_symptom_onset <- unique(df$index_symptom_onset)
  n_tests <- unique(n_tests)
  
  # Integrate under the curve
  df_inf <- df %>%
    group_by(j) %>%
    summarise(area_all = sum(y))
  
  if(scenario == "rapid"){
    
    # Probability of testing positive
    df$test_prob <- sensitivity_function(df$y, rapid = TRUE)
    
    if(n_tests==1){
      test_days <- round(index_symptom_onset/t_interval)*t_interval
    } else {
      test_days <- seq(round(index_symptom_onset/t_interval)*t_interval, 
                       round(index_symptom_onset/t_interval)*t_interval + n_tests - 1, 
                       by = 1)
    }
    
    #Repeat daily testing until positive
    test_df <- df %>%
      arrange(j, x) %>%
      group_by(j) %>%
      rowwise() %>% 
      filter(x %in% test_days) %>%
      group_by(j) %>% 
      mutate(test_day_index = row_number()) %>%
      rowwise() %>% 
      mutate(result = rbinom(n(), 1, test_prob),
             isolation_start = ifelse(result == 1, x, Inf)) %>% 
      group_by(j) %>% 
      slice_min(isolation_start, with_ties = F) %>% 
      ungroup() %>% 
      mutate(test_day_index = ifelse(isolation_start == Inf, n_tests, test_day_index)) %>%
      select(j, n_tests_used = test_day_index, isolation_start) %>% 
      arrange(j)
    
  } else if (scenario == "pcr") {
    # Apply one test per trajectory on testing day
    
    df$test_prob <- sensitivity_function(df$y, rapid = FALSE)
    
    test_df <- df %>%
      group_by(j) %>%
      slice(which.min(x >= total_delay)) %>%
      rowwise() %>%
      mutate(result = rbinom(n(), 1, test_prob)) %>%
      mutate(isolation_start = ifelse(result == 1, x, Inf)) %>%
      mutate(n_tests_used = 1) %>% 
      select(j, n_tests_used, isolation_start) %>% 
      ungroup()
    
  } else if (scenario == "quarantine_fast") {
    
    # isolate upon index cases
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= index_symptom_onset) %>%
      slice(1) %>%
      mutate(n_tests_used = 0) %>% 
      select(j, n_tests_used, isolation_start = x) %>% 
      ungroup()
    
  } else if (scenario == "quarantine_slow") {
    
    # isolate upon index cases
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= total_delay) %>%
      slice(1) %>%
      mutate(n_tests_used = 0) %>% 
      select(j, n_tests_used, isolation_start = x) %>% 
      ungroup()
    
  } else if (scenario == "no_test_iso") {
    
    test_df <- df %>%
      group_by(j) %>% 
      summarise(n_tests_used = 0,
             isolation_start = Inf) %>% 
      select(j, n_tests_used, isolation_start) %>% 
      ungroup()
    
  }
  
  # Simulate isolation due to symptom onset if symptom_isolation is TRUE
  if (symptom_isolation) {
    
    #
    df_inf_test_symptom <- df %>%
      mutate(sec_symptom_onset = ifelse(symptomatic==1,
                                        sec_symptom_onset, 
                                        Inf)) %>%
      left_join(test_df) %>%
      group_by(j,n_tests_used) %>%
      mutate(in_iso = x >= pmin(isolation_start,sec_symptom_onset)) %>%
      group_by(j,n_tests_used,in_iso) %>%
      filter(in_iso == F) %>% 
      summarise(area_infectious = sum(y))
    
    # Calculate the prevention in infections due to testing, isolation, and symptom onset
    df_res <- df_inf %>%
      left_join(df_inf_test_symptom) %>%
      mutate(area_infectious = ifelse(is.na(area_infectious), 1, area_infectious)) %>% 
      mutate(inf_prevented = ((area_all-area_infectious) / area_all))
    
  } else {
    
    # Join df and test df
    df_inf_test <- df %>%
      right_join(test_df) %>%
      group_by(j,n_tests_used) %>%
      mutate(in_iso = x >= isolation_start) %>%
      group_by(j,n_tests_used,in_iso) %>%
      filter(in_iso == F) %>% 
      summarise(area_infectious = sum(y))
    
    # Calculate the prevention in infections due to testing and isolation
    df_res <- df_inf %>%
      left_join(df_inf_test) %>%
      mutate(area_infectious = ifelse(is.na(area_infectious), 1, area_infectious)) %>% 
      mutate(inf_prevented = ((area_all-area_infectious) / area_all))
    
  }
  
  df_res %>% mutate(scenario = scenario, 
                    symptom_isolation = symptom_isolation, 
                    n_tests = n_tests,
                    n_tests_used = n_tests_used)
  
}

plotting_theme <- theme_minimal()+
  theme(axis.ticks = element_line(),
        axis.title = element_text(),
        axis.text = element_text(),
        strip.text = element_text(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        #panel.border = element_rect(fill=NA,),
        #panel.grid = element_blank(),
        legend.position = "bottom",
        strip.placement = "outside",
        axis.line = element_line(),
        line = element_line(),
        text = element_text())
