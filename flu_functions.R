# Define the logistic function
sensitivity_function <- function(x, rapid = F) {
  if(rapid){
    p <- 1 / (1 + exp(-4 * (x - 1)))
  } else {
    p <- 1 / (1 + exp(-4 * (x - 2)))
  }
  return(p)
}

#plot sensitivity function with ggplot2
x <- seq(0, 5, 0.1)
y_rapid <- sensitivity_function(x, rapid = T)
y_pcr <- sensitivity_function(x, rapid = F)

data.frame(x = x, 
           y_rapid = y_rapid,
           y_pcr = y_pcr) %>%
  pivot_longer(-x, names_to = "test_type", values_to = "y") %>% 
  ggplot(aes(x = x, y = y, colour=test_type)) +
  geom_line() +
  labs(x = "Viral titres", y = "Sensitivity") +
  theme_minimal()

# Generate index cases
generate_index_cases <- function(num_cases) {
  
  flu_generation_time <- rweibull(n=num_cases, shape = 2.4, scale = 3.2)
  
  symptom_onset <- rweibull(n=num_cases, shape = 2.1, scale = 3.8)
  
  data.frame(
    index = 1:num_cases,
    symptom_onset = symptom_onset,
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
      total_delay = symptom_onset + onset_to_test_delay + test_to_result_delay + test_to_tracing_delay
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
    ungroup()
  
}

# Function to interpolate viral load trajectories
interpolate_trajectories <- function(m, x) {
  
  #browser()
  
  y = m(x)
  
  data.frame(x=x, y = y) %>% 
    drop_na(y)

}

# Function to simulate testing, isolation due to symptoms, and isolation due to positive test
simulate_isolation <- function(df, scenario, symptom_isolation = FALSE, n_tests) {
  
  browser()
  
  scenario <-  unique(scenario)
  symptom_isolation <- unique(symptom_isolation)
  n_tests <- unique(n_tests)
  
  # Integrate under the curve
  df_inf <- df %>%
    group_by(j) %>%
    summarise(area_all = sum(y))
  
  # Find the peak viral load for each trajectory
  peak_x <- df %>%
    group_by(j) %>%
    filter(y == max(y)) %>% 
    select(j,peak_x = x)
  
  if(scenario == "rapid"){
    
    # Probability of testing positive
    df$test_prob <- sensitivity_function(df$y, rapid = TRUE)
    
    #Repeat daily testing until positive
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= symptom_onset) %>%
      group_by(j) %>%
      slice(1:n_tests) %>%
      mutate(result = rbinom(n(), 1, test_prob)) %>%
      mutate(isolation_start = ifelse(result == 1, x, NA)) %>% 
      group_by(j) %>%
      filter(result == 1) %>%
      arrange(j, x) %>% 
      slice(1) %>%
      select(j, isolation_start) %>% 
      ungroup()
    
  } else if (scenario == "pcr") {
    # Apply one test per trajectory on testing day
    
    df$test_prob <- sensitivity_function(df$y, rapid = FALSE)
    
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= total_delay) %>%
      slice(1) %>%
      rowwise() %>%
      mutate(result = rbinom(n(), 1, test_prob)) %>%
      mutate(isolation_start = ifelse(result == 1, x, NA)) %>%
      select(j, isolation_start) %>% 
      ungroup()
    
  } else if (scenario == "quarantine_fast") {
    
    # isolate upon index cases
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= symptom_onset) %>%
      slice(1) %>%
      select(j, isolation_start = x) %>% 
      ungroup()
    
  } else if (scenario == "quarantine_slow") {
    
    # isolate upon index cases
    test_df <- df %>%
      group_by(j) %>%
      filter(x >= total_delay) %>%
      slice(1) %>%
      select(j, isolation_start = x) %>% 
      ungroup()
    
  }
  
  # Simulate isolation due to symptom onset if symptom_isolation is TRUE
  if (symptom_isolation) {
    
    # Find the time of symptom onset for each trajectory
    symptom_df <- df %>%
      left_join(peak_x) %>%
      group_by(j) %>%
      filter(x >= peak_x) %>%
      slice(1) %>% 
      group_by(j) %>%
      select(j, isolation_start_symptom = x)
    
    # Join df, test df, and symptom df
    df_inf_test_symptom <- symptom_df %>%
      full_join(test_df) %>%
      right_join(df) %>%
      group_by(j) %>%
      filter(x >= ifelse(symptomatic==1,
                         pmin(isolation_start, isolation_start_symptom, na.rm = TRUE),
                         Inf
                         )) %>%
      summarise(area_test_symptom = sum(y))
    
    # Calculate the prevention in infections due to testing, isolation, and symptom onset
    df_res <- df_inf %>%
      left_join(df_inf_test_symptom) %>%
      mutate(area_test_symptom = ifelse(is.na(area_test_symptom), 0, area_test_symptom)) %>% 
      mutate(inf_prevented = area_test_symptom / area_all) 
    
  } else {
    
    # Join df and test df
    df_inf_test <- df %>%
      right_join(test_df) %>%
      group_by(j) %>%
      filter(x >= isolation_start) %>%
      summarise(area_test_symptom = sum(y)) 
    
    # Calculate the prevention in infections due to testing and isolation
    df_res <- df_inf %>%
      left_join(df_inf_test) %>%
      mutate(area_test_symptom = ifelse(is.na(area_test_symptom), 0, area_test_symptom)) %>% 
      mutate(inf_prevented = area_test_symptom / area_all) 
    
  }
  
  df_res %>% mutate(scenario = scenario, symptom_isolation = symptom_isolation, n_tests = n_tests)
  
}

