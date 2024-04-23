pacman::p_load(tidyverse,tidytable,update = T)

# Define the logistic function
sensitivity_function <- function(x) {
  p <- 1 / (1 + exp(-0.006 * (x - 500)))
  return(p)
}

# Function to generate piecewise linear data
generate_piecewise_linear <- function(num_points, peak_x_mean, peak_x_sd, peak_y_mean, peak_y_sd, max_x_mean, max_x_sd) {
  # Duration
  duration <- rnorm(1, mean = max_x_mean, sd = max_x_sd)
  
  # Generate breakpoints
  breakpoints <- c(0, rnorm(1, mean = peak_x_mean, sd = peak_x_sd), duration)
  
  # Generate values
  values <- c(0, rnorm(1, mean = peak_y_mean, sd = peak_y_sd), 0)
  
  # Interpolate between breakpoints
  piecewise_linear <- function(x) {
    ifelse(x <= breakpoints[2],
           values[1] + (values[2] - values[1]) / (breakpoints[2] - breakpoints[1]) * (x - breakpoints[1]),
           values[2] + (values[3] - values[2]) / (breakpoints[3] - breakpoints[2]) * (x - breakpoints[2]))
  }
  
  # Generate x values
  x_values <- seq(0, duration, length.out = num_points)
  
  # Generate y values using the piecewise linear function
  y_values <- piecewise_linear(x_values)
  
  return(data.frame(x = x_values, y = y_values))
}

# Function to simulate viral load trajectories with testing and symptom onset
simulate_trajectories <- function(params) {
  num_trajectories <- params$num_trajectories
  num_points <- params$num_points
  peak_x_mean <- params$peak_x_mean
  peak_x_sd <- params$peak_x_sd
  peak_y_mean <- params$peak_y_mean
  peak_y_sd <- params$peak_y_sd
  max_x_mean <- params$max_x_mean
  max_x_sd <- params$max_x_sd
  
  # Generate multiple trajectories
  trajectories <- lapply(1:num_trajectories, function(i) {
    generate_piecewise_linear(num_points = num_points, peak_x_mean, peak_x_sd, peak_y_mean, peak_y_sd, max_x_mean, max_x_sd)
  })
  
  # Combine all trajectories into one data frame
  trajectories <- do.call(rbind, lapply(1:length(trajectories), function(i) {
    cbind(trajectories[[i]], trajectory = i)
  }))
  
  return(trajectories)
}

# Function to simulate testing, isolation due to symptoms, and isolation due to positive test
simulate_isolation <- function(trajectories, delay, symptom_isolation) {
 
  df <- trajectories
  
  # Integrate under the curve
  df_inf <- df %>%
    group_by(trajectory) %>%
    summarise(area_all = sum(y))
  
  # Probability of testing positive
  df$test_prob <- sensitivity_function(df$y)
  
  # Simulate testing and symptom onset
  testing_day <- delay
  
  # Find the peak viral load for each trajectory
  peak_x <- df %>%
    group_by(trajectory) %>%
    filter(y == max(y)) %>% 
    select(trajectory,peak_x = x)
  
  # Apply one test per trajectory on testing day
  test_df <- df %>%
    group_by(trajectory) %>%
    filter(x >= testing_day) %>%
    slice(1) %>%
    rowwise() %>%
    mutate(result = rbinom(n(), 1, test_prob)) %>%
    mutate(isolation_start_test = ifelse(result == 1, x, NA)) %>%
    select(-x, -test_prob, -y)
  
  # Simulate isolation due to symptom onset if symptom_isolation is TRUE
  if (symptom_isolation) {
    symptom_df <- df %>%
      left_join(peak_x, by = "trajectory") %>%
      group_by(trajectory) %>%
      mutate(isolation_start_symptom = ifelse(x >= peak_x, peak_x, NA))
    
    # Join df, test df, and symptom df
    df_inf_test_symptom <- df %>%
      left_join(test_df, by = c("trajectory")) %>%
      left_join(symptom_df) %>%
      group_by(trajectory) %>%
      filter(x >= pmin(isolation_start_test, isolation_start_symptom, na.rm = TRUE)) %>%
      summarise(area_test_symptom = sum(y))
    
    # Calculate the prevention in infections due to testing, isolation, and symptom onset
    df_res <- df_inf %>%
      left_join(df_inf_test_symptom, by = "trajectory") %>%
      mutate(inf_prevented = (area_all - area_test_symptom) / area_all) %>%
      mutate(inf_prevented = ifelse(is.na(inf_prevented), 0, inf_prevented))
    
  } else {
    
    # Join df and test df
    df_inf_test <- df %>%
      left_join(test_df, by = c("trajectory")) %>%
      group_by(trajectory) %>%
      filter(x >= isolation_start_test) %>%
      summarise(area_test = sum(y))
    
    # Calculate the prevention in infections due to testing and isolation
    df_res <- df_inf %>%
      left_join(df_inf_test, by = "trajectory") %>%
      mutate(inf_prevented = (area_all - area_test) / area_all) %>%
      mutate(inf_prevented = ifelse(is.na(inf_prevented), 0, inf_prevented))
    
  }
  
  return(df_res)
}

# Example usage
flu_params <- list(
  num_trajectories = 100,
  num_points = 1000,
  peak_x_mean = 2.5,
  peak_x_sd = 1,
  peak_y_mean = 1000,
  peak_y_sd = 20,
  max_x_mean = 10,
  max_x_sd = 2
)

trajectories <- simulate_trajectories(flu_params)

# Simulate viral load trajectories with testing, isolation due to symptoms, and delays of 1, 2, and 3 days
data.frame(delay = c(c(1:7), Inf)) %>% 
  crossing(symptom_isolation = c(TRUE,FALSE)) %>%
  mutate(simulation = map2(delay, symptom_isolation, ~ simulate_isolation(trajectories, delay = .x, symptom_isolation = .y))) %>%
  unnest(simulation) %>%
  ggplot(aes(x = factor(delay), y = 1 - inf_prevented, fill = factor(symptom_isolation))) +
  geom_violin() +
  facet_wrap(~symptom_isolation)+
  lims(y = c(0, 1)) +
  labs(x = "Delay (days)", y = "Transmission potential averted", title = "Effect of Testing, Isolation, and Symptom Onset on Infection Prevention",
       fill = "Symptom Isolation") +
  theme_minimal()
