library(ggplot2)

# Define the logistic function
sensitivity_function <- function(x) {
  p <- 1 / (1 + exp(-0.006 * (x - 500)))
  return(p)
}

# Generate x values
x_values <- seq(0, 1000, length.out = 1000)

# Generate y values using the logistic function
y_values <- sensitivity_function(x_values)

# Plot the logistic function using ggplot
ggplot(data.frame(x = x_values, y = y_values), aes(x = x, y = y)) +
  geom_line() +
  labs(x = "Days since infection", y = "Sensitivity") +
  theme_minimal()

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

num_trajectories = 1000
num_points = 1000
peak_x_mean = 2.5
peak_x_sd = 1
peak_y_mean = 1000
peak_y_sd = 20
max_x_mean = 10
max_x_sd = 2

# Generate multiple trajectories
trajectories <- lapply(1:1000, function(i) {
  generate_piecewise_linear(num_points = 1000, peak_x_mean, peak_x_sd, peak_y_mean, peak_y_sd, max_x_mean, max_x_sd)
})

# Combine all trajectories into one data frame
trajectories <- do.call(rbind, lapply(1:length(trajectories), function(i) {
  cbind(trajectories[[i]], trajectory = i)
}))

# Plot using ggplot
ggplot(df, aes(x = x, y = y,group = trajectory)) +
  geom_line(alpha = 0.05) +
  ylim(0, 1200) +
  labs(x = "x", y = "y", title = "Multiple Piecewise Linear Curves")


# Function to simulate viral load trajectories with testing
simulate_testing <- function(trajectories, delay) {
  #browser()
  
  #integrate under the curve
  df_inf <- df %>% 
    group_by(trajectory) %>%
    summarise(area_all = sum(y))
  
  #probability of testing positive
  df$test_prob <- sensitivity_function(df$y)
  
  # Simulate testing
  testing_day <- delay
  
  # Apply one test per trajectory on testing day
  test_df <- df %>% 
    group_by(trajectory) %>%
    filter(x >= testing_day) %>%
    slice_min(x) %>% 
    rowwise() %>%
    mutate(result = rbinom(n(), 1, test_prob)) %>% 
    mutate(isolation_start = ifelse(result == 1, x, NA)) %>% 
    select(-x, - test_prob, -y)
  
  # Join df and test df 
 df_inf_test <- df %>% 
   left_join(test_df, by = c("trajectory")) %>% 
   group_by(trajectory) %>%
   filter(x >= isolation_start) %>%
   summarise(area_test = sum(y))
 
 # Calculate the prevention in infections due to testing and isolation
 df_res <- df_inf %>% 
   left_join(df_inf_test, by = "trajectory") %>% 
   mutate(inf_prevented = (area_all - area_test) / area_all) %>% 
   mutate(inf_prevented = ifelse(is.na(inf_prevented), 0, inf_prevented))
  
  return(df_res)
 
}

# Simulate viral load trajectories with testing and delays of 1, 2, and 3 days
data.frame(delay = c(c(1:7),Inf)) %>% 
  mutate(simulation = map(delay, ~simulate_testing(trajectories,
                                                   delay = .x)))%>% 
  unnest(simulation) %>% 
  ggplot(aes(x = factor(delay), y = 1-inf_prevented)) +
    geom_violin() +
    lims(y = c(0, 1)) +
    labs(x = "Delay (days)", y = "Transmission potential averted", title = "Effect of Testing and Isolation on Infection Prevention")



