gen_screening_draws <- function(x){
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x){
 # browser()
  # generate required times for screening 
  
  # what's the probability of detection at each test time given a value of CT?
  x_ <- x %>%
   inner_join(trajectories$models, by=c("sec_idx"="idx","type")) %>% 
    select(-data) %>% 
    mutate(test_t  = pmax(test_t,sec_traced_t),
           test_q  = test_t-sec_exposed_t,
           ct      = pmap_dbl(.f=calc_sensitivity,list(model=m,x=test_q))) %>% 
    mutate(detection_range = case_when(sens_LFA=="higher" ~ as.numeric(as.character(cut(ct,
                                                                          breaks=c(-Inf,27,30,35,Inf),
                                                                          labels=c("0.95","0.65","0.3","0")))),
                                       sens_LFA=="lower" ~ as.numeric(as.character(cut(ct,
                                                                                       breaks = c(-Inf,20,25,30,35,Inf),
                                                                                       labels=c("0.824","0.545","0.083","0.053","0")))))) %>% 
    mutate(test_p=case_when(assay=="PCR"&ct<35~1,
                            assay=="PCR"&ct>=35~0,
                            assay=="LFA"~detection_range,
                            TRUE~NA_real_)) %>% 
    mutate(screen      = runif(n(), 0, 1)) %>% 
    mutate(test_label  = detector(pcr = test_p,  u = screen)) 
 
  #browser()
  return(x_)
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_trajectories <- function(n_cases){
  #simulate CT trajectories
  #browser()
  traj <- data.frame(idx=1:n_cases) %>% 
    crossing(start=0,type=c("symptomatic","asymptomatic") %>% 
               factor(x = .,
                      levels = .,
                      ordered = T)) %>% 
    mutate(u = runif(n(),0,1)) %>%
    # duration from: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30172-5/fulltext
    # scaling of asymptomatics taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(end=case_when(type == "symptomatic"  ~ qnormTrunc(p = u, mean=17, 
                                                             sd=approx_sd(15.5,18.6), min = 0),
                         type == "asymptomatic" ~ qnormTrunc(p = u, mean=17*0.6, 
                                                             sd=approx_sd(15.5,18.6), min = 0))) %>% 
    # incubation period from https://bmjopen.bmj.com/content/10/8/e039652.info
    mutate(onset_t=qlnormTrunc(p = u,
                               meanlog=1.63,
                               sdlog=0.5,
                               min = 0,
                               max=end)) %>% 
    pivot_longer(cols = -c(idx,type,u),
                 values_to = "x") %>% 
    # peak CT taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(y=case_when(name=="start"   ~ 40,
                       name=="end"     ~ 40,
                       name=="onset_t" ~ rnorm(n=n(),mean=22.3,sd=4.2))) 
  
  models <- traj %>%
    nest(data = -c(idx,type,u)) %>%  
    dplyr::mutate(
      # Perform loess calculation on each individual 
      m  = purrr::map(data, ~splinefunH(x = .x$x, y = .x$y,
                                       m = c(0,0,0))),
      rx = purrr::map(data, ~range(.x$x)),
      ry = purrr::map(data, ~range(.x$y)),
      
      inf_period=purrr::pmap(.f = infectious_period,
                             .l = list(model = m, 
                                       rx = rx, # range in x needed for when we're outside the range of the sampled end time
                                       ry = ry))) %>% # range in y is for debugging
    unnest_wider(inf_period)
  
  return(list(traj=traj,models=models))
}

# Create viral load trajectories for a given number of sims
make_trajectories <- function(
    n_sims = 100,
    asymp_parms = asymp_fraction,
    variant_info, 
    browsing = FALSE
){
  
  if (browsing) browser()
  
  set.seed(seed)
  #simulate CT trajectories
  
  inf <- rbbinom(n = n_sims,
                 size=1,
                 alpha = asymp_parms$shape1,
                 beta  = asymp_parms$shape2) %>%
    as_tidytable() %>%
    rename.("asymptomatic"=x) %>%
    mutate.(sim=row_number.(),
            asymptomatic=FALSE)#as.logical(asymptomatic))
  
  traj <- inf %>% 
    crossing.(start=0) %>% 
    mutate.(
      prolif=case_when.(heterogen_vl~rnormTrunc(n=n(),mean=mean_prolif,
                                                sd=sd_prolif,min = 1,max=14),
                        TRUE~median(rnormTrunc(n=n(),mean=mean_prolif,
                                               sd=sd_prolif,min = 1,max=14))),
      clear=case_when.(heterogen_vl~rnormTrunc(n=n(), mean=mean_clear, 
                                               sd=sd_clear, min = 1,max=30),
                       TRUE~median(rnormTrunc(n=n(), mean=mean_clear, 
                                              sd=sd_clear, min = 1,max=30))),
      end=prolif+clear,
      onset_t=prolif+rnorm(n=n(),mean = 2,sd=1.5)
    ) %>%
    select.(-c(mean_prolif, sd_prolif, mean_clear, sd_clear,clear)) %>%
    pivot_longer.(cols = -c(sim,variant,onset_t, asymptomatic, heterogen_vl,
                            mean_peakvl,sd_peakvl),
                  values_to = "x") %>%
    mutate.(y=case_when.(name=="start" ~ 40,
                         name=="end"   ~ 40,
                         name=="prolif"~case_when.(heterogen_vl~rnormTrunc(n=n(),
                                                                           mean=mean_peakvl,
                                                                           sd=sd_peakvl,min=0,max=40),
                                                   TRUE~median(rnormTrunc(n=n(),
                                                                          mean=mean_peakvl,
                                                                          sd=sd_peakvl,min=0,max=40))))) %>% 
    select.(-c(mean_peakvl,sd_peakvl))
  
  
  models <- traj %>%
    nest.(data = -c(sim,variant,onset_t,asymptomatic,heterogen_vl)) %>%  
    mutate.(
      # Perform approxfun on each set of points
      m  = map.(data, ~approxfun(x=.x$x,y=.x$y))) 
  
  #cannot pivot wider with "m" column - extract and rejoin
  x_model <- models %>% 
    select.(-data)
  
  models <- models %>% 
    select.(-m) %>% 
    unnest.(data,.drop=F) %>%  
    select.(-c(y)) %>% 
    pivot_wider.(names_from=name,values_from = x) %>% 
    left_join.(x_model) %>%
    select.(c(sim, variant, heterogen_vl, onset_t, prolif, start, end, m)) %>% 
    arrange.(sim)
  
}

inf_curve_func <- function(m,start=0,end=30,interval=1,trunc_t){
  #browser()
  x <- tidytable(t=seq(start,end,by=interval)) %>% 
    mutate.(ct=m(t),
            vl=convert_Ct_logGEML(ct))
  
  return(infectiousness=x)
}

calc_sensitivity <- function(model, x){
  #browser()
  if(!is.na(x)){
    s <- model(x)
  } else {
    s <- NA_real_
  }
  
  return(s)
}

#### Main Model ----
run_model <- function(testing_scenarios, scenarios, contact_dat=contact_data, browsing=F){
  
  if(browsing){browser()}
  
  #### Generate infections of hh (household) contacts ####
  indiv_params <- traj %>%  
    select.(-m) %>%
    crossing.(heterogen_contacts = unique(scenarios$heterogen_contacts),
              period=unique(scenarios %>% 
                              mutate.(period=fct_drop(period)) %>% 
                              pull.(period))) %>%   
    mutate.(hh_contacts=ifelse(heterogen_contacts,
                               sample_filter(condition = period, 
                                             df = contact_dat, 
                                             col="e_home", 
                                             n=n()),
                               round(mean_filter(condition = period, 
                                                 df = contact_data, 
                                                 col="e_home"))),
            .by=c(period)) 
  
  indiv_params_long <- indiv_params %>% 
    left_join.(traj_)
  
  #simulate infections (and keep first instance)
  hh_infections <- indiv_params_long %>% 
    uncount.(hh_contacts,.id="id",.remove = F) %>% 
    mutate.(hh_duration = case_when.(heterogen_contacts ~ sample(contacts_hh_duration,
                                                                 size=n(),replace=T),
                                     TRUE               ~ median(contacts_hh_duration)),
            infected    = rbernoulli(n(),p=culture_p*hh_duration)) %>% 
    filter.(infected==T) %>% 
    slice.(min(t), .by=c(all_of(key_grouping_var),hh_contacts,id)) %>% 
    count.(t,all_of(key_grouping_var),hh_contacts,name = "hh_infected") %>% 
    arrange.(sim)
  
  #### Calculate nhh infections ####
  
  nhh_infections <- indiv_params_long %>% 
    
    # Sample daily contacts
    mutate.(nhh_contacts = ifelse(heterogen_contacts,
                                  sample_filter(condition = period,
                                                df=contact_dat,
                                                col="e_other",n=n()),
                                  round(mean_filter(condition = period,
                                                    df=contact_dat,
                                                    col="e_other"))),
            .by=c(period)) %>% 
    
    # Simulate infections 
    uncount.(nhh_contacts,.remove = F) %>% 
    mutate.(nhh_duration = case_when.(heterogen_contacts ~ sample(contacts_nhh_duration,
                                                                  size=n(),replace=T),
                                      TRUE               ~ median(contacts_nhh_duration)),
            nhh_infected = rbernoulli(n=n(),p = culture_p*nhh_duration)) %>% 
    summarise.(nhh_infected=sum(nhh_infected),.by=c(t,all_of(key_grouping_var),nhh_contacts,test)) %>% 
    
    # Testing: determine if and when testing + isolating by specified sampling frequency, adherence  
    right_join.(testing_scenarios) %>% 
    mutate.(
      test_day = case_when.((t - begin_testing) %% sampling_freq == 0 ~ TRUE,
                            nhh_contacts > event_size ~ TRUE,
                            TRUE ~ FALSE)) %>% 
    mutate.(
      earliest_pos = min(t[test&test_day]),
      test_iso = t>=earliest_pos & self_iso_test,
      .by=c(all_of(key_grouping_var),prop_self_iso_test,self_iso_test,begin_testing,sampling_freq,event_size)) %>%
    filter.(test_iso==F) %>% 
    select.(everything(),-test_iso,-test,-earliest_pos,-test_day) 
  
  # Join nhh and hh contacts and summarise
  processed_infections <- indiv_params_long %>% 
    right_join.(testing_scenarios) %>% 
    left_join.(hh_infections) %>% 
    left_join.(nhh_infections) %>% 
    replace_na.(list(hh_infected=0,nhh_infected=0,nhh_contacts=0)) %>% 
    arrange.(period,lower_inf_thresh) %>% 
    mutate.(
      total_contacts = nhh_contacts+hh_contacts,
      total_infections=nhh_infected+hh_infected)
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, trajectories,n_sec_cases){
  
  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  res <- lapply(names(props), 
                function(x){
                  filter(trajectories, type == x) %>%
                    sample_frac(., size = props[[x]])
                })
  
  res <- do.call("rbind",res) %>% 
    sample_n(n_sec_cases)
  
}


make_released_quantiles <- function(x, vars){
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  
  dots <- append(dots1, dots2)
  
  x_count <- x %>%
    dplyr::ungroup(.) %>%
    dplyr::select(!!! dots) %>%
    dplyr::group_by_all(.) %>%
    dplyr::count(.) 
  
  x_count %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup %>%
    as.list %>%
    map(unique) %>%
    expand.grid %>%
    dplyr::left_join(x_count) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    tidyr::nest(data = c(sim, n)) %>%
    dplyr::mutate(
      Q = purrr::map(
        .x = data,
        .f = ~quantile(.x$n, probs = probs)),
      M = purrr::map_dbl(
        .x = data, 
        .f = ~mean(.x$n))) %>%
    tidyr::unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup(.)
}

make_quantiles <- function(x, y_var,vars, sum = TRUE){
  #browser()
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  y_var <- as.name(y_var)
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, y_var) %>%
      group_by_at(.vars = vars(-y_var)) %>% 
      summarise(!!y_var := sum(!!y_var, na.rm=T),n=n()) %>% 
      mutate(!!y_var:=!!y_var/n)
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, !! y_var) 
  
  x_days %>%
    nest(data = c(!!y_var, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x[[y_var]],
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x[[y_var]]))) %>%
    unnest_wider(Q) %>% 
    # mutate(Q = purrr::map(.x = data, ~bayestestR::hdi( .x[[y_var]],
    #                                             ci = c(0.95,0.5))),
    #        Mean = map_dbl(.x = data, ~mean(.x[[y_var]])),
    #        MAP = map_dbl(.x=data, ~as.numeric(bayestestR::point_estimate(.x[[y_var]],centrality="MAP")))) %>% 
    #unnest(Q) %>% 
    #pivot_wider(names_from = CI,values_from=c(CI_low,CI_high)) %>% 
    dplyr::select(-data) 
  
}


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  gamma2mv(ans$estimate[["shape"]],
           ans$estimate[["rate"]])
}


rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  #browser()
  p_b <- pgamma(q = b, shape = shape, rate = rate)
  p_a <- pgamma(q = a, shape = shape, rate = rate)
  
  u   <- runif(n = n, min = p_a, max = p_b)
  q   <- qgamma(p = u, shape = shape, rate = rate)
  
  return(q)
}


check_unique_values <- function(df, vars){
  # given a data frame and a vector of variables to be used to facet or group, 
  # which ones have length < 1?
  
  l <- lapply(X = vars, 
              FUN =function(x){
                length(unique(df[, x]))
              })
  
  vars[l > 1]
  
}

test_times <- function(multiple_tests,tests,tracing_t,sec_exposed_t,quar_dur,sampling_freq = 1, max_tests = 14, n_tests){
  #browser()
  
  if(multiple_tests){
    test_timings <- data.frame(test_t = seq(from=tracing_t,to=(tracing_t+max_tests-1),by=sampling_freq)) %>% 
    mutate(test_no = paste0("test_", row_number())) %>% 
    mutate(have_test = row_number()<=n_tests) %>% 
    filter(have_test)
  } else {
    test_timings <- data.frame(test_t=sec_exposed_t+quar_dur) %>% 
      mutate(test_no = paste0("test_", row_number())) %>% 
      mutate(have_test = tests)
  }
  
  return(test_timings)
}

earliest_pos <- function(df){
  #browser()
  x_q <- unlist(df[(df$test_label),"test_q"])
  
  if (length(x_q) == 0L){
    return(Inf)
  } else {
    return(min(x_q))
  }
}

earliest_pos2 <- function(df){
  #browser()
  x_q <- df %>% filter(test_label==TRUE) 
  
  if (nrow(x_q) == 0L){
    return(data.frame(test_no="None",test_p=0,test_q=Inf))
  } else {
    return((x_q %>% select(test_no,test_p,test_q) %>% slice_min(test_q)))
  }
}


infectious_period <- function(model, rx, ry){
  #browser()
  newdata <- data.frame(x = seq(from = 0, to = 30, 0.1))
  
  newdata$y_pred <- model(newdata$x)
  
  
  newdata <- mutate(newdata,
                    y_pred = ifelse(x > rx[2], 40, y_pred))
  
  newdata %>%
    mutate(infectious = y_pred <= 30) %>%
    filter(infectious) %>%
    summarise(inf_start = min(x),
              inf_end   = max(x))
  
}

calc_overlap <- function(x){
#browser()
  x <-
    x %>% 
    rowwise() %>% 
    mutate(
      symp_overlap = case_when(
        type == "symptomatic" ~ Overlap(
          x = c(inf_start,
                inf_end),
          y = c(onset_q,
                symp_end_q)
        ) * adhering_iso,
        type == "asymptomatic" ~ 0
      ),
      symp_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, symp_overlap),
      quar_overlap = case_when(
        #If symptomatic before tracing, no quarantine
        !multiple_tests & onset_q < traced_q & type == "symptomatic" ~ 0, 
        #If symptomatic during quarantine, quarantine ends at onset
        !multiple_tests & onset_q > traced_q & onset_q < quar_end_q & type == "symptomatic" ~ 
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(traced_q,
                onset_q)
        ) * adhering_quar,
        #If symptomatic after quarantine, quarantine is completed in full
        !multiple_tests & onset_q > traced_q & onset_q > quar_end_q & type == "symptomatic" ~ 
          Overlap(
            x = c(inf_start,
                  inf_end),
            y = c(traced_q,
                  quar_end_q)
          ) * adhering_quar,
        #If asymptomatic, quarantine is completed in full
        !multiple_tests & type == "asymptomatic" ~ Overlap(
          x = c(inf_start,
                inf_end),
          y = c(traced_q,
                quar_end_q)
        ) * adhering_quar,
        multiple_tests ~ 0
      ),
      quar_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, quar_overlap),
      test_overlap = case_when(
        #If testing and symptoms occur before the test, no-post test isolation
        tests & onset_q < earliest_q &type =="symptomatic" ~ 0,
        #If testing and symptoms occur after a positive test, end post-positive test isolation and begin symptomatic
        tests & onset_q > earliest_q & onset_q < test_iso_end_q & type == "symptomatic" ~
          Overlap(
        x = c(inf_start,
              inf_end),
        y = c(earliest_q,
              onset_q)
      ) * adhering_iso,
      #If testing and symptoms occur after the end of post-positive test isolation, complete full self-isolation
      tests & onset_q > earliest_q & onset_q > test_iso_end_q & type == "symptomatic" ~
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(earliest_q,
                test_iso_end_q)
        ) * adhering_iso,
      #If testing and asymptomatic, complete full self-isolation
      tests & type == "asymptomatic" ~
        Overlap(
          x = c(inf_start,
                inf_end),
          y = c(earliest_q,
                test_iso_end_q)
        ) * adhering_iso,
      !tests ~ 0),
      test_overlap = ifelse(is.infinite(inf_start) &
                              is.infinite(inf_end), 0, test_overlap),
      max_overlap = case_when(
        #If quarantine and testing, time in quarantine is sum of quarantine  symptomatic and test self-isolation durations
        tests &
          !multiple_tests ~ symp_overlap + quar_overlap + test_overlap,
        #If daily testing, time in isolation is sum of symptomatic and test self-isolation durations
        tests &
          multiple_tests  ~ symp_overlap + test_overlap,
        #If quarantine only, sum of symptomatic and quarantine overlap
        !tests            ~ symp_overlap + quar_overlap
      )
    )
  
}

approx_sd <- function(x1, x2){
      (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
   }
