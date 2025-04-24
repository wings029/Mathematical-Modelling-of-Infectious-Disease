rm(list = ls())
library(deSolve)
library(ggplot2)
library(cowplot)
library(reshape2)

### SIRS model with interventions
SIRS_int <- function(t, y, parms) {
  Sp <- y[1] #primary infection
  Ip <- y[2] #primary infection
  R <- y[3]
  Ss <- y[4] #secondary infection
  Is <- y[5] #secondary infection
  Spv <- y[6] #primary infection
  Ipv <- y[7] #primary infection
  Rv <- y[8]
  Ssv <- y[9] #secondary infection
  Isv <- y[10] #secondary infection
  
  mu <- parms["mu"]
  mu_v <- parms["mu_v"]
  mu_c <- parms["mu_c"]
  alpha <- parms["alpha"]
  sigma <- parms["sigma"]
  sigma_v <- parms["sigma_v"]
  beta <- parms["beta"]
  beta_v <- parms["beta_v"]
  omega <- parms["omega"] 
  omega_v <- parms["omega_v"]
  theta <- parms["theta"] #level of immune protection against death
  
  # Intervention-related params
  day_interv <- parms["day_interv"]
  interv_duration <- parms["interv_duration"]
  day_lift_interv <- day_interv + interv_duration 
  
  # Coverage and efficacy params
  cov_v <- parms["cov_v"]
  cov_nv <- parms["cov_nv"]
  efficacy_v <- parms["efficacy_v"]
  efficacy_nv <- parms["efficacy_nv"]
  
  # Time-dependent intervention effect
  int_cov_v <- (t >= day_interv & t < day_lift_interv) * cov_v
  int_cov_nv <- (t >= day_interv & t < day_lift_interv) * cov_nv
  intervention_v <- (1 - efficacy_v * int_cov_v)
  intervention_nv <- (1 - efficacy_nv * int_cov_nv)
  
  # Define birth rate & total population size
  N <- Sp + Ip + R + Ss + Is + Spv + Ipv + Rv + Ssv + Isv
  b <- (mu * (Sp + Ip + R + Ss + Is) + mu_v * (Spv + Ipv + Rv + Ssv + Isv))/N
  
  # Intervention affects transmission
  dSp <- b * N - intervention_nv * beta * Sp * (Ip + Ipv + Is + Isv) / N - (mu + alpha) * Sp  
  dIp <- intervention_nv * beta * Sp * (Ip + Ipv + Is + Isv) / N - (mu + sigma) * Ip
  dR <- sigma * (Ip + Is) - (mu + alpha + omega) * R
  dSs <- omega * R - (mu + alpha) * Ss - intervention_nv * beta * Ss * (Ip + Ipv + Is + Isv) / N
  dIs <- intervention_nv * beta * Ss * (Ip + Ipv + Is + Isv) / N- (mu + sigma) * Is
  dSpv <- alpha * Sp - mu_v * Spv - intervention_v * beta_v * Spv * (Ip + Ipv + Is + Isv) / N
  dIpv <- intervention_v * beta_v * Spv * (Ip + Ipv + Is + Isv) / N - (mu_c + mu_v + sigma_v) * Ipv
  dRv <- sigma_v * (Ipv + Isv) + alpha * R - (mu_v + omega_v) * Rv
  dSsv <- omega_v * Rv + alpha * Ss - mu_v * Ssv - intervention_v * beta_v * Ssv * (Ip + Ipv + Is + Isv) / N
  dIsv <- intervention_v * beta_v * Ssv * (Ip + Ipv + Is + Isv) / N - ((1 - theta) * mu_c + mu_v + sigma_v) * Isv
  
  res <- c(dSp, dIp, dR, dSs, dIs, dSpv, dIpv, dRv, dSsv, dIsv)
  list(res)
}

times <- seq(0, 1500, by = 1) 

# Ranges of intervention coverage and efficacy
# For vulnerable pool, efficacy is always 1
cov_v_values <- seq(0, 1, by = 0.05)  
# For non-vulnerable pool, efficacy is a function of coverage
cov_nv_values <- seq(0, 1, by = 0.05)  
eff_nv_values <- (cov_nv_values*100)^2/10000
combination <- data.frame(coverage_v = cov_v_values, coverage_nv = cov_nv_values,
                          efficacy_nv = eff_nv_values)

########################################### Population 1 - A ##########################################################
## Population 1 params
parms1 <- c(alpha = 0.01134021 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.5 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365, 
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start1 <- c(Sp = 0.978 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.022, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop1 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms1["cov_v"] <- combination$coverage_v[i]
  parms1["efficacy_v"] <- 1
  parms1["day_interv"] <- 80
  parms1["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms1["cov_nv"] <- combination$coverage_nv[j]
    parms1["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start1, times = times, func = SIRS_int, parms = parms1)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms1["mu_c"] + out[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop1 <- rbind(results_pop1, data.frame(
      coverage_v = parms1["cov_v"],
      coverage_nv = parms1["cov_nv"],
      efficacy_v = parms1["efficacy_v"],
      efficacy_nv = parms1["efficacy_nv"],
      cov_eff_product_nv = parms1["cov_nv"] * parms1["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

# Calculate percentage change in peak Iv_total & cumulative excess mortality
results_pop1$perc_reduction_Iv <- (results_pop1[1, "peak_Iv"]-results_pop1$peak_Iv)/results_pop1[1, "peak_Iv"]*100
results_pop1$perc_reduction_mortality <- (results_pop1[1, "cum_exc_mor"] - results_pop1$cum_exc_mor)/results_pop1[1, "cum_exc_mor"] * 100
results_pop1$perc_change_Iv <- -results_pop1$perc_reduction_Iv 
results_pop1$perc_change_mortality <- -results_pop1$perc_reduction_mortality 

########################################### Population 2 - A ##########################################################
## Population 2 params
parms2 <- c(alpha = 0.01571429 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.05 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365,
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start2 <- c(Sp = 0.761 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.239, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop2 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms2["cov_v"] <- combination$coverage_v[i]
  parms2["efficacy_v"] <- 1
  parms2["day_interv"] <- 80
  parms2["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms2["cov_nv"] <- combination$coverage_nv[j]
    parms2["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start2, times = times, func = SIRS_int, parms = parms2)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms2["mu_c"] + out[, "Isv"] * (1 - parms2["theta"]) * parms2["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop2 <- rbind(results_pop2, data.frame(
      coverage_v = parms2["cov_v"],
      coverage_nv = parms2["cov_nv"],
      efficacy_v = parms2["efficacy_v"],
      efficacy_nv = parms2["efficacy_nv"],
      cov_eff_product_nv = parms2["cov_nv"] * parms2["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

results_pop2$perc_reduction_Iv <- (results_pop2[1, "peak_Iv"]-results_pop2$peak_Iv)/results_pop2[1, "peak_Iv"]*100
results_pop2$perc_reduction_mortality <- (results_pop2[1, "cum_exc_mor"] - results_pop2$cum_exc_mor)/results_pop2[1, "cum_exc_mor"] * 100
results_pop2$perc_change_Iv <- -results_pop2$perc_reduction_Iv 
results_pop2$perc_change_mortality <- -results_pop2$perc_reduction_mortality 

########################################### Population 3 - A ##########################################################
## Population 3 params
parms3 <- c(alpha = 0.01189189 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.2 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365,
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start3 <- c(Sp = 0.944 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.056, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop3 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms3["cov_v"] <- combination$coverage_v[i]
  parms3["efficacy_v"] <- 1
  parms3["day_interv"] <- 80
  parms3["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms3["cov_nv"] <- combination$coverage_nv[j]
    parms3["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start3, times = times, func = SIRS_int, parms = parms3)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms3["mu_c"] + out[, "Isv"] * (1 - parms3["theta"]) * parms3["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop3 <- rbind(results_pop3, data.frame(
      coverage_v = parms3["cov_v"],
      coverage_nv = parms3["cov_nv"],
      efficacy_v = parms3["efficacy_v"],
      efficacy_nv = parms3["efficacy_nv"],
      cov_eff_product_nv = parms3["cov_nv"] * parms3["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

results_pop3$perc_reduction_Iv <- (results_pop3[1, "peak_Iv"]-results_pop3$peak_Iv)/results_pop3[1, "peak_Iv"]*100
results_pop3$perc_reduction_mortality <- (results_pop3[1, "cum_exc_mor"] - results_pop3$cum_exc_mor)/results_pop3[1, "cum_exc_mor"] * 100
results_pop3$perc_change_Iv <- -results_pop3$perc_reduction_Iv 
results_pop3$perc_change_mortality <- -results_pop3$perc_reduction_mortality 

######################################### Population 1 - A+C ##############################################
parms4 <- c(alpha = 0.0185 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.037 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365,
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start4 <- c(Sp = 0.667 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.333, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop4 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms4["cov_v"] <- combination$coverage_v[i]
  parms4["efficacy_v"] <- 1
  parms4["day_interv"] <- 80
  parms4["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms4["cov_nv"] <- combination$coverage_nv[j]
    parms4["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start4, times = times, func = SIRS_int, parms = parms4)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms4["mu_c"] + out[, "Isv"] * (1 - parms4["theta"]) * parms4["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop4 <- rbind(results_pop4, data.frame(
      coverage_v = parms4["cov_v"],
      coverage_nv = parms4["cov_nv"],
      efficacy_v = parms4["efficacy_v"],
      efficacy_nv = parms4["efficacy_nv"],
      cov_eff_product_nv = parms4["cov_nv"] * parms4["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

results_pop4$perc_reduction_Iv <- (results_pop4[1, "peak_Iv"]-results_pop4$peak_Iv)/results_pop4[1, "peak_Iv"]*100
results_pop4$perc_reduction_mortality <- (results_pop4[1, "cum_exc_mor"] - results_pop4$cum_exc_mor)/results_pop4[1, "cum_exc_mor"] * 100
results_pop4$perc_change_Iv <- -results_pop4$perc_reduction_Iv 
results_pop4$perc_change_mortality <- -results_pop4$perc_reduction_mortality 

######################################### Population 2 - A+C ##############################################
parms5 <- c(alpha = 0.02475 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.027 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365,
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start5 <- c(Sp = 0.522 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.478, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop5 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms5["cov_v"] <- combination$coverage_v[i]
  parms5["efficacy_v"] <- 1
  parms5["day_interv"] <- 80
  parms5["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms5["cov_nv"] <- combination$coverage_nv[j]
    parms5["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start5, times = times, func = SIRS_int, parms = parms5)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms5["mu_c"] + out[, "Isv"] * (1 - parms5["theta"]) * parms5["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop5 <- rbind(results_pop5, data.frame(
      coverage_v = parms5["cov_v"],
      coverage_nv = parms5["cov_nv"],
      efficacy_v = parms5["efficacy_v"],
      efficacy_nv = parms5["efficacy_nv"],
      cov_eff_product_nv = parms5["cov_nv"] * parms5["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

results_pop5$perc_reduction_Iv <- (results_pop5[1, "peak_Iv"]-results_pop5$peak_Iv)/results_pop5[1, "peak_Iv"]*100
results_pop5$perc_reduction_mortality <- (results_pop5[1, "cum_exc_mor"] - results_pop5$cum_exc_mor)/results_pop5[1, "cum_exc_mor"] * 100
results_pop5$perc_change_Iv <- -results_pop5$perc_reduction_Iv 
results_pop5$perc_change_mortality <- -results_pop5$perc_reduction_mortality 

######################################### Population 3 - A+C ##############################################
parms6 <- c(alpha = 0.01419149 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.0667 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365,
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 30, 
            interv_duration = 100,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start6 <- c(Sp = 0.825 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.175, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

### Start day = 80, duration = 100 
results_pop6 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                           efficacy_v = numeric(), efficacy_nv = numeric(), 
                           cov_eff_product_nv = numeric(),
                           peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms6["cov_v"] <- combination$coverage_v[i]
  parms6["efficacy_v"] <- 1
  parms6["day_interv"] <- 80
  parms6["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms6["cov_nv"] <- combination$coverage_nv[j]
    parms6["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start6, times = times, func = SIRS_int, parms = parms6)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms6["mu_c"] + out[, "Isv"] * (1 - parms6["theta"]) * parms6["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_pop6 <- rbind(results_pop6, data.frame(
      coverage_v = parms6["cov_v"],
      coverage_nv = parms6["cov_nv"],
      efficacy_v = parms6["efficacy_v"],
      efficacy_nv = parms6["efficacy_nv"],
      cov_eff_product_nv = parms6["cov_nv"] * parms6["efficacy_nv"],
      peak_Iv = peak_Iv,
      cum_exc_mor = cum_exc_mor
    ))
  }
}

results_pop6$perc_reduction_Iv <- (results_pop6[1, "peak_Iv"]-results_pop6$peak_Iv)/results_pop6[1, "peak_Iv"]*100
results_pop6$perc_reduction_mortality <- (results_pop6[1, "cum_exc_mor"] - results_pop6$cum_exc_mor)/results_pop6[1, "cum_exc_mor"] * 100
results_pop6$perc_change_Iv <- -results_pop6$perc_reduction_Iv 
results_pop6$perc_change_mortality <- -results_pop6$perc_reduction_mortality 

######################################### Heatmaps ##############################################
# Iv
Iv_scale_min <- min(results_pop1$perc_change_Iv,
                    results_pop2$perc_change_Iv,
                    results_pop3$perc_change_Iv,
                    results_pop4$perc_change_Iv,
                    results_pop5$perc_change_Iv,
                    results_pop6$perc_change_Iv)
Iv_scale_max <- max(results_pop1$perc_change_Iv,
                    results_pop2$perc_change_Iv,
                    results_pop3$perc_change_Iv,
                    results_pop4$perc_change_Iv,
                    results_pop5$perc_change_Iv,
                    results_pop6$perc_change_Iv)

Iv_plot_pop1 <- ggplot(results_pop1, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop1 

Iv_plot_pop2 <- ggplot(results_pop2, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop2 

Iv_plot_pop3 <- ggplot(results_pop3, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop3 

Iv_plot_pop4 <- ggplot(results_pop4, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop4

Iv_plot_pop5 <- ggplot(results_pop5, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop5 

Iv_plot_pop6 <- ggplot(results_pop6, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_pop6

# Mortality
Mortality_scale_min <- min(results_pop1$perc_change_mortality,
                           results_pop2$perc_change_mortality,
                           results_pop3$perc_change_mortality,
                           results_pop4$perc_change_mortality,
                           results_pop5$perc_change_mortality,
                           results_pop6$perc_change_mortality)
Mortality_scale_max <- max(results_pop1$perc_change_mortality,
                           results_pop2$perc_change_mortality,
                           results_pop3$perc_change_mortality,
                           results_pop4$perc_change_mortality,
                           results_pop5$perc_change_mortality,
                           results_pop6$perc_change_mortality)

Mortality_plot_pop1 <- ggplot(results_pop1, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop1

Mortality_plot_pop2 <- ggplot(results_pop2, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop2

Mortality_plot_pop3 <- ggplot(results_pop3, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop3

Mortality_plot_pop4 <- ggplot(results_pop4, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop4

Mortality_plot_pop5 <- ggplot(results_pop5, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop5

Mortality_plot_pop6 <- ggplot(results_pop6, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_pop6



