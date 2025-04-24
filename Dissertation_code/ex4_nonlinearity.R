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

# Ranges of intervention coverage and efficacy
# For vulnerable pool, efficacy is always 1
cov_v_values <- seq(0, 1, by = 0.05)  
# For non-vulnerable pool, efficacy is a function of coverage
cov_nv_values <- seq(0, 1, by = 0.05)  
eff_nv_values <- (cov_nv_values*100)^2/10000
combination <- data.frame(coverage_v = cov_v_values, coverage_nv = cov_nv_values,
                          efficacy_nv = eff_nv_values)

######################################### Start day = 60, duration = 100 #############################################
results_60_100 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                             efficacy_v = numeric(), efficacy_nv = numeric(), 
                             cov_eff_product_nv = numeric(),
                             peak_Iv = numeric(), cum_exc_mor = numeric())

for (i in seq_len(nrow(combination))) {
  parms1["cov_v"] <- combination$coverage_v[i]
  parms1["efficacy_v"] <- 1
  parms1["day_interv"] <- 60
  parms1["interv_duration"] <- 100
  
  for (j in seq_len(nrow(combination))) {
    parms1["cov_nv"] <- combination$coverage_nv[j]
    parms1["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start1, times = times, func = SIRS_int, parms = parms1)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms1["mu_c"] + out[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_60_100 <- rbind(results_60_100, data.frame(
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

# Percentage reduction and percentage change
results_60_100$perc_reduction_Iv <- (results_60_100[1, "peak_Iv"]-results_60_100$peak_Iv)/results_60_100[1, "peak_Iv"]*100
results_60_100$perc_reduction_mortality <- (results_60_100[1, "cum_exc_mor"] - results_60_100$cum_exc_mor)/results_60_100[1, "cum_exc_mor"] * 100
results_60_100$perc_change_Iv <- -results_60_100$perc_reduction_Iv
results_60_100$perc_change_mortality <- -results_60_100$perc_reduction_mortality

### Time series
# Parameters of best intervention for reducing mortality
index <- which.max(results_60_100$perc_reduction_mortality)
results_60_100[index, 1] #cov_v = 1
results_60_100[index, 2] #cov_nv = 0
results_60_100[index, 4] #efficacy_nv = 0

# Parameters of best intervention for reducing peak Iv
index <- which.max(results_60_100$perc_reduction_Iv)
results_60_100[index, 1] #cov_v = 0.9
results_60_100[index, 2] #cov_nv = 0.55
results_60_100[index, 4] #efficacy_nv = 0.3025

# Baseline scenario - no intervention
parms_no_int <- parms1
parms_no_int["cov_v"] <- 0
parms_no_int["cov_nv"] <- 0
out_0 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_no_int)
excess_mortality_0 <- out_0[,"Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_0 <- cumsum(excess_mortality_0)
total_Iv_0 <- out_0[,"Ipv"] + out_0[, "Isv"]

# Intervention scenario
parms_int_60_100 <- parms1

# # params of mortality-reducing intervention
# parms_int_60_100["cov_v"] <- 1
# parms_int_60_100["efficacy_v"] <- 1
# parms_int_60_100["cov_nv"] <- 0
# parms_int_60_100["efficacy_nv"] <- 0

# params of peak Iv-reducing intervention
parms_int_60_100["cov_v"] <- 0.9
parms_int_60_100["efficacy_v"] <- 1
parms_int_60_100["cov_nv"] <- 0.55
parms_int_60_100["efficacy_nv"] <- 0.3025

parms_int_60_100["day_interv"] <- 60
parms_int_60_100["interv_duration"] <- 100
out_int_60_100 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int_60_100)
excess_mortality_int_60_100 <- out_int_60_100[,"Ipv"] * parms1["mu_c"] + out_int_60_100[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_int_60_100 <- cumsum(excess_mortality_int_60_100)
total_Iv_int_60_100 <- out_int_60_100[,"Ipv"] + out_int_60_100[, "Isv"]

# Shade the intervention period
shade_df <- data.frame(xmin = 60, xmax = 160, ymin = -Inf, ymax = Inf) 

time_series_data_60_100 <- data.frame(
  time = out_0[,"time"],
  cum_exc_mor_0 = cum_exc_mor_0*1000000, 
  cum_exc_mor_int = cum_exc_mor_int_60_100*1000000,
  total_Iv_0 = total_Iv_0,
  total_Iv_int_60_100 = total_Iv_int_60_100
)

# Plotting - mortality plot
time_series_60_100 <- ggplot(time_series_data_60_100, aes(x=time)) +
  geom_line(aes(y = cum_exc_mor_0, color = "0% coverage")) +
  geom_line(aes(y = cum_exc_mor_int, color = "100% coverage_v + 0% coverage_nv")) +
  labs(
    x = "Day",
    y = "Cumulative Disease Deaths per 1,000,000 Population",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_60_100

# Plotting - peak Iv
time_series_Iv_60_100 <- ggplot(time_series_data_60_100, aes(x=time)) +
  geom_line(aes(y = total_Iv_0, color = "0% coverage")) +
  geom_line(aes(y = total_Iv_int_60_100, color = "90% coverage_v + 55% coverage_nv")) +
  labs(
    x = "Day",
    y = "Iv",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_Iv_60_100

# Phase plot
phase_plot_data_0 <- as.data.frame(out_0)
phase_plot_data_0$S_total <- phase_plot_data_0$Sp + phase_plot_data_0$Spv + phase_plot_data_0$Ss + phase_plot_data_0$Ssv
phase_plot_data_0$I_total <- phase_plot_data_0$Ip + phase_plot_data_0$Ipv + phase_plot_data_0$Is + phase_plot_data_0$Isv
phase_plot_data_0$Scenario <- "0% coverage"

phase_plot_60_100 <- as.data.frame(out_int_60_100)
phase_plot_60_100$S_total <- phase_plot_60_100$Sp + phase_plot_60_100$Spv + phase_plot_60_100$Ss + phase_plot_60_100$Ssv
phase_plot_60_100$I_total <- phase_plot_60_100$Ip + phase_plot_60_100$Ipv + phase_plot_60_100$Is + phase_plot_60_100$Isv
phase_plot_60_100$Scenario <- "90% cov_v + 55% cov_nv"

combined <- rbind(phase_plot_data_0, phase_plot_60_100)

phase_plot_60_100 <- ggplot(combined, aes(x = S_total, y = I_total, color = Scenario, linetype = Scenario)) +
  geom_path(linewidth = 1) +
  labs(x = "S (Total)",
       y = "I (Total)",
       color = "Scenario",
       linetype = "Scenario") +
  theme_bw()
phase_plot_60_100

######################################### Start day = 80, duration = 100 #############################################
results_80_100 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
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
    
    results_80_100 <- rbind(results_80_100, data.frame(
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

# Percentage reduction & percentage change
results_80_100$perc_reduction_Iv <- (results_80_100[1, "peak_Iv"]-results_80_100$peak_Iv)/results_80_100[1, "peak_Iv"]*100
results_80_100$perc_reduction_mortality <- (results_80_100[1, "cum_exc_mor"] - results_80_100$cum_exc_mor)/results_80_100[1, "cum_exc_mor"] * 100
results_80_100$perc_change_Iv <- -results_80_100$perc_reduction_Iv
results_80_100$perc_change_mortality <- -results_80_100$perc_reduction_mortality

### Time series
# Parameters of best intervention in reducing mortality
index <- which.max(results_80_100$perc_reduction_mortality)
results_80_100[index, 1] #cov_v = 1
results_80_100[index, 2] #cov_nv = 0
results_80_100[index, 4] #efficacy_nv = 0

# Parameters of best intervention in reducing peak Iv
index <- which.max(results_80_100$perc_reduction_Iv)
results_80_100[index, 1] #cov_v = 1
results_80_100[index, 2] #cov_nv = 0.3
results_80_100[index, 4] #efficacy_nv = 0.09

# Baseline scenario (no intervention)
parms_no_int <- parms1
parms_no_int["cov_v"] <- 0
parms_no_int["cov_nv"] <- 0
out_0 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_no_int)
excess_mortality_0 <- out_0[,"Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_0 <- cumsum(excess_mortality_0)
total_Iv_0 <- out_0[,"Ipv"] + out_0[, "Isv"]

# Intervention scenario 
parms_int_80_100 <- parms1

# # params of mortality-reducing intervention
# parms_int_80_100["cov_v"] <- 1
# parms_int_80_100["efficacy_v"] <- 1
# parms_int_80_100["cov_nv"] <- 0

# params of peak Iv-reducing intervention
parms_int_80_100["cov_v"] <- 1
parms_int_80_100["efficacy_v"] <- 1
parms_int_80_100["cov_nv"] <- 0.3
parms_int_80_100["efficacy_nv"] <- 0.09

parms_int_80_100["day_interv"] <- 80
parms_int_80_100["interv_duration"] <- 100
out_int_80_100 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int_80_100)
excess_mortality_int_80_100 <- out_int_80_100[,"Ipv"] * parms1["mu_c"] + out_int_80_100[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_int_80_100 <- cumsum(excess_mortality_int_80_100)
total_Iv_int_80_100 <- out_int_80_100[,"Ipv"] + out_int_80_100[, "Isv"] 

# Shade the intervention period
shade_df <- data.frame(xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf) 

time_series_data_80_100 <- data.frame(
  time = out_0[,"time"],
  cum_exc_mor_0 = cum_exc_mor_0*1000000, 
  cum_exc_mor_int = cum_exc_mor_int_80_100*1000000,
  total_Iv_0 = total_Iv_0,
  total_Iv_int_80_100 = total_Iv_int_80_100
)

# Plotting - mortality
time_series_80_100 <- ggplot(time_series_data_80_100, aes(x=time)) +
  geom_line(aes(y = cum_exc_mor_0, color = "0% coverage")) +
  geom_line(aes(y = cum_exc_mor_int, color = "100% coverage_v + 0% coverage_nv")) +
  labs(
    x = "Day",
    y = "Cumulative Disease Deaths per 1,000,000 Population",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_80_100

# Plotting - peak Iv
time_series_Iv_80_100 <- ggplot(time_series_data_80_100, aes(x=time)) +
  geom_line(aes(y = total_Iv_0, color = "0% coverage")) +
  geom_line(aes(y = total_Iv_int_80_100, color = "100% coverage_v + 0% coverage_nv")) +
  labs(
    x = "Day",
    y = "Iv",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_Iv_80_100

# Phase plot
phase_plot_data_0 <- as.data.frame(out_0)
phase_plot_data_0$S_total <- phase_plot_data_0$Sp + phase_plot_data_0$Spv + phase_plot_data_0$Ss + phase_plot_data_0$Ssv
phase_plot_data_0$I_total <- phase_plot_data_0$Ip + phase_plot_data_0$Ipv + phase_plot_data_0$Is + phase_plot_data_0$Isv
phase_plot_data_0$Scenario <- "0% coverage"

phase_plot_80_100 <- as.data.frame(out_int_80_100)
phase_plot_80_100$S_total <- phase_plot_80_100$Sp + phase_plot_80_100$Spv + phase_plot_80_100$Ss + phase_plot_80_100$Ssv
phase_plot_80_100$I_total <- phase_plot_80_100$Ip + phase_plot_80_100$Ipv + phase_plot_80_100$Is + phase_plot_80_100$Isv
phase_plot_80_100$Scenario <- "100% cov_v + 0% cov_nv"

combined <- rbind(phase_plot_data_0, phase_plot_80_100)

phase_plot_80_100 <- ggplot(combined, aes(x = S_total, y = I_total, color = Scenario, linetype = Scenario)) +
  geom_path(linewidth = 1) +
  labs(x = "S (Total)",
       y = "I (Total)",
       color = "Scenario",
       linetype = "Scenario") +
  theme_bw()
phase_plot_80_100

######################################### Heatmaps - changing start day (60 vs 80) ##############################################
# Iv
Iv_scale_min <- min(results_60_100$perc_change_Iv, 
                    results_80_100$perc_change_Iv)
Iv_scale_max <- max(results_60_100$perc_change_Iv, 
                    results_80_100$perc_change_Iv)

Iv_plot_60_100 <- ggplot(results_60_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_60_100 

Iv_plot_80_100 <- ggplot(results_80_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_80_100 

# Mortality
Mortality_scale_min <- min(results_60_100$perc_change_mortality, 
                           results_80_100$perc_change_mortality)
Mortality_scale_max <- max(results_60_100$perc_change_mortality, 
                           results_80_100$perc_change_mortality)

Mortality_plot_60_100 <- ggplot(results_60_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_60_100

Mortality_plot_80_100 <- ggplot(results_80_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_80_100

######################################### Start day = 60, duration = 200 #############################################
results_60_200 <- data.frame(coverage_v = numeric(), coverage_nv = numeric(),
                             efficacy_v = numeric(), efficacy_nv = numeric(), 
                             cov_eff_product_nv = numeric(),
                             peak_Iv = numeric(), cum_exc_mor = numeric())

# Loop through different coverage values for both vulnerable and non-vulnerable pool
for (i in seq_len(nrow(combination))) {
  parms1["cov_v"] <- combination$coverage_v[i]
  parms1["efficacy_v"] <- 1
  parms1["day_interv"] <- 60
  parms1["interv_duration"] <- 200
  
  for (j in seq_len(nrow(combination))) {
    parms1["cov_nv"] <- combination$coverage_nv[j]
    parms1["efficacy_nv"] <- combination$efficacy_nv[j]
    
    out <- ode(y = start1, times = times, func = SIRS_int, parms = parms1)
    
    peak_Iv <- max(out[,"Ipv"] + out[,"Isv"])
    excess_mortality <- out[,"Ipv"] * parms1["mu_c"] + out[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
    cum_exc_mor <- tail(cumsum(excess_mortality), n = 1)
    
    results_60_200 <- rbind(results_60_200, data.frame(
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

results_60_200$perc_reduction_Iv <- (results_60_200[1, "peak_Iv"]-results_60_200$peak_Iv)/results_60_200[1, "peak_Iv"]*100
results_60_200$perc_reduction_mortality <- (results_60_200[1, "cum_exc_mor"] - results_60_200$cum_exc_mor)/results_60_200[1, "cum_exc_mor"] * 100
results_60_200$perc_change_Iv <- -results_60_200$perc_reduction_Iv
results_60_200$perc_change_mortality <- -results_60_200$perc_reduction_mortality

### Time series
# Parameters of best intervention in reducing mortality
index <- which.max(results_60_200$perc_reduction_mortality)
results_60_200[index, 1] #cov_v = 1
results_60_200[index, 2] #cov_nv = 0
results_60_200[index, 4] #efficacy_nv = 0

# Parmaeters of best intervention in reducing peak Iv
index <- which.max(results_60_200$perc_reduction_Iv)
results_60_200[index, 1] #cov_v = 0.8
results_60_200[index, 2] #cov_nv = 0.65
results_60_200[index, 4] #efficacy_nv = 0.4225

# Baseline scenario
parms_no_int <- parms1
parms_no_int["cov_v"] <- 0
parms_no_int["cov_nv"] <- 0
out_0 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_no_int)
excess_mortality_0 <- out_0[,"Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_0 <- cumsum(excess_mortality_0)
total_Iv_0 <- out_0[,"Ipv"] + out_0[, "Isv"]

# Intervention scenario
parms_int_60_200 <- parms1

# params of mortality-reducing intervention
# parms_int_60_200["cov_v"] <- 1
# parms_int_60_200["efficacy_v"] <- 1
# parms_int_60_200["cov_nv"] <- 0

# params of peak Iv-reducing intervention
parms_int_60_200["cov_v"] <- 0.8
parms_int_60_200["efficacy_v"] <- 1
parms_int_60_200["cov_nv"] <- 0.65
parms_int_60_200["efficacy_nv"] <- 0.4225

parms_int_60_200["day_interv"] <- 60
parms_int_60_200["interv_duration"] <- 200
out_int_60_200 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int_60_200)
excess_mortality_int_60_200 <- out_int_60_200[,"Ipv"] * parms1["mu_c"] + out_int_60_200[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"]
cum_exc_mor_int_60_200 <- cumsum(excess_mortality_int_60_200)
total_Iv_int_60_200 <- out_int_60_200[,"Ipv"] + out_int_60_200[, "Isv"] 

# Shade the intervention period
shade_df <- data.frame(xmin = 60, xmax = 260, ymin = -Inf, ymax = Inf) 

time_series_data_60_200 <- data.frame(
  time = out_0[,"time"],
  cum_exc_mor_0 = cum_exc_mor_0*1000000, 
  cum_exc_mor_int = cum_exc_mor_int_60_200*1000000,
  total_Iv_0 = total_Iv_0,
  total_Iv_int_60_200 = total_Iv_int_60_200
)

# Plotting - mortality
time_series_60_200 <- ggplot(time_series_data_60_200, aes(x=time)) +
  geom_line(aes(y = cum_exc_mor_0, color = "0% coverage")) +
  geom_line(aes(y = cum_exc_mor_int, color = "100% cov_v + 0% cov_nv")) +
  labs(
    x = "Day",
    y = "Cumulative Disease Deaths per 1,000,000 Population",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_60_200

# Plotting - peak Iv
time_series_Iv_60_200 <- ggplot(time_series_data_60_200, aes(x=time)) +
  geom_line(aes(y = total_Iv_0, color = "0% coverage")) +
  geom_line(aes(y = total_Iv_int_60_200, color = "85% coverage_v + 25% coverage_nv")) +
  labs(
    x = "Day",
    y = "Iv",
    color = "0% coverage vs. best intervention"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
time_series_Iv_60_200

# Phase plot
phase_plot_data_0 <- as.data.frame(out_0)
phase_plot_data_0$S_total <- phase_plot_data_0$Sp + phase_plot_data_0$Spv + phase_plot_data_0$Ss + phase_plot_data_0$Ssv
phase_plot_data_0$I_total <- phase_plot_data_0$Ip + phase_plot_data_0$Ipv + phase_plot_data_0$Is + phase_plot_data_0$Isv
phase_plot_data_0$Scenario <- "0% coverage"

phase_plot_60_200 <- as.data.frame(out_int_60_200)
phase_plot_60_200$S_total <- phase_plot_60_200$Sp + phase_plot_60_200$Spv + phase_plot_60_200$Ss + phase_plot_60_200$Ssv
phase_plot_60_200$I_total <- phase_plot_60_200$Ip + phase_plot_60_200$Ipv + phase_plot_60_200$Is + phase_plot_60_200$Isv
phase_plot_60_200$Scenario <- "100% cov_v + 0% cov_nv"

combined <- rbind(phase_plot_data_0, phase_plot_60_200)

phase_plot_60_200 <- ggplot(combined, aes(x = S_total, y = I_total, color = Scenario, linetype = Scenario)) +
  geom_path(linewidth = 1) +
  labs(x = "S (Total)",
       y = "I (Total)",
       color = "Scenario",
       linetype = "Scenario") +
  theme_bw()
phase_plot_60_200
######################################### Heatmaps - changing duration (100 vs 200) #####################
### Peak Iv
Iv_scale_min <- min(results_60_100$perc_change_Iv,
                    results_60_200$perc_change_Iv)
Iv_scale_max <- max(results_60_100$perc_change_Iv,
                    results_60_200$perc_change_Iv)

Iv_plot_60_100 <- ggplot(results_60_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_60_100 

Iv_plot_60_200 <- ggplot(results_60_200, aes(x = coverage_v, y = coverage_nv, fill = perc_change_Iv)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", name = "% Change in Peak Iv",
                      limits = c(Iv_scale_min, Iv_scale_max)) +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Iv_plot_60_200 

### Mortality
Mortality_scale_min <- min(results_60_100$perc_change_mortality,
                           results_60_200$perc_change_mortality)
Mortality_scale_max <- max(results_60_100$perc_change_mortality,
                           results_60_200$perc_change_mortality)

Mortality_plot_60_100 <- ggplot(results_60_100, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_60_100

Mortality_plot_60_200 <- ggplot(results_60_200, aes(x = coverage_v, y = coverage_nv, fill = perc_change_mortality)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(Mortality_scale_min, 0, Mortality_scale_max)),
                       limits = c(Mortality_scale_min, Mortality_scale_max),
                       name = "% Change in Cumulative Disease Mortality") +
  labs(x = "Coverage (Vulnerable)", y = "Coverage (Non-vulnerable)") +
  theme_bw()
Mortality_plot_60_200
