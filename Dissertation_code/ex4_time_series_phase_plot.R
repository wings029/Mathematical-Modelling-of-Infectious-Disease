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
  list(res, intervention_v)
}

times <- seq(0, 1500, by = 1) 

parms1 <- c(alpha = 0.01134021 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.5 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365,
            beta = 100 / 365, beta_v = 100 / 365, 
            omega = 1 / 600,
            omega_v =  1 / 600,
            theta = 1,
            day_interv = 0, 
            interv_duration = 0,
            cov_v = 0,   
            cov_nv = 0,
            efficacy_v = 0,
            efficacy_nv = 0)  
start1 <- c(Sp = 0.978 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.022, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

# 0% coverage
parms_0 <- parms1
out_0 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_0)
cum_exc_mor_0 <- cumsum(out_0[,"Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"])
peak_Iv_0 <- out_0[,"Ipv"] + out_0[, "Isv"]
R_total_0 <- out_0[,"R"] + out_0[,"Rv"]

# Any scenario to be tested
parms_int <- parms1
parms_int["cov_nv"] <- 0.7
parms_int["efficacy_nv"] <- (0.7*100)^2/10000
parms_int["cov_v"] <- 0
parms_int["efficacy_v"] <- 1
parms_int["day_interv"] <- 80
parms_int["interv_duration"] <- 100
out_int <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int)
cum_exc_mor_int <- cumsum(out_int[,"Ipv"] * parms1["mu_c"] + out_int[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"])
peak_Iv_int <- out_int[,"Ipv"] + out_int[, "Isv"]
R_total_int <- out_int[,"R"] + out_int[,"Rv"]

# The intervention period
shade_df <- data.frame(xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf)

############################################ Peak Iv
prevalence <- data.frame(
  time = out_0[,"time"],
  peak_Iv_0 = peak_Iv_0, 
  scenario_tested = peak_Iv_int
)

Iv_plot <- ggplot(prevalence, aes(x=time)) +
  geom_line(aes(y = peak_Iv_0, color = "0")) +
  geom_line(aes(y = scenario_tested, color = "1")) +
  labs(
    x = "Day",
    y = "Iv",
    color = "Intervention Scenarios"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey3", alpha = 0.2) +
  theme_bw()
Iv_plot

############################################ Cumulative excess mortality
cum_excess_mor <- data.frame(
  time = out_0[,"time"],
  cum_exc_mor_0 = cum_exc_mor_0*1000000, 
  scenario_tested = cum_exc_mor_int*1000000
)

mortality_plot <- ggplot(cum_excess_mor, aes(x=time)) +
  geom_line(aes(y = cum_exc_mor_0, color = "0")) +
  geom_line(aes(y = scenario_tested, color = "1")) +
  labs(
    x = "Day",
    y = "Cumulative Disease Deaths per 1,000,000 Population",
    color = "Intervention Scenario"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey3", alpha = 0.2) +
  theme_bw()
mortality_plot

####################################### R total
R_total <- data.frame(
  time = out_0[,"time"],
  R_total_0 = R_total_0, 
  R_total_int = R_total_int
)

R_plot <- ggplot(R_total, aes(x=time)) +
  geom_line(aes(y = R_total_0, color = "0")) +
  geom_line(aes(y = R_total_int, color = "1")) +
  labs(
    x = "Day",
    y = "R (Total)",
    color = "Intervention Scenario"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey30", alpha = 0.2) +
  theme_bw()
R_plot

####################################### Phase plot
phase_plot_data_0 <- as.data.frame(out_0)
phase_plot_data_0$S_total <- phase_plot_data_0$Sp + phase_plot_data_0$Spv + phase_plot_data_0$Ss + phase_plot_data_0$Ssv
phase_plot_data_0$I_total <- phase_plot_data_0$Ip + phase_plot_data_0$Ipv + phase_plot_data_0$Is + phase_plot_data_0$Isv
phase_plot_data_0$Scenario <- "0"

phase_plot_data_int <- as.data.frame(out_int)
phase_plot_data_int$S_total <- phase_plot_data_int$Sp + phase_plot_data_int$Spv + phase_plot_data_int$Ss + phase_plot_data_int$Ssv
phase_plot_data_int$I_total <- phase_plot_data_int$Ip + phase_plot_data_int$Ipv + phase_plot_data_int$Is + phase_plot_data_int$Isv
phase_plot_data_int$Scenario <- "1"

combined <- rbind(phase_plot_data_0, phase_plot_data_int)

phase_plot <- ggplot(combined, aes(x = S_total, y = I_total, color = Scenario, linetype = Scenario)) +
  geom_path(linewidth = 1) +
  labs(x = "S (Total)",
       y = "I (Total)",
       color = "Scenario",
       linetype = "Scenario") +
  theme_bw()
phase_plot
