rm(list = ls())
library(deSolve)
library(ggplot2)
library(cowplot)
library(rootSolve)

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
  cov_i <- parms["cov_i"]
  efficacy <- parms["efficacy"]
  
  # Time-dependent intervention effect
  int_cov <- (t >= day_interv & t < day_lift_interv) * cov_i
  intervention <- (1 - efficacy * int_cov)
  
  # Define birth rate & total population size
  N <- Sp + Ip + R + Ss + Is + Spv + Ipv + Rv + Ssv + Isv
  b <- (mu * (Sp + Ip + R + Ss + Is) + mu_v * (Spv + Ipv + Rv + Ssv + Isv))/N
  
  # Intervention affects transmission
  dSp <- b * N - intervention * beta * Sp * (Ip + Ipv + Is + Isv) / N - (mu + alpha) * Sp  
  dIp <- intervention * beta * Sp * (Ip + Ipv + Is + Isv) / N - (mu + sigma) * Ip
  dR <- sigma * (Ip + Is) - (mu + alpha + omega) * R
  dSs <- omega * R - (mu + alpha) * Ss - intervention * beta * Ss * (Ip + Ipv + Is + Isv) / N
  dIs <- intervention * beta * Ss * (Ip + Ipv + Is + Isv) / N- (mu + sigma) * Is
  dSpv <- alpha * Sp - mu_v * Spv - intervention * beta_v * Spv * (Ip + Ipv + Is + Isv) / N
  dIpv <- intervention * beta_v * Spv * (Ip + Ipv + Is + Isv) / N - (mu_c + mu_v + sigma_v) * Ipv
  dRv <- sigma_v * (Ipv + Isv) + alpha * R - (mu_v + omega_v) * Rv
  dSsv <- omega_v * Rv + alpha * Ss - mu_v * Ssv - intervention * beta_v * Ssv * (Ip + Ipv + Is + Isv) / N
  dIsv <- intervention * beta_v * Ssv * (Ip + Ipv + Is + Isv) / N - ((1 - theta) * mu_c + mu_v + sigma_v) * Isv
  
  results <- c(dSp, dIp, dR, dSs, dIs, dSpv, dIpv, dRv, dSsv, dIsv)
  list(results)
}

#times <- seq(0, 365*3, by=1) #SIR
times <- seq(0, 1500, by=1) #SIRS

### Population 1
parms1 <- c(alpha = 0.01134021 / 365, mu_c = 0.005 / 7.3,
            mu = 0.004 / 365, mu_v = 0.5 / 365,
            sigma = 50 / 365, sigma_v = 50 / 365, 
            beta = 100 / 365, beta_v = 100 / 365, 
            omega = 1 / 600, #change omega to reflect SIR/SIRS dynamics
            omega_v = 1 / 600, 
            theta = 0, #change theta to reflect different immune protection against death
            day_interv = 0, 
            interv_duration = 0,
            cov_i = 0,   
            efficacy = 1)  
start1 <- c(Sp = 0.978 - 0.000001, Ip = 0.000001, R = 0, Ss = 0, Is = 0, Spv = 0.022, Ipv = 0, Rv = 0, Ssv = 0, Isv = 0)

# Equilibrium solutions
equilibrium <- runsteady(y = start1, times = c(0, 36500), func = SIRS_int, parms = parms1) #simulate for 100 years
equilibrium

# 0% coverage
parms_0 <- parms1
parms_0["cov_i"] <- 0  
out_0 <- ode(y = start1, times = times, func = SIRS_int, parms = parms_0)
cum_exc_mor_0 <- cumsum(out_0[,"Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"])
peak_Iv_0 <- out_0[,"Ipv"] + out_0[, "Isv"]
R_total_0 <- out_0[,"R"] + out_0[,"Rv"]
death_0 <- out_0[, "Ipv"] * parms1["mu_c"] + out_0[, "Isv"] * parms1["mu_c"] * (1 - parms1["theta"])

# Any scenario to be tested
parms_int <- parms1
parms_int["cov_i"] <- 0.25
parms_int["day_interv"] <- 400
parms_int["interv_duration"] <- 100
out_int <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int)
cum_exc_mor_int <- cumsum(out_int[,"Ipv"] * parms1["mu_c"] + out_int[, "Isv"] * (1 - parms1["theta"]) * parms1["mu_c"])
peak_Iv_int <- out_int[,"Ipv"] + out_int[, "Isv"]
R_total_int <- out_int[,"R"] + out_int[,"Rv"]
death_int <- out_int[, "Ipv"] * parms1["mu_c"] + out_int[, "Isv"] * parms1["mu_c"] * (1 - parms1["theta"])

plot(out_int[, "Sp"] + out_int[, "Ss"] + out_int[, "Spv"] + out_int[, "Ssv"])
plot(out_int[, "Ip"] + out_int[, "Is"] + out_int[, "Ipv"] + out_int[, "Isv"])
plot(out_int[, "R"] + out_int[, "Rv"])

# Define a dataframe for the intervention period
shade_df <- data.frame(xmin = 400, xmax = 500, ymin = -Inf, ymax = Inf)

############################################ Peak Iv
prevalence <- data.frame(
  time = out_0[,"time"],
  peak_Iv_0 = peak_Iv_0, 
  scenario_tested = peak_Iv_int
)

Iv_plot <- ggplot(prevalence, aes(x=time)) +
  geom_line(aes(y = peak_Iv_0, color = "0% coverage")) +
  geom_line(aes(y = scenario_tested, color = "20% coverage at day 40")) +
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

############################################ Peak deaths
death_data <- data.frame(
  time = out_0[,"time"],
  death_0 = death_0 * 1000000,
  death_int = death_int * 1000000
)

death_plot <- ggplot(death_data, aes(x=time)) +
  geom_line(aes(y = death_0, color = "0% coverage")) +
  geom_line(aes(y = death_int, color = "20% coverage at day 40")) +
  labs(
    x = "Day",
    y = "Disease Deaths per 1,000,000 Population",
    color = "Intervention Scenarios"
  ) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12)) +
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey3", alpha = 0.2) +
  theme_bw()
death_plot

############################################ Cumulative excess mortality
cum_excess_mor <- data.frame(
  time = out_0[,"time"],
  cum_exc_mor_0 = cum_exc_mor_0*1000000, 
  scenario_tested = cum_exc_mor_int*1000000
)

mortality_plot <- ggplot(cum_excess_mor, aes(x=time)) +
  geom_line(aes(y = cum_exc_mor_0, color = "0% coverage")) +
  geom_line(aes(y = scenario_tested, color = "30% coverage at day 80")) +
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
  geom_line(aes(y = R_total_0, color = "0% coverage")) +
  geom_line(aes(y = R_total_int, color = "40% coverage at day 75")) +
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
phase_plot_data_0$Scenario <- "0% coverage"

phase_plot_data_int <- as.data.frame(out_int)
phase_plot_data_int$S_total <- phase_plot_data_int$Sp + phase_plot_data_int$Spv + phase_plot_data_int$Ss + phase_plot_data_int$Ssv
phase_plot_data_int$I_total <- phase_plot_data_int$Ip + phase_plot_data_int$Ipv + phase_plot_data_int$Is + phase_plot_data_int$Isv
phase_plot_data_int$Scenario <- "40% coverage at day 75"

combined <- rbind(phase_plot_data_0, phase_plot_data_int)

phase_plot <- ggplot(combined, aes(x = S_total, y = I_total, color = Scenario, linetype = Scenario)) +
  geom_path(linewidth = 1) +
  labs(x = "S (Total)",
       y = "I (Total)",
       color = "Scenario",
       linetype = "Scenario") +
  theme_bw()
phase_plot

