rm(list = ls())
library(deSolve)

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
  N <- (Sp + Ip + R + Ss + Is + Spv + Ipv + Rv + Ssv + Isv)
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
  list(results, intervention, N)
}

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

# Any intervention scenario to be tested
parms_int <- parms1
parms_int["cov_i"] <- 0.3
parms_int["day_interv"] <- 80
parms_int["interv_duration"] <- 100
out_int <- ode(y = start1, times = times, func = SIRS_int, parms = parms_int)
colnames(out_int)[12] <- "intervention"
colnames(out_int)[13] <- "N"

# Calculate real-time Rt (effective reproduction number)
out_df <- as.data.frame(out_int)
out_df$Rt <- parms1["beta"] / parms1["sigma"] * ((out_df$Sp + out_df$Spv + out_df$Ss + out_df$Ssv)/out_df$N)
View(out_df)
