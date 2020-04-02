
require("ggplot2")
library("ggplot2")

# Define stochastic model
simulate_stochastic_SIR <- function(y, t_eval, dt, parms) {
  
  t_min = min(t_eval)
  t_max = max(t_eval)
  
  k_max = 1 + floor((t_max-t_min)/dt)
  
  S = numeric(k_max)
  I = numeric(k_max)
  R = numeric(k_max)
  
  S0 = y[1]
  I0 = y[2]
  R0 = y[3]
  
  S[1] = S0
  I[1] = I0
  R[1] = R0
  
  beta  = parms[1]
  gamma = parms[2]
  
  for (k in 1:(k_max-1)) {
    s = S[k]
    i = I[k]
    r = R[k]
    
    S[k+1] = s
    I[k+1] = i
    R[k+1] = r
    
    random = runif(2)
    
    if (random[1] < beta*s*i/N*dt) {
      S[k+1] = s-1
      I[k+1] = i+1
    }
    if (random[2] < gamma*i*dt) {
      I[k+1] = i-1
      R[k+1] = r+1
    }
  }
  
  full_time_series = data.frame(time=seq(t_min,t_max,dt))
  
  full_time_series$S = S
  full_time_series$I = I
  full_time_series$R = R
  
  time_series <- subset(full_time_series, time %in% t_eval)
  return(time_series)
}

# Initial condition
N = 1000
S0 = 999
I0 = 1
R0 = 0
x0 = c(S=S0, I=I0, R=R0)

# Parameters
beta=1/2
gamma=1/5
mu = c(beta, gamma)

# Time interval
t_min = 0
t_max = 100
t_eval <- seq(t_min,t_max)

stochastic_SIR_time_series <- simulate_stochastic_SIR(y=x0, t_eval=seq(0,100), dt=0.01, parms=mu)

# Plot results
ggplot(stochastic_SIR_time_series, aes(x=time)) + 
  labs(title='Stochastic SIR Model Time Series', x='Time', y='Population count') +
  geom_line(aes(y=S, col='Susceptible'), size=1) +
  geom_line(aes(y=I, col='Infected'), size=1) +
  geom_line(aes(y=R, col='Recovered'), size=1)

