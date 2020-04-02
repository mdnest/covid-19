# covid-19
Source files for COVID-19 project for Math 699 Stochastic Processes, Spring 2020




```R
require("deSolve","ggplot2")
library("deSolve")
library("ggplot2")
```

# Deterministic SIR Model (ODEs)

The source code for the deterministic model is SIR.R. It uses 4th order Runge-Kutta ODE solver.

$$ \frac{dS}{dT} $$

```R
#Define ODE system
f <- function(t, y, parms) {
  S = y[1]
  I = y[2]
  R = y[3]
  
  beta = parms[1]
  gamma = parms[2]
  
  dS = -beta*S*I/N
  dI = beta*S*I/N - gamma*I
  dR = gamma*I

  return(list(c(dS, dI, dR)))
}
```


```R
# Initial condition
N = 1000
S0 = 999
I0 = 1
R0 = 0
x0 = c(S=S0, I=I0, R=R0)
```


```R
# Parameters
beta=1/2
gamma=1/5
mu = c(beta, gamma)
```


```R
# Time interval
t_min = 0
t_max = 100
t_eval <- seq(t_min,t_max)
```


```R
# Run Runge-Kutta algorithm to solve ODE
soln = rk4(y=x0, times=t_eval, func=f, parms=mu)
```


```R
# Plot results
ggplot(data.frame(soln), aes(x=time)) + 
  labs(title='SIR Model Time Series', x='Time', y='Population count') +
  geom_line(aes(y=S, col='Susceptible'), size=1) +
  geom_line(aes(y=I, col='Infected'), size=1) +
  geom_line(aes(y=R, col='Recovered'), size=1)
```


![png](/SIR_plot.png)


# Stochastic SIR Model (Markov Chain)

Uses same parameters and initial conditions as the deterministic model.


```R
# Define Markov chain
simulate_stochastic_SIR <- function(y, times, dt, parms) {
    
    t_min = min(times)
    t_max = max(times)
    
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
    
    df = data.frame(time=seq(t_min,t_max,dt))
    
    df$S = S
    df$I = I
    df$R = R
    
    time_series <- subset(df, time %in% times)
    return(time_series)
}
```


```R
# Run model
stochastic_SIR_time_series <- simulate_stochastic_SIR(y=x0, times=t_eval, dt=0.01, parms=mu)
```


```R
# Plot results
ggplot(data=stochastic_SIR_time_series, aes(x=time)) + 
  labs(title='Stochastic SIR Model Time Series', x='Time', y='Population count') +
  geom_line(aes(y=S, col='Susceptible'), size=1) +
  geom_line(aes(y=I, col='Infected'), size=1) +
  geom_line(aes(y=R, col='Recovered'), size=1)
```


![png](/stochastic_SIR_plot.png)

