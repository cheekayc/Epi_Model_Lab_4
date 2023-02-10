Lab 4
================
2023-02-08

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2

    ## Warning: package 'ggplot2' was built under R version 4.2.2

    ## Warning: package 'readr' was built under R version 4.2.2

    ## Warning: package 'purrr' was built under R version 4.2.2

    ## Warning: package 'dplyr' was built under R version 4.2.2

    ## Warning: package 'stringr' was built under R version 4.2.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    ## Warning: package 'deSolve' was built under R version 4.2.2

# PART 1 SIR model with demography

The SIR model with demography provides a way to study endemic diseases,
in particular what proportion of the population are susceptible or
infected when the disease become endemic. Do the following simulations
and report your findings:

Write a basic SIR function.

``` r
SIR = function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    dS = -beta*S*I/N;
    dI = beta*S*I/N - gamma*I;
    
    # return the rate of change
    list(c(dS,dI))  
  }) # end with(as.list...)
}
```

Write a function for SIR model with demography.

``` r
SIRdem = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    dS = mu*N - beta*S*I/N - mu*S;
    dI = beta*S*I/N - gamma*I - mu*I;
    
    # return the rate of change
    list(c(dS,dI))
  }) # end with(as.list...)
}
```

Set initial conditions and parameters.

``` r
N = 1e5 
I0 = 20
S0 = 0.2*N
state = c(S = S0, I = I0)
# parameters
beta = 520/365; gamma = 1/7; mu = 1/(70*365)
```

Assemble parameters set for each model

``` r
if(T){
  
  parameters.SIR = c(beta = beta, gamma = gamma);
  
  parameters.SIRdem = c(beta = beta, gamma = gamma, mu = mu)
}

times = seq(0, 365*100, by = 1)
```

Call the functions.

``` r
# for SIR model
sim.SIR = ode(y = state, times = times, func = SIR, parms = parameters.SIR);
```

    ## DLSODA-  At T (=R1), too much accuracy requested  
    ##       for precision of machine..  See TOLSF (=R2) 
    ## In above message, R1 = 8673.96, R2 = nan
    ## 

``` r
s.SIR = sim.SIR[,'S']/N  # %S
i.SIR = sim.SIR[,'I']/N  # %I

# for SIRdem model
sim.SIRdem = ode(y = state, times = times, func = SIRdem, parms = parameters.SIRdem);
s.SIRdem = sim.SIRdem[,'S']/N
i.SIRdem = sim.SIRdem[,'I']/N
```

Plot %I vs time for both models and compare.

``` r
if(T){
  par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(1.8, .5, 0))
  plot(sim.SIR[, 'time'], i.SIR, ylab = '% I', xlab = 'Time (day)', main = 'SIR model', type = 'l', col = 'blue', lwd = 1)
  plot(sim.SIRdem[, 'time'], i.SIRdem, ylab = '% I', xlab = 'Time (day)', main = 'SIR with demography', type = 'l', col = 'red', lwd = 1)
  
  # to see just the *first* few years: 
  # set the limit for x-axis using xlim=c(xmin,xmax)
  plot(sim.SIR[, 'time'], i.SIR, ylab = '% I', xlab = 'Time (day)', main = 'SIR model', xlim = c(1, 365*10), type = 'l', col = 'blue')
  plot(sim.SIRdem[, 'time'], i.SIRdem, ylab = '% I', xlab = 'Time (day)', main = 'SIR with demography', xlim = c(1, 365*10), type = 'l', col = 'red')
  
  # to see just the *last* few years: 
  # set the limit for x-axis using xlim=c(xmin,xmax)
  plot(sim.SIR[, 'time'], i.SIR, ylab = '% I', xlab = 'Time (day)', main = 'SIR model', xlim = c(365*90+1, 365*100), type = 'l', col = 'blue')
  plot(sim.SIRdem[, 'time'], i.SIRdem, ylab = '% I', xlab = 'Time (day)', main = 'SIR with demography', xlim = c(365*90+1, 365*100), type = 'l', col = 'red')
}
```

![](Lab-4_files/figure-gfm/plot%20graphs-1.png)<!-- -->![](Lab-4_files/figure-gfm/plot%20graphs-2.png)<!-- -->![](Lab-4_files/figure-gfm/plot%20graphs-3.png)<!-- -->

What are the %S (S/N), and %I (I/N) at the end of the simulation?

``` r
s.SIRdem[36501]
```

    ## [1] 0.1001976

``` r
i.SIRdem[36501]
```

    ## [1] 0.0002458337

# PART 2: SIS model

Write a function for SIS model.

``` r
SIS = function(t, state, parameters){
  with(as.list(c(state,parameters)),{
  
    dS = -beta*S*I/N + gamma*I
    dI = beta*S*I/N - gamma*I
    
    list(c(dS,dI))
  })
}
```

Set initial conditions and parameters.

``` r
N = 1e5
I0 = 10
S0 = N - I0
state = c(S = S0, I = I0)
parameters = c(beta = 0.5, gamma = 0.3)

times = seq(0, 100, by = 1)
```

Call the function:

``` r
sim = ode(y = state, times = times, func = SIS, parms = parameters)
s = sim[, 'S']/N
i = sim[, 'I']/N
```

Plot S vs time and I vs time.

``` r
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(1.8, .5, 0))
  plot(sim[, 'time'], sim[, 'S'], ylab = '# of susceptible', xlab = 'Time (day)', main = 'Susceptible vs. Time graph', type = 'l', col = 'green', lwd = 1)
  plot(sim[, 'time'], sim[, 'I'], ylab = '# of infectious', xlab = 'Time (day)', main = 'Infectious vs. Time graph', type = 'l', col = 'red', lwd = 1)
```

![](Lab-4_files/figure-gfm/plot%20SI%20graphs-1.png)<!-- -->
