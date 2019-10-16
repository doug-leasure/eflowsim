# eflowsim

Simulate population dynamics and survey data that respond to stream hydrology.  

The R package does the following steps:
1. Acquire hydrology time series data from a user-defined stream with a USGS stream gauge
2. Calculate four hydrologic metrics describing mean, variance, high and low flow characteristics.
3. Allow user to define correlations between hydrologic metrics and population parameters: carrying capacity and intrinsic growth rate.
4. Simulate a population time series in the selected stream using a multi-population viability analysis model (MPVA; Leasure et al 2019 Ecology).
5. Simulate field survey data from this population using the sampling and observation models from MPVA.

The user can try to recover the known flow-ecology relationships from the simulated survey data using any analysis method.  

Simulations can be repeated many times to calculate type I, II, and III error rates.