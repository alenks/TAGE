# TAGE
Implementation of TAGE branch predictor in Sniper simulator is done as part of course assignment for CS5222: Advanced Computer Architecture. Branch predictor showed prediction accuracies ranging from 89.88% to 95.07% for different high MPKI benchmarks, whereas the prediction accuracies for other benchmarks showed ranges from 92.94% to 99.56% from SPEC CPU2006 benchmark suite. Need to make a new config file (or edit an existing one) to call the TAGE branch predictor.
```
[perf_model/branch_predictor]
type = tage
mispredict_penalty = 8 # Reflects just the front-end portion (approx) of the penalty for Interval Simulation
```
