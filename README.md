# TAGE
Implementation of TAGE branch predictor in Sniper simulator is done as part of course assignment for CS5222: Advanced Computer Architecture. Branch predictor showed prediction accuracies ranging from 91% to 99% for different benchmarks. Need to make a new config file (or edit an existing one) to call the TAGE branch predictor.
```[perf_model/branch_predictor]
type = tage
mispredict_penalty = 8 # Reflects just the front-end portion (approx) of the penalty for Interval Simulation```
