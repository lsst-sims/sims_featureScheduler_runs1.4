# sims_featureScheduler_runs1.4
The continuing excitement of more simulated LSST schedulers

Updated things in this release:

* A more realistic database of the expected seeing conditions. This mean the median atmospheric seeing is increased from ~0.6" to 0.7". Seasonal seeing variation should be more realistic as well
* We restrict the camera rotation slightly more to prevent tracking beyond rotator limits




## Note--these runs are still under development and subject to change. We will announce on https://community.lsst.org/c/sci when they are officially released.






## baseline

The baseline experiment. XXX--maybe a brief description of relevant points


## pair_strat

We experiment with different strategies of how we take observations in pairs. In the baseline, most observations are taken as g+r, r+i, i+z, while u, y, and sometimes z are taken unpaired. 

## satellite_dodge

There is the potential that LSST will want to avoid large satellite constellations. In this experiment, we randomly assign ~15% of observations as potentially including a satellite and pre-preemptively drop them from the observing queue. This is a test to see the impact of implamenting a simple satellite mitigation strategy.

## twilight_filters

We test limiting the number of filters available in/near twilight

## weather

Mostly illustrative tests where we look at performance with different levels of weather downtime and maintenance downtime


