# sims_featureScheduler_runs1.4
The continuing excitement of more simulated LSST surveys

Updated in this release:

* A more realistic database of the expected seeing conditions. This means the median atmospheric seeing is increased from ~0.6" to 0.7". Seasonal seeing variation should be more realistic as well
* We restrict the camera rotation to prevent tracking beyond rotator limits.
* Refactored the configuration scripts to include more documentation and better readability


## AGN_DDF

Scheduling deep drilling fields so they are optimized for AGN studies.

## DCR

This includes taking u and g filter observations at high airmass so DCR can be measured.

## DESC_DDF

Scheduling deep drilling fields so they are optimized for the Dark Energy Science Collaboration.


## alt_roll_dust

Simulations using a footprint that masks out high extinction regions. The sims use a rolling cadence that emphasizes half the sky, as well we a version that alternates observing north and south each day.

## baseline

The baseline survey that all the other experiments use as a template starting point. We include baseline simulations with visits with 1x30s exposures and 2x15s exposures.

Unlike previous baseline simulations, we now take more observations in pairs (previously u and y filter observations were unpaired).

Note the baseline survey with 2 snaps does not quite read the SRD requirements for number of visits in the WFD area. 

## euclid_DDF

Scheduling deep drilling fields so they include the planned Euclid fields

## footprints

Running the baseline strategy with different survey footprints.

## pair_strat

We experiment with different strategies of how we take observations in pairs. 

## rolling

Testing rolling cadences where we divide the WFD area into 2, 3, or 6 bands. These all perform regular full-sky coverage at the start and end of the survey.

The rolling cadences has a larger mean slewtime, even though the number of fitler changes seems about the same. The basis function weights probably need more tuning for rolling cadence strategy, or the rolling footprint basis function might not be quite right since the slewtime problem appears in the first year.

## short_exp

In addition to the usual baseline, we cover the survey footprint with shorter exposures.

## spiders

We set the camera rotator to +/- 45 degrees so diffraction spikes lie along rows and columns

## twilight_filters

We test limiting the number of filters available in/near twilight

## twilight_neo

Executing a NEO search strategy in twilight time

## var_expt

Uses variable exposure times so observations have similar limiting depths.

## weather

Mostly illustrative tests where we look at performance with different levels of weather downtime and maintenance downtime

## wfd_depth

Mostly illustrative to see how performance varies as we change the requested WFD depth.

