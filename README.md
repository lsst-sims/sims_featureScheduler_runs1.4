# sims_featureScheduler_runs1.4
The continuing excitement of more simulated LSST surveys

Updated in this release:

* A more realistic database of the expected seeing conditions. This means the median atmospheric seeing is increased from ~0.6" to 0.7". Seasonal seeing variation should be more realistic as well
* We restrict the camera rotation to prevent tracking beyond rotator limits
* Refactored the configuration scripts to include more documentation and better readability
* Incorporated a new class to ensure cross-platform repeatability (always compare floats at a fixed precision rather than machine precision which can vary)
* Improved default DDF behavior (lower airmass observations)
* Refactored rolling cadence basis functions, resulting in more uniform final a survey
* Adjusted the default behavior for when filters are loaded into the filter changer


## AGN_DDF

Scheduling deep drilling fields so they are optimized for AGN studies.

## DCR

This includes taking u and g filter observations at high airmass so DCR can be measured.

## DDF_eperiment

We try some different DDF strategies, making longer DDF observing seasons, execuiting only in dark time, and leaving the u filter mounted for longer.

## DESC_DDF

Scheduling deep drilling fields so they are optimized for the Dark Energy Science Collaboration.

## alt_roll_dust

Simulations using a footprint that masks out high extinction regions. The sims use a rolling cadence that emphasizes half the sky, as well we a version that alternates observing north and south each day.

## baseline

The baseline survey that all the other experiments use as a template starting point. We include baseline simulations that have visits with 1x30s exposures and 2x15s exposures.

Unlike previous baseline simulations, we now take more observations in pairs (previously u and y filter observations were unpaired).

Note the baseline survey with 2 snaps does not quite read the SRD requirements for number of visits in the WFD area. 

## bulge

Experiments where we cover the bulge to different depths and use different cadences.

## euclid_DDF

Scheduling deep drilling fields so they include the planned Euclid fields.

## footprints

Running the baseline strategy with different survey footprints.

## pair_strat

We experiment with different strategies of how we take observations in pairs. 

## rolling

Testing rolling cadences where we divide the WFD area into 2, 3, or 6 bands. These all perform regular full-sky coverage at the start and end of the survey.

## short_exp

In addition to the usual baseline, we cover the survey footprint with shorter exposures.

## spiders

We set the camera rotator to +/- 45 degrees so diffraction spikes lie along rows and columns of the CCDs

## third_obs

We take some time at the end of the night to attempt to observe a third observations of regions that have already been observed twice.

## twilight_filters

We test limiting the number of filters available in/near twilight

## twilight_neo

Executing a NEO search strategy in twilight time.

## u60

Testing taking 60s u-filter exposures (to lift u observations above readnoise dominated regime)

## u_pairs

Testing different u-pairing strategies and varying how long the u filter is mounted.

## var_expt

Uses variable exposure times so observations have similar limiting depths.

## weather

Mostly illustrative tests where we look at performance with different levels of weather downtime and maintenance downtime.

## wfd_depth

Mostly illustrative to see how performance varies as we change the requested WFD depth.

## wfd_vary

More illustrative runs to show that adding observations to the galactic plane does not change the WFD depth much.

