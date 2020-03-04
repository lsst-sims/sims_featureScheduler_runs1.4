import numpy as np
from lsst.sims.featureScheduler.surveys import BaseSurvey
import copy
import lsst.sims.featureScheduler.basis_functions as basis_functions
from lsst.sims.featureScheduler.utils import empty_observation
from lsst.sims.featureScheduler import features
import logging
import random
from lsst.sims.featureScheduler.surveys import Deep_drilling_survey


def dd_bfs(RA, dec, survey_name, ha_limits, frac_total=0.0185/2., aggressive_frac=0.011/2.):
    """
    Convienence function to generate all the feasibility basis functions
    """
    sun_alt_limit = -18.
    time_needed = 62.
    fractions = [0.00, aggressive_frac, frac_total]
    bfs = []
    bfs.append(basis_functions.Not_twilight_basis_function(sun_alt_limit=sun_alt_limit))
    bfs.append(basis_functions.Time_to_twilight_basis_function(time_needed=time_needed))
    bfs.append(basis_functions.Hour_Angle_limit_basis_function(RA=RA, ha_limits=ha_limits))
    bfs.append(basis_functions.Moon_down_basis_function())
    bfs.append(basis_functions.Fraction_of_obs_basis_function(frac_total=frac_total, survey_name=survey_name))
    bfs.append(basis_functions.Look_ahead_ddf_basis_function(frac_total, aggressive_frac,
                                                             sun_alt_limit=sun_alt_limit, time_needed=time_needed,
                                                             RA=RA, survey_name=survey_name,
                                                             ha_limits=ha_limits))
    bfs.append(basis_functions.Soft_delay_basis_function(fractions=fractions, delays=[0., 0.5, 1.5],
                                                         survey_name=survey_name))

    return bfs


def generate_dd_surveys(nside=None, nexp=2, detailers=None, reward_value=100):
    """Utility to return a list of standard deep drilling field surveys.

    XXX-Someone double check that I got the coordinates right!

    """

    surveys = []

    # ELAIS S1
    RA = 9.45
    dec = -44.
    survey_name = 'DD:ELAISS1'
    ha_limits = ([0., 3], [20, 24.])
    bfs = dd_bfs(RA, dec, survey_name, ha_limits)
    surveys.append(Deep_drilling_survey(bfs, RA, dec, sequence='urgizy',
                                        nvis=[8, 20, 10, 20, 26, 20],
                                        survey_name=survey_name, reward_value=reward_value,
                                        nside=nside, nexp=nexp, detailers=detailers))

    # XMM-LSS
    survey_name = 'DD:XMM-LSS'
    RA = 35.708333
    dec = -4-45/60.
    ha_limits = ([0., 3], [20, 24.])
    bfs = dd_bfs(RA, dec, survey_name, ha_limits)

    surveys.append(Deep_drilling_survey(bfs, RA, dec, sequence='urgizy',
                                        nvis=[8, 20, 10, 20, 26, 20], survey_name=survey_name, reward_value=reward_value,
                                        nside=nside, nexp=nexp, detailers=detailers))

    # Extended Chandra Deep Field South
    RA = 53.125
    dec = -28.-6/60.
    survey_name = 'DD:ECDFS'
    ha_limits = [[0.5, 3.0], [20., 22.5]]
    bfs = dd_bfs(RA, dec, survey_name, ha_limits)
    surveys.append(Deep_drilling_survey(bfs, RA, dec, sequence='urgizy',
                                        nvis=[8, 20, 10, 20, 26, 20],
                                        survey_name=survey_name, reward_value=reward_value, nside=nside,
                                        nexp=nexp, detailers=detailers))

    # COSMOS
    RA = 150.1
    dec = 2.+10./60.+55/3600.
    survey_name = 'DD:COSMOS'
    ha_limits = ([0., 3], [20, 24.])
    bfs = dd_bfs(RA, dec, survey_name, ha_limits)
    surveys.append(Deep_drilling_survey(bfs, RA, dec, sequence='urgizy',
                                        nvis=[8, 20, 10, 20, 26, 20],
                                        survey_name=survey_name, reward_value=reward_value, nside=nside,
                                        nexp=nexp, detailers=detailers))

    # Euclid Fields
    # I can use the sequence kwarg to do two positions per sequence
    filters = 'urgizy'
    nexps = [8, 5, 7, 19, 24, 5]
    survey_name = 'DD:EDFS'
    # Note the sequences need to be in radians since they are using observation objects directly
    RAs = np.radians([58.97, 63.6])
    decs = np.radians([-49.28, -47.60])
    sequence = []
    exptime = 30
    visit_nexp = 1
    for filtername, nexp in zip(filters, nexps):
        for ra, dec in zip(RAs, decs):
            for num in range(nexp):
                obs = empty_observation()
                obs['filter'] = filtername
                obs['exptime'] = exptime
                obs['RA'] = ra
                obs['dec'] = dec
                obs['nexp'] = visit_nexp
                obs['note'] = survey_name
                sequence.append(obs)

    ha_limits = ([0., 3], [20, 24.])
    # And back to degrees for the basis function
    bfs = dd_bfs(np.degrees(RAs[0]), np.degrees(decs[0]), survey_name, ha_limits)
    surveys.append(Deep_drilling_survey(bfs, RA, dec, sequence=sequence,
                                        survey_name=survey_name, reward_value=reward_value, nside=nside,
                                        nexp=nexp, detailers=detailers))

    return surveys
