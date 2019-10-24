import numpy as np
import matplotlib.pylab as plt
import healpy as hp
from lsst.sims.featureScheduler.modelObservatory import Model_observatory
from lsst.sims.featureScheduler.schedulers import Core_scheduler
from lsst.sims.featureScheduler.utils import standard_goals, create_season_offset
import lsst.sims.featureScheduler.basis_functions as bf
from lsst.sims.featureScheduler.surveys import (generate_dd_surveys, Greedy_survey,
                                                Blob_survey)
from lsst.sims.featureScheduler import sim_runner
import lsst.sims.featureScheduler.detailers as detailers
import sys
import subprocess
import os
import argparse


def gen_greedy_surveys(nside, nexp=1, m5_weight=3., footprint_weight=0.3, slewtime_weight=3.,
                       stayfilter_weight=3., footprints=None, season_modulo=None, day_offset=None,
                       all_footprints_sum=None, all_rolling_sum=None):
    """
    Make a quick set of greedy surveys
    """

    # Let's remove the bluer filters since this should only be near twilight
    filters = ['r', 'i', 'z', 'y']
    surveys = []

    detailer = detailers.Camera_rot_detailer(min_rot=-87., max_rot=87.)

    for filtername in filters:
        bfs = []
        bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight))
        fps = [footprint[filtername] for footprint in footprints]
        bfs.append((bf.Footprint_rolling_basis_function(filtername=filtername,
                                                        footprints=fps,
                                                        out_of_bounds_val=np.nan, nside=nside,
                                                        all_footprints_sum=all_footprints_sum, all_rolling_sum=all_rolling_sum,
                                                        day_offset=day_offset, season_modulo=season_modulo), footprint_weight))
        bfs.append((bf.Slewtime_basis_function(filtername=filtername, nside=nside), slewtime_weight))
        bfs.append((bf.Strict_filter_basis_function(filtername=filtername), stayfilter_weight))
        # Masks, give these 0 weight
        bfs.append((bf.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.), 0))
        bfs.append((bf.Moon_avoidance_basis_function(nside=nside, moon_distance=30.), 0))

        bfs.append((bf.Filter_loaded_basis_function(filternames=filtername), 0))
        bfs.append((bf.Planet_mask_basis_function(nside=nside), 0))

        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        surveys.append(Greedy_survey(basis_functions, weights, block_size=1, filtername=filtername,
                                     dither=True, nside=nside, ignore_obs='DD', nexp=nexp,
                                     detailers=[detailer], survey_name='greedy'))

    return surveys


def generate_blobs(nside, mixed_pairs=False, nexp=1, offset=None,
                   m5_weight=6., footprint_weight=0.6, slewtime_weight=3.,
                   stayfilter_weight=3., template_weight=12., footprints=None,
                   season_modulo=None, day_offset=None,
                   all_footprints_sum=None, all_rolling_sum=None):

    # List to hold all the surveys (for easy plotting later)
    surveys = []

    # Set up observations to be taken in blocks
    filter1s = ['u', 'g', 'r', 'i', 'z', 'y']
    if mixed_pairs:
        filter2s = [None, 'r', 'i', 'z', None, None]
    else:
        filter2s = [None, 'g', 'r', 'i', None, None]

    # Ideal time between taking pairs
    pair_time = 22.
    times_needed = [pair_time, pair_time*2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(detailers.Camera_rot_detailer(min_rot=-87., max_rot=87.))
        detailer_list.append(detailers.Close_alt_detailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight/2.))
            bfs.append((bf.M5_diff_basis_function(filtername=filtername2, nside=nside), m5_weight/2.))

        else:
            bfs.append((bf.M5_diff_basis_function(filtername=filtername, nside=nside), m5_weight))

        if filtername2 is not None:
            fps = [footprint[filtername] for footprint in footprints]
            bfs.append((bf.Footprint_rolling_basis_function(filtername=filtername,
                                                            footprints=fps,
                                                            out_of_bounds_val=np.nan, nside=nside,
                                                            all_footprints_sum=all_footprints_sum,
                                                            all_rolling_sum=all_rolling_sum,
                                                            season_modulo=season_modulo,
                                                            day_offset=day_offset), footprint_weight/2.))
            fps = [footprint[filtername2] for footprint in footprints]
            bfs.append((bf.Footprint_rolling_basis_function(filtername=filtername2,
                                                            footprints=fps,
                                                            out_of_bounds_val=np.nan, nside=nside,
                                                            all_footprints_sum=all_footprints_sum,
                                                            all_rolling_sum=all_rolling_sum,
                                                            season_modulo=season_modulo,
                                                            day_offset=day_offset), footprint_weight/2.))
        else:
            fps = [footprint[filtername] for footprint in footprints]
            bfs.append((bf.Footprint_rolling_basis_function(filtername=filtername,
                                                            footprints=fps,
                                                            out_of_bounds_val=np.nan, nside=nside,
                                                            all_footprints_sum=all_footprints_sum,
                                                            all_rolling_sum=all_rolling_sum,
                                                            season_modulo=season_modulo,
                                                            day_offset=day_offset), footprint_weight))

        bfs.append((bf.Slewtime_basis_function(filtername=filtername, nside=nside), slewtime_weight))
        bfs.append((bf.Strict_filter_basis_function(filtername=filtername), stayfilter_weight))

        if filtername2 is not None:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints[-1][filtername],
                                                         n_obs=3, season=300.), template_weight/2.))
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername2, nside=nside,
                                                         footprint=footprints[-1][filtername2],
                                                         n_obs=3, season=300.), template_weight/2.))
        else:
            bfs.append((bf.N_obs_per_year_basis_function(filtername=filtername, nside=nside,
                                                         footprint=footprints[-1][filtername],
                                                         n_obs=3, season=300.), template_weight))
        # Masks, give these 0 weight
        bfs.append((bf.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.), 0.))
        bfs.append((bf.Moon_avoidance_basis_function(nside=nside, moon_distance=30.), 0.))
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.Filter_loaded_basis_function(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.Time_to_twilight_basis_function(time_needed=time_needed), 0.))
        bfs.append((bf.Not_twilight_basis_function(), 0.))
        bfs.append((bf.Planet_mask_basis_function(nside=nside), 0.))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = 'blob, %s' % filtername
        else:
            survey_name = 'blob, %s%s' % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.Take_as_pairs_detailer(filtername=filtername2))
        surveys.append(Blob_survey(basis_functions, weights, filtername1=filtername, filtername2=filtername2,
                                   ideal_pair_time=pair_time, nside=nside,
                                   survey_note=survey_name, ignore_obs='DD', dither=True,
                                   nexp=nexp, detailers=detailer_list))

    return surveys


def run_sched(surveys, survey_length=365.25, nside=32, fileroot='baseline_', verbose=False,
              extra_info=None):
    years = np.round(survey_length/365.25)
    scheduler = Core_scheduler(surveys, nside=nside)
    n_visit_limit = None
    observatory = Model_observatory(nside=nside)
    observatory, scheduler, observations = sim_runner(observatory, scheduler,
                                                      survey_length=survey_length,
                                                      filename=fileroot+'%iyrs.db' % years,
                                                      delete_past=True, n_visit_limit=n_visit_limit,
                                                      verbose=verbose, extra_info=extra_info)


def slice_wfd_area(nslice, target_map, scale_down_factor=0.2):
    """
    Slice the WFD area into even dec bands
    """
    # Make it so things still sum to one.
    scale_up_factor = nslice - scale_down_factor*(nslice-1)

    wfd = target_map['r'] * 0
    wfd_indices = np.where(target_map['r'] == 1)[0]
    wfd[wfd_indices] = 1
    wfd_accum = np.cumsum(wfd)
    split_wfd_indices = np.floor(np.max(wfd_accum)/nslice*(np.arange(nslice)+1)).astype(int)
    split_wfd_indices = split_wfd_indices.tolist()
    split_wfd_indices = [0] + split_wfd_indices

    all_scaled_down = {}
    for filtername in target_map:
        all_scaled_down[filtername] = target_map[filtername]+0
        all_scaled_down[filtername][wfd_indices] *= scale_down_factor

    scaled_maps = []
    for i in range(len(split_wfd_indices)-1):
        new_map = {}
        indices = wfd_indices[split_wfd_indices[i]:split_wfd_indices[i+1]]
        for filtername in all_scaled_down:
            new_map[filtername] = all_scaled_down[filtername] + 0
            new_map[filtername][indices] = target_map[filtername][indices]*scale_up_factor
        scaled_maps.append(new_map)

    return scaled_maps


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.set_defaults(pairs=True)
    parser.add_argument("--verbose", dest='verbose', action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default=365.25*10)
    parser.add_argument("--outDir", type=str, default="")
    parser.add_argument("--maxDither", type=float, default=0.7, help="Dither size for DDFs (deg)")
    parser.add_argument("--splits", type=int, default=2)
    parser.add_argument("--scale_down_factor", type=float, default=0.2)

    args = parser.parse_args()
    survey_length = args.survey_length  # Days
    outDir = args.outDir
    verbose = args.verbose
    max_dither = args.maxDither
    mod_year = args.splits
    scale_down_factor = args.scale_down_factor

    nside = 32
    per_night = True
    nexp = 1

    extra_info = {}
    exec_command = ''
    for arg in sys.argv:
        exec_command += ' ' + arg
    extra_info['exec command'] = exec_command
    extra_info['git hash'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
    extra_info['file executed'] = os.path.realpath(__file__)

    fileroot = 'testrolling_'
    file_end = 'v1.4_'

    observatory = Model_observatory(nside=nside)
    conditions = observatory.return_conditions()

    # Mark position of the sun at the start of the survey. Usefull for rolling cadence.
    sun_ra_0 = conditions.sunRA  # radians
    offset = create_season_offset(nside, sun_ra_0) + 365.25
    # Set up the DDF surveys to dither
    dither_detailer = detailers.Dither_detailer(per_night=per_night, max_dither=max_dither)
    details = [detailers.Camera_rot_detailer(min_rot=-87., max_rot=87.), dither_detailer]
    ddfs = generate_dd_surveys(nside=nside, nexp=nexp, detailers=details)

    sg = standard_goals()
    roll_maps = slice_wfd_area(mod_year, sg, scale_down_factor=scale_down_factor)
    footprints = roll_maps + [sg]

    all_footprints_sum = 0
    all_rolling_sum = 0

    wfd_indx = np.where(sg['r'] == 1)
    for fp in sg:
        all_footprints_sum += np.sum(sg[fp])
        all_rolling_sum += np.sum(sg[fp][wfd_indx])

    greedy = gen_greedy_surveys(nside, nexp=nexp, footprints=footprints, season_modulo=mod_year, day_offset=offset,
                                all_footprints_sum=all_footprints_sum, all_rolling_sum=all_rolling_sum)
    blobs = generate_blobs(nside, nexp=nexp, mixed_pairs=True, day_offset=offset, footprints=footprints, season_modulo=mod_year,
                           all_footprints_sum=all_footprints_sum, all_rolling_sum=all_rolling_sum)
    surveys = [ddfs, blobs, greedy]
    run_sched(surveys, survey_length=survey_length, verbose=verbose,
              fileroot=os.path.join(outDir, fileroot+file_end), extra_info=extra_info,
              nside=nside)
