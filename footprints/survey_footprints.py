import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
import lsst.sims.featureScheduler.utils as utils
from lsst.sims.featureScheduler.utils import generate_goal_map
from lsst.sims.featureScheduler.utils import standard_goals

# OK, what are the footprints we'd like to try?

def bluer_footprint(nside=32):
    """Try a bluer filter dist. May want to turn this into a larger parameter search.
    """
    result = {}
    result['u'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.31,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    # Turn up the g in WFD
    result['g'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.2,
                                    WFD_fraction=0.9,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    result['r'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    result['i'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    # Turn down the z and y in WFD
    result['z'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.4,
                                    WFD_fraction=0.7,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    result['y'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.7,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15)
    return result


def gp_smooth(nside=32):
    # Treat the galactic plane as just part of the WFD
    result = {}
    result['u'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.31,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.31)
    result['g'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.2,
                                    WFD_fraction=0.44,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.44)
    result['r'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=1.0)
    result['i'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=1.0)
    result['z'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.4,
                                    WFD_fraction=0.9,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.9)
    result['y'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.9,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.9)
    return result


def no_gp_north(nside=32):
    result = {}
    gl1 = 290.
    gl2 = 40.
    result['u'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.31,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    result['g'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.2,
                                    WFD_fraction=0.44,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    result['r'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    result['i'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.46,
                                    WFD_fraction=1.0,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    result['z'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.4,
                                    WFD_fraction=0.9,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    result['y'] = generate_goal_map(nside=nside,
                                    NES_fraction=0.,
                                    WFD_fraction=0.9,
                                    SCP_fraction=0.15,
                                    GP_fraction=0.15,
                                    gp_long1=gl1, gp_long2=gl2)
    return result


def add_mag_clouds(inmap=None, nside=32):
    if inmap is None:
        result = standard_goals(nside=nside)
    else:
        result = inmap
    mag_clouds_hpix = utils.magellanic_clouds_healpixels(nside)
    for key in result:
        result[key][np.where(mag_clouds_hpix == 1)[0]] = np.max(result[key])
    return result


def big_sky(nside=32, weights={'u': [0.31, 0.15, False], 'g': [0.44, 0.15],
                               'r': [1., 0.3], 'i': [1., 0.3], 'z': [0.9, 0.3],
                               'y': [0.9, 0.3, False]}):
    """
    Based on the Olsen et al Cadence White Paper
    """
    wfd_north = 12.4
    wfd_south = -72.25
    gal_lat_limit = 15.
    full_north = 30.

    # WFD in big sky = dec range -72.5 to 12.5, avoiding galactic plane |b| < 15. deg.
    bigsky = utils.WFD_no_gp_healpixels(nside, dec_min=wfd_south, dec_max=wfd_north,
                                        center_width=gal_lat_limit, gal_long1=0, gal_long2=360)
    # Add extention to the north, up to 30 deg.
    ra, dec = utils.ra_dec_hp_map(nside=nside)
    bigsky = np.where((dec > np.radians(wfd_north)) & (dec < np.radians(full_north)), 1.e-6, bigsky)

    # Now let's break it down by filter
    result = {}
    for key in weights:
        result[key] = bigsky + 0.
        result[key][np.where(result[key] == 1)] = weights[key][0]
        result[key][np.where(result[key] == 1e-6)] = weights[key][1]
        if len(weights[key]) == 3:
            result[key][np.where(dec > np.radians(wfd_north))] = 0.
    return result


def big_sky_nouiy(nside=32, weights={'u': [0.31, 0., False], 'g': [0.44, 0.15],
                                     'r': [1., 0.3], 'i': [1., 0.0, False], 'z': [0.9, 0.3],
                                     'y': [0.9, 0.0, False]}):
    result = big_sky(nside=nside, weights=weights)
    return result


def big_sky_dust(nside=32, weights={'u': [0.31, 0.15, False], 'g': [0.44, 0.15],
                 'r': [1., 0.3], 'i': [1., 0.3], 'z': [0.9, 0.3],
                 'y': [0.9, 0.3, False]}, dust_limit=0.19):
    """
    Based on the Olsen et al Cadence White Paper
    """
    wfd_north = 12.4
    wfd_south = -72.25
    full_north = 30.

    # WFD in big sky = dec range -72.5 to 12.5, avoiding galactic plane |b| < 15. deg.
    bigsky = utils.WFD_no_dust_healpixels(nside, dec_min=wfd_south, dec_max=wfd_north,
                                          dust_limit=dust_limit)
    # Add extention to the north, up to 30 deg.
    ra, dec = utils.ra_dec_hp_map(nside=nside)
    bigsky = np.where((dec > np.radians(wfd_north)) & (dec < np.radians(full_north)), 1.e-6, bigsky)

    # Now let's break it down by filter
    result = {}
    for key in weights:
        result[key] = bigsky + 0.
        result[key][np.where(result[key] == 1)] = weights[key][0]
        result[key][np.where(result[key] == 1e-6)] = weights[key][1]
        if len(weights[key]) == 3:
            result[key][np.where(dec > wfd_north)] = 0.
    return result


def new_regions(nside=32, north_limit=2.25):
    ra, dec = utils.ra_dec_hp_map(nside=nside)
    coord = SkyCoord(ra=ra*u.rad, dec=dec*u.rad)
    g_long, g_lat = coord.galactic.l.radian, coord.galactic.b.radian

    # OK, let's just define the regions
    north = np.where((dec > np.radians(north_limit)) & (dec < np.radians(30.)))[0]
    wfd = np.where(utils.WFD_healpixels(dec_min=-72.25, dec_max=north_limit, nside=nside) > 0)[0]
    nes = np.where(utils.NES_healpixels(dec_min=north_limit, min_EB=-30., max_EB=10.) > 0)[0]
    scp = np.where(utils.SCP_healpixels(nside=nside, dec_max=-72.25) > 0)[0]

    new_gp = np.where((dec < np.radians(north_limit)) & (dec > np.radians(-72.25)) &
                      (np.abs(g_lat) < np.radians(15.)) &
                      ((g_long < np.radians(90.)) | (g_long > np.radians(360.-70.))))[0]

    anti_gp = np.where((dec < np.radians(north_limit)) & (dec > np.radians(-72.25)) &
                       (np.abs(g_lat) < np.radians(15.)) &
                       (g_long < np.radians(360.-70.)) & (g_long > np.radians(90.)))[0]

    footprints = {'north': north, 'wfd': wfd, 'nes': nes, 'scp': scp, 'gp': new_gp, 'gp_anti': anti_gp}

    return footprints


def newA(nside=32):
    """
    From https://github.com/rhiannonlynne/notebooks/blob/master/New%20Footprints.ipynb

    XXX--this seems to have very strange u-band distributions
    """
    zeros = np.zeros(hp.nside2npix(nside), dtype=float)
    footprints = new_regions(north_limit=12.25)

    # Define how many visits per field we want
    obs_region = {'gp': 750, 'wfd': 839, 'nes': 255, 'scp': 200, 'gp_anti': 825, 'north': 138}

    wfd_ratio = {'u': 0.06796116504854369, 'g': 0.0970873786407767,
                 'r': 0.22330097087378642, 'i': 0.22330097087378642, 'z': 0.1941747572815534, 'y': 0.1941747572815534}
    uniform_ratio = {'u': 0.16666666666666666, 'g': 0.16666666666666666,
                     'r': 0.16666666666666666, 'i': 0.16666666666666666, 'z': 0.16666666666666666, 'y': 0.16666666666666666}

    filter_ratios = {'gp': wfd_ratio,
                     'gp_anti': wfd_ratio,
                     'nes': {'u': 0.0, 'g': 0.2, 'r': 0.4, 'i': 0.4, 'z': 0.0, 'y': 0.0},
                     'north': uniform_ratio,
                     'scp': uniform_ratio,
                     'wfd': wfd_ratio}

    results = {}
    norm_val = obs_region['wfd']*filter_ratios['wfd']['r']
    for filtername in filter_ratios['wfd']:
        results[filtername] = zeros + 0
        for region in footprints:
            if np.max(filter_ratios[region][filtername]) > 0:
                results[filtername][footprints[region]] = filter_ratios[region][filtername]*obs_region[region]/norm_val

    return results


def newB(nside=32):
    """
    From https://github.com/rhiannonlynne/notebooks/blob/master/New%20Footprints.ipynb

    XXX--this seems to have very strange u-band distributions
    """
    zeros = np.zeros(hp.nside2npix(nside), dtype=float)
    footprints = new_regions(north_limit=12.25)

    # Define how many visits per field we want
    obs_region = {'gp': 650, 'wfd': 830, 'nes': 255, 'scp': 200, 'gp_anti': 100, 'north': 207}

    wfd_ratio = {'u': 0.06796116504854369, 'g': 0.0970873786407767,
                 'r': 0.22330097087378642, 'i': 0.22330097087378642, 'z': 0.1941747572815534, 'y': 0.1941747572815534}
    uniform_ratio = {'u': 0.16666666666666666, 'g': 0.16666666666666666,
                     'r': 0.16666666666666666, 'i': 0.16666666666666666, 'z': 0.16666666666666666, 'y': 0.16666666666666666}

    filter_ratios = {'gp': wfd_ratio,
                     'gp_anti': wfd_ratio,
                     'nes': {'u': 0.0, 'g': 0.2, 'r': 0.4, 'i': 0.4, 'z': 0.0, 'y': 0.0},
                     'north': uniform_ratio,
                     'scp': uniform_ratio,
                     'wfd': wfd_ratio}

    results = {}
    norm_val = obs_region['wfd']*filter_ratios['wfd']['r']
    for filtername in filter_ratios['wfd']:
        results[filtername] = zeros + 0
        for region in footprints:
            if np.max(filter_ratios[region][filtername]) > 0:
                results[filtername][footprints[region]] = filter_ratios[region][filtername]*obs_region[region]/norm_val

    return results


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


def stuck_rolling(nside=32, scale_down_factor=0.2):
    """A bit of a trolling footprint. See what happens if we use a rolling footprint, but don't roll it. 
    """
    sg = standard_goals()
    footprints = slice_wfd_area(2, sg, scale_down_factor=scale_down_factor)
    # Only take the first set
    footprints = footprints[0]
    return footprints







