import pandas as pd
import numpy as np

import third_party.keplersplinev2.keplersplinev2 as ksp


TESS_EPOCH = 2457000


def phase_fold(time_series, mid_transit, period):
    return ((time_series - mid_transit) / period + 0.5) % 1. - 0.5


def estimate_white_noise(light_curve, detrend_key, is_magnitude, planet_mask, cut=3.):
    mask = light_curve['mask']
    noise_mask = np.logical_and(mask, np.logical_not(planet_mask))
    if is_magnitude:
        mag = light_curve[detrend_key][noise_mask]
        zero_point, mean_std, _ = ksp.robust_mean(mag, cut=cut)
        unmasked_flux = 10 ** (-0.4 * (light_curve[detrend_key] - zero_point))
        flux = unmasked_flux[noise_mask]
    else:
        unmasked_flux = light_curve[detrend_key]
        flux = unmasked_flux[noise_mask]
        _, mean_std, _ = ksp.robust_mean(flux, cut=cut) 
    stddev = mean_std  * np.sqrt(len(flux) - 1)
    return stddev, unmasked_flux


def estimate_white_noise_v2(flux, planet_mask, cut=3.):
    masked_flux = flux[np.logical_not(planet_mask)]
    mean, mean_std, _ = ksp.robust_mean(masked_flux, cut)
    return mean_std * np.sqrt(len(masked_flux) + 1), mean_std


def mask_planet(time, period, mid_transit_time, transit_duration, duration_factor):
    """Only use in-transit points within `duration_factor` times of the transit duration."""
    phase_limit = transit_duration / period / 2. * duration_factor
    phase = (time - mid_transit_time) / period
    phase += 0.5
    phase %= 1.
    phase -= 0.5
    planet_mask = np.abs(phase) < phase_limit
    return planet_mask


def read_tess_lc(filename, periods, t0s, tdurs):
    lc_table = pd.read_csv(filename, index_col=None)
    cols = lc_table.columns
    light_curve = np.zeros(
        len(lc_table),
        dtype=[('mask', '?'), ('BJD', 'f8'), ('KSPSAP_FLUX', 'f8'), ('KSPSAP_FLUX_UNC', 'f8')])
    light_curve['mask'][:] = True
    light_curve['BJD'][:] = lc_table['time']
    light_curve['KSPSAP_FLUX'][:] = lc_table['kspflux']
    prelim_planet_mask = np.zeros(len(lc_table), dtype=bool)
    for period, t0, tdur in zip(periods, t0s, tdurs):
        prelim_planet_mask = np.logical_or(
            prelim_planet_mask, mask_planet(
                light_curve['BJD'], period, t0, tdur, 1.5))
    if 'kspflux_unc' in cols:
        light_curve['KSPSAP_FLUX_UNC'][:] = lc_table['kspflux_unc']
    else:
        tess_noise, _ = estimate_white_noise(light_curve, 'KSPSAP_FLUX', False, prelim_planet_mask, cut=3.)

    planet_mask = np.zeros(len(lc_table), dtype=bool)
    for period, t0, tdur in zip(periods, t0s, tdurs):
        planet_mask = np.logical_or(
            planet_mask, mask_planet(
                light_curve['BJD'], period, t0, tdur, 5.))
    tess_mask = np.logical_and(light_curve['mask'], planet_mask)

    tess_time = np.ascontiguousarray(light_curve['BJD'][tess_mask])
    tess_flux = np.ascontiguousarray(light_curve['KSPSAP_FLUX'][tess_mask])
    tess_dflux = tess_flux - 1
    if 'kspflux_unc' in cols:
        tess_noise = np.ascontiguousarray(light_curve['KSPSAP_FLUX_UNC'][tess_mask])
    return tess_time, tess_dflux, tess_noise


def read_generic_lc(filename, time_key, flux_key, unc_key,
                    detrend_keys=None, time_epoch=0., input_dflux=False,
                    **kwargs):
    lc_table = pd.read_csv(filename, **kwargs)
    glc_time = np.array(lc_table[time_key]) - (TESS_EPOCH - time_epoch)
    glc_dflux = np.array(lc_table[flux_key]) - (0 if input_dflux else 1)
    glc_unc = np.array(lc_table[unc_key])
    if detrend_keys is None:
        return glc_time, glc_dflux, glc_unc
    else:
        glc_detrend = np.array(lc_table.loc[:, detrend_keys])
        glc_detrend -= np.average(glc_detrend, axis=0)
        glc_detrend /= np.max(np.abs(glc_detrend), axis=0)
        return glc_time, glc_dflux, glc_unc, glc_detrend


def read_generic_rv(filename):
    df = pd.read_csv(
        filename,
        comment='#',
        )
    return df


def select_rv_by_instrument(table, instrument):
    select_result = table.query(f'Tel == "{instrument}" & f_RV == 0').sort_values('BJD')
    return (
        np.array(select_result['BJD']),
        np.array(select_result['RV']),
        np.array(select_result['e_RV']),
    )
