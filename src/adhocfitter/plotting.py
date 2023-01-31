import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

import exoplanet as xo

from . import timeseries as aftimeseries



def plot_model_light_curve(
    ldlc_obj, orbit, radii_planets,
    texp, supersampling_factor,
    max_phases, num_points=2000):
    model_phases = []
    model_light_curve = []
    for i, max_phase in enumerate(max_phases):
        model_phase = np.linspace(-max_phase, max_phase, num_points)
        model_phases.append(model_phase)
        model_time = model_phase * orbit.period[i] + orbit.t0[i]
        model_light_curve.append(ldlc_obj.get_light_curve(
            orbit=orbit,
            r=radii_planets,
            t=model_time,
            texp=texp,
            oversample=supersampling_factor,
        ).eval()[:, i])
    return model_phases, model_light_curve


def plot_multi_planet_folded_light_curve(
    num_planets, orbit, rp,
    lc_times, lc_dfluxes, filters, limb_dark_params, mean_fluxes,
    exposure_times, supersampling_factors,
    nbins, max_phases, alphas, dpi=300,
    detrend_series=None,
    detrend_coeffs=None):

    num_light_curves = len(lc_times)
    if detrend_series is None:
        detrend_design = [None] * num_light_curves
    else:
        detrend_design = []
        for detrend in detrend_series:
            if detrend is None:
                detrend_design.append(None)
                continue
            detrend_design.append(np.hstack((
                np.ones((len(detrend), 1)),
                detrend)))
    if detrend_coeffs is None:
        detrend_coeffs = [None] * num_light_curves

    periods = orbit.period.eval()
    epochs = orbit.t0.eval()
    radii_planets = rp * orbit.r_star.eval()
    ldlc_objs = dict()
    for filter_name in set(filters):
        ldlc_objs[filter_name] = xo.LimbDarkLightCurve(
                limb_dark_params[f'u_{filter_name}'])
    
    fig, axs = plt.subplots(
        num_light_curves,
        num_planets,
        sharex='col',
        sharey='row',
        figsize=(5*num_planets, 3*num_light_curves),
        dpi=dpi,
        squeeze=False,
    )

    for row, lc_time, lc_dflux, filter_name, texp, supersample, nbin, alpha, design, coeffs in zip(
        axs, lc_times, lc_dfluxes, filters,
        exposure_times, supersampling_factors,
        nbins, alphas, detrend_design, detrend_coeffs):

        ldlc_obj = ldlc_objs[filter_name]
        model_dflux = ldlc_obj.get_light_curve(
            orbit=orbit,
            r=radii_planets,
            t=lc_time,
            texp=texp,
            oversample=supersample,
        ).eval()
        model_dflux_sum = np.sum(model_dflux, axis=1)
        plot_model_phases, plot_model_dflux = plot_model_light_curve(
            ldlc_obj, orbit, radii_planets, texp, supersample, max_phases)
        
        if design is None:
                lc_trend = 0.
        else:
            if coeffs is not None:
                lc_trend = design @ coeffs
            else:
                raise ValueError(f'Unknown detrending coeffs for light curve')

        for i, (ax, period, epoch, max_phase) in enumerate(
            zip(row, periods, epochs, max_phases)):

            lc_phase = aftimeseries.phase_fold(lc_time, epoch, period)
            lc_dflux_only_planet = lc_dflux - model_dflux_sum + model_dflux[:, i] - mean_fluxes[i] - lc_trend
            ax.plot(
                lc_phase,
                lc_dflux_only_planet,
                '.', color='gray', alpha=alpha, rasterized=True)
            ax.plot(
                plot_model_phases[i], plot_model_dflux[i])

            phase_mask = np.abs(lc_phase) < max_phase
            select_phase = lc_phase[phase_mask]
            select_dflux = lc_dflux_only_planet[phase_mask]
            if nbin is not None:
                binned_mean, bins, _ = stats.binned_statistic(
                    select_phase, select_dflux, statistic='mean', bins=nbin)
                binned_err, _, _ = stats.binned_statistic(
                    select_phase, select_dflux, statistic=lambda a: np.std(a, ddof=1), bins=bins)
                binned_count, _ = np.histogram(select_phase, bins=bins)
                print(binned_count)
                mid_bin = (bins[:-1] + np.diff(bins) / 2.)
                ax.errorbar(mid_bin, binned_mean,
                    yerr=binned_err/np.sqrt(binned_count),
                    fmt='.', color='C1')
            ax.set_xlim(-max_phase, max_phase)
    fig.tight_layout()
    return fig, axs



rv_data_default_style = {
    'ecolor': 'gray',
    'elinewidth': 1,
    'alpha': 0.7,
    'fmt': 'o',
    'markersize': 3,
}
rv_model_default_style = {
    'color': 'slateblue',
    'zorder': 0,
}
rv_trend_unc_default_style = {
    'color': 'slateblue',
    'alpha': 0.2,
}


def plot_model_rv(num_planets, orbit, rv_semiamps, num_points=500):
    periods = orbit.period.eval()
    epochs = orbit.period.eval()
    model_phase = np.linspace(-0.5, 0.5, num_points)
    model_rvs = []
    for i in range(num_planets):
        model_time = model_phase * orbit.period[i] + orbit.t0[i]
        model_rvs.append(orbit.get_radial_velocity(
            model_time, K=rv_semiamps,
        ).eval()[:, i])
    return model_phase, model_rvs


def plot_multi_planet_folded_rv(
    folded_axs,
    unfolded_ax,
    residual_ax,
    num_planets, orbit,
    rv_semiamps, gammas, jitters,
    rv_times, rv_data, rv_uncs, rv_names,
    trends=None, model_trend_func=None,
    model_trend_unc_func=None,
    rv_data_style=rv_data_default_style,
    rv_model_style=rv_model_default_style,
    rv_inst_styles=None,
    rv_trend_unc_style=rv_trend_unc_default_style,
):
    if rv_inst_styles is None:
        rv_inst_styles = [dict()] * len(rv_times)

    min_obs_time = np.min([np.min(t) for t in rv_times])
    max_obs_time = np.max([np.max(t) for t in rv_times])
    unfold_model_time = np.linspace(min_obs_time, max_obs_time, num=int((max_obs_time-min_obs_time)//0.05))
    unfold_model = np.sum(orbit.get_radial_velocity(
            unfold_model_time, K=rv_semiamps,
    ).eval(), axis=1)
    if unfolded_ax is not None:
        if model_trend_func is not None:
            unfold_model_trend = model_trend_func(unfold_model_time)
            unfold_model += unfold_model_trend
            unfolded_ax.plot(unfold_model_time, unfold_model_trend, color='gray', linestyle='--')
            if model_trend_unc_func is not None:
                unfold_model_trend_unc = model_trend_unc_func(unfold_model_time)
                trend_unc = unfold_model_trend + unfold_model_trend_unc
                unfolded_ax.fill_between(unfold_model_time, trend_unc[0], trend_unc[1], **rv_trend_unc_style)
        unfolded_ax.plot(unfold_model_time, unfold_model, **rv_model_style)
        if residual_ax is not None:
            residual_ax.axhline(0, **rv_model_style)
            if model_trend_unc_func is not None:
                residual_ax.fill_between(unfold_model_time, unfold_model_trend_unc[0], unfold_model_trend_unc[1], **rv_trend_unc_style)

    model_phase, model_rvs = plot_model_rv(num_planets, orbit, rv_semiamps)
    for i, (ax, model_rv) in enumerate(zip(folded_axs, model_rvs)):
        ax.plot(model_phase, model_rv, **rv_model_style)

    epochs = orbit.t0.eval()
    periods = orbit.period.eval()
    if trends is None:
        trends = [0.] * len(rv_times)
    rv_data_style = rv_data_style.copy()
    for rv_time, rv, rv_unc, gamma, jitter, trend, label, plot_style in zip(
        rv_times, rv_data, rv_uncs, gammas, jitters, trends, rv_names, rv_inst_styles):
        rv_errorbar = np.sqrt(rv_unc**2+jitter**2)
        rv_shifted = rv - gamma
        rv_detrend = rv_shifted - trend

        model_rv = orbit.get_radial_velocity(
            rv_time, K=rv_semiamps,
        ).eval()
        model_rv_sum = np.sum(model_rv, axis=1)

        rv_data_style.update(plot_style)
        if unfolded_ax is not None:
            unfolded_ax.errorbar(rv_time, rv_shifted, rv_errorbar, label=label, **rv_data_style)
            if residual_ax is not None:
                residual_ax.errorbar(rv_time, rv_detrend-model_rv_sum, rv_errorbar, label=label, **rv_data_style)

        for i, (ax, period, epoch) in enumerate(zip(folded_axs, periods, epochs)):
            rv_phase = aftimeseries.phase_fold(rv_time, epoch, period)
            rv_only_planet = rv_shifted - model_rv_sum + model_rv[:, i] - trend
            ax.errorbar(rv_phase, rv_only_planet, rv_errorbar, label=label, **rv_data_style)


def make_multi_planet_rv_axes(num_planets, unfolded=True, residuals=True, figure_kwargs={'dpi': 300}):
    figure_kwargs['figsize'] = (7, 8)
    heights = [5, 3]
    fig = plt.figure(constrained_layout=False, **figure_kwargs)
    if unfolded:
        gs = gridspec.GridSpec(2, num_planets, figure=fig, wspace=0.25, hspace=0.25, height_ratios=heights)
        folded_row = 1
        if residuals:
            gs0 = gs[0, :].subgridspec(2, 1, hspace=0, height_ratios=[2, 1])
            unfolded_ax = fig.add_subplot(gs0[0])
            residual_ax = fig.add_subplot(gs0[1], sharex=unfolded_ax)
        else:
            unfolded_ax = fig.add_subplot(gs[0, :])
            residual_ax = None
    else:
        folded_row = 0
        unfolded_ax = None
        residual_ax = None
    folded_axs = []
    for i in range(num_planets):
        folded_axs.append(fig.add_subplot(gs[folded_row, i]))
    return fig, folded_axs, unfolded_ax, residual_ax
