import os
from pathlib import Path

import corner
import gc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import seaborn as sns
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from dynesty import plotting as dyplot
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

from .formatter import FormatParams
from .model import calc_astro
from .sort_nicely import sort_nicely as sn
from ..lib import util

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction": "in", "ytick.direction": "in"})
matplotlib.use('Agg')
matplotlib.rcParams.update({'lines.markeredgewidth': 0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset': False})


# #00
def mjd_to_utc(time):
    """Converts a list of MJDs to a list of dates in years."""
    res = []
    for i in time:
        res.append(Time(i, format='mjd').jyear)
    return np.array(res)


def utc_to_mjd(time):
    """Converts a list of dates in years to a list of MJDs."""
    res = []
    for i in time:
        res.append(Time(i, format='jyear').mjd)
    return np.array(res)


def mjd_to_isot(time):
    """Converts a list of MJDs to a list of dates in YYYY-MM-DD."""
    res = []
    for i in time:
        res.append(str(Time(i, format='mjd').isot).split('T')[0])
    return np.array(res)


def obs_times(meta, times, ivisits, iorbits, updated=False):
    """Plot of the visit index as a function of observed time for the
    observations. Includes a table with the number of orbits in each
    visit and a zoom into each visit.

    Parameters
    ----------
    updated : bool
        If the user decided to not use all visits but set some "which_visits"
        in the pcf, this bool is need to save a plot for all files in the data
        directory and a plot for the onces defined with "which_visits".
        It prevents that when the function is being called again, the previous
        plot isnt overwritten.
    """
    fig = plt.figure(figsize=(12, 5))

    Nvisits = np.max(ivisits)+1  # Number of visits
    gs = GridSpec(nrows=Nvisits, ncols=2)
    ax1 = fig.add_subplot(gs[:, 0])

    visit_start = []
    for i in range(Nvisits):
        axn = fig.add_subplot(gs[i, 1])
        if i == 0:
            axn.set_title('each visit')
        axn.text(0.95, 0.41, f'v{i}',
                 horizontalalignment='center',
                 verticalalignment='center', transform=axn.transAxes)
        visit_mask = ivisits == i
        times_i = times[visit_mask]
        visit_start.append(times_i[0])
        axn.scatter(times_i, np.ones(len(times_i)),
                    marker='s', c='r', s=5, alpha=0.5)
        axn.get_xaxis().set_ticks([])
        axn.get_yaxis().set_ticks([])
        axn.axes.xaxis.set_visible(False)
        axn.axes.yaxis.set_visible(False)

    ax1.scatter(times, ivisits, marker='s', c='r', s=5, alpha=0.5)
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_ylabel('visit number')
    ax1.set_xlabel('time (MJD)')
    ax1.set_ylim(np.max(ivisits)+0.9, np.min(ivisits)-0.9)  # Turn around y axis
    ax1sec = ax1.secondary_xaxis('top', functions=(mjd_to_utc, utc_to_mjd))
    ax1sec.set_xlabel('time (year)')

    # It checks the orbit index right before it goes back to 0 (because a new visit starts).
    # We also have to append the very last orbit number
    norbits = iorbits[np.where(np.diff(iorbits) < 0)]
    norbits = np.append(norbits, iorbits[-1:])+1
    table_data = np.array([np.arange(len(norbits), dtype=int), norbits, mjd_to_isot(visit_start)])
    col_names = ['ivisit', '#orbits', 'start times']
    table = ax1.table(cellText=table_data.T,
                      cellLoc='center', colLabels=col_names,
                      loc='right', bbox=[-0.7, 0., 0.55, 1])
    table.scale(0.5, 3.0)
    plt.subplots_adjust(wspace=0.04)

    if meta.save_obs_times_plot:
        obs_date_dir = meta.workdir / 'figs' / 's00_obs_dates'
        if not obs_date_dir.exists():
            obs_date_dir.mkdir(parents=True)

        if not updated:
            plt.savefig(obs_date_dir / 'obs_dates_all.png',
                        bbox_inches='tight', pad_inches=0.05, dpi=120)
        else:
            plt.savefig(obs_date_dir / 'obs_dates.png',
                        bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


# 02
def barycorr(x, y, z, time, obsx, obsy, obsz, coordtable, meta):
    """This function plots the vectorfile positions of HST and where
    the observations where taken

    Parameters
    ----------
    x: array
        X position from vectorfile.
    y: array
        Y position from vectorfile.
    z: array
        Z position from vectorfile.
    time: array
        times from the vectorfile.
    obsx: array
        X position of observations.
    obsy: array
        Y position of observations.
    obsz: array
        Z position of observations.
    coordtable
        a list of files containing the vector information of HST
        downloaded in s01.
    meta
        the name of the metadata file.

    Returns
    -------
    Saves and/or Shows a plot.

    Notes
    -----
    History:
        Written by Sebastian Zieba      December 2021
    """
    plt.rcParams["figure.figsize"] = (8, 6)
    fig = plt.figure(1001)
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(x, y, z, c=time, s=10, cmap=cm.inferno)
    ax.scatter(x[0], y[0], z[0], s=200,
               marker='x', c='r', label='horizons start')
    ax.scatter(x[-1], y[-1], z[-1], s=200,
               marker='x', c='b', label='horizons end')
    ax.scatter(obsx, obsy, obsz, s=50, marker='x', c='k', label='observations')
    # [ax.text(x[i],y[i], z[i],'{0}'.format(astropy.time.Time(val=time, format='jd', scale='utc').iso[i])) for i in range(len(x))[::10]]
    # ax.text(x[0],y[0], z[0],'{0}'.format(astropy.time.Time(val=time[0], format='jd', scale='utc').iso))
    # ax.text(x[-1],y[-1], z[-1],'{0}'.format(astropy.time.Time(val=time[-1], format='jd', scale='utc').iso))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    cbar = fig.colorbar(cax, orientation='vertical', label='Time (JD)')
    plt.legend()
    plt.tight_layout()
    if meta.save_barycorr_plot:
        barycorr_dir = meta.workdir / 'figs' / 's02_barycorr'
        if not barycorr_dir.exists():
            barycorr_dir.mkdir(parents=True)
        plt.savefig(barycorr_dir / f'bjdcorr_{coordtable.split('/')[-1].split('.')[0]}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


# 03
def smooth(meta, x, y, x_smoothed, y_smoothed):
    """Plots the raw stellar spectrum and the smoothed spectrum."""
    plt.plot(x, y, label='raw spectrum')
    plt.plot(x_smoothed, y_smoothed, label='smoothed spectrum')
    plt.legend(loc=1)
    plt.xscale('log')
    plt.xlabel('wavelength (m)')
    plt.ylabel('normalized intensity')
    if meta.save_smooth_plot:
        smooth_dir = meta.workdir / 'figs' / 's03_smooth'
        if not smooth_dir.exists():
            smooth_dir.mkdir(parents=True)
        plt.savefig(smooth_dir / 'smooth.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def refspec(bp_wvl, bp_val, sm_wvl, sm_flux, ref_wvl, ref_flux, meta):
    """Plots the bandpass, the stellar spectrum and the product of the both."""
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.figure(1003)
    plt.clf()
    plt.plot(bp_wvl, bp_val, label='bandpass')
    plt.plot(sm_wvl, sm_flux, label='stellar spectrum ({0})'.format(meta.sm))
    plt.plot(ref_wvl, ref_flux, label='stellar spectrum * bandpass')
    plt.xscale('log')
    if meta.grism == 'G102':
        plt.xlim(0.6e-6, 1.4e-6)
    elif meta.grism == 'G141':
        plt.xlim(0.8e-6, 2e-6)
    plt.xlabel('wavelength (m)')
    plt.ylabel('norm. intensity')
    plt.legend(loc=4)
    plt.tight_layout()
    if meta.save_refspec_plot:
        refspec_dir = meta.workdir / 'figs' / 's03_refspec'
        if not refspec_dir.exists():
            refspec_dir.mkdir(parents=True)
        plt.savefig(refspec_dir / 'refspec.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


# 10
def image_quick(ima, i, meta):
    """This plots the full direct image."""
    plt.figure(10044)
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    plt.rcParams['image.cmap'] = 'viridis'
    nrow = len(ima[1].data[:, 0])
    ncol = len(ima[1].data[0, :])

    plt.suptitle(f"Direct image #{i}, visit #{meta.ivisit_di[i]}, orbit #{meta.iorbit_di[i]}", y=0.94)

    ax.title.set_text('Full Direct Image')
    im = ax.imshow(ima[1].data * ima[0].header['exptime'], origin='lower', vmin=0, vmax=500)

    cmin = meta.di_cmin
    cmax = meta.di_cmax
    rmin = meta.di_rmin
    rmax = meta.di_rmax
    ax.plot([cmin, cmin, cmax, cmax, cmin],
            [rmin, rmax, rmax, rmin, rmin],
            lw=2, c='r', alpha=0.8)

    ax.set_xlabel('columns')
    ax.set_ylabel('rows')

    # make colorbar same size as subplot
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.colorbar(im, cax=cax, label='flux')

    plt.tight_layout()
    if meta.save_image_plot:
        s10_images_dir = meta.workdir / 'figs' / 's10_images'
        if not s10_images_dir.exists():
            s10_images_dir.mkdir(parents=True)
        plt.savefig(s10_images_dir / f'quick_di{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=180)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def image(dat, ima, results, i, meta):
    """This plots the full direct image with the guess of the target (defined
    using di_rmin, etc.) marked as a red box. It also plots a zoom into the
    guess position of the target with the gaussian fit solution marked with
    a cross.
    """
    plt.figure(1004)
    plt.clf()
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))

    plt.rcParams['image.cmap'] = 'viridis'

    nrow = len(ima[1].data[:, 0])
    ncol = len(ima[1].data[0, :])

    plt.suptitle(f"Direct image #{i}, visit #{meta.ivisit_di[i]}, orbit #{meta.iorbit_di[i]}", y=0.885)

    ax[0].title.set_text('Full Direct Image')
    ax[0].imshow(ima[1].data * ima[0].header['exptime'],
                 origin='lower', vmin=0, vmax=500)
    cmin = meta.di_cmin
    cmax = meta.di_cmax
    rmin = meta.di_rmin
    rmax = meta.di_rmax
    ax[0].plot([cmin, cmin, cmax, cmax, cmin],
               [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.8)
    ax[0].set_xlabel('columns')
    ax[0].set_ylabel('rows')

    ax[1].title.set_text('Cutout')
    image_cutout = dat * ima[0].header['exptime']
    vmax = np.max(image_cutout)
    im = ax[1].imshow(image_cutout, origin='lower', vmin=0, vmax=vmax/10)
    ax[1].plot(results[2], results[3], marker='x',
               color='orange', markeredgewidth=3.,
               ms=10, label='centroid', linestyle="none")
    ax[1].legend(numpoints=1)

    # make colorbar same size as subplot
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax[1])
    width = axes_size.AxesY(ax[1], aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.colorbar(im, cax=cax, label='flux')

    plt.tight_layout()
    if meta.save_image_plot:
        s10_images_dir = meta.workdir / 'figs' / 's10_images'
        if not s10_images_dir.exists():
            s10_images_dir.mkdir(parents=True)
        plt.savefig(s10_images_dir / f'quick_di{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


# 20
def sp2d(d, meta, i):
    """Plot the spectrum with a low vmax to make the background better visible."""
    plt.imshow(d[1].data, origin='lower', vmin=0, vmax=500)
    plt.colorbar()
    plt.tight_layout()
    plt.title(f'Spectrum w/ low vmax, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    if meta.save_sp2d_plot:
        s20_sp2d_dir = meta.workdir / 'figs' / 's20_sp2d'
        if not s20_sp2d_dir.exists():
            s20_sp2d_dir.mkdir(parents=True)
        plt.savefig(s20_sp2d_dir / f'sp2d_{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def badmask_2d(array1, array2, array3, meta, i):
    """Plots the badmask arrays which are used by the optimal
    extraction routine."""
    fig, ax = plt.subplots(1, 4, figsize=(9, 6), sharey=True, sharex=True)
    ax[0].set_title('bpix==4')
    ax[1].set_title('bpix==512')
    ax[2].set_title('flatfield')
    ax[3].set_title('combined')
    ax[0].imshow(array1, cmap='Greys')
    ax[1].imshow(array2, cmap='Greys')
    ax[2].imshow(array3, cmap='Greys')
    ax[3].imshow(array1 | array2 | array3, cmap='Greys')
    ax[0].set_xlabel('ncol')
    ax[0].set_ylabel('nrow')
    s20_badmask_dir = meta.workdir + 'figs' / 's20_badmask'
    if not s20_badmask_dir.exists():
        s20_badmask_dir.mkdir(parents=True)
    plt.savefig(s20_badmask_dir / f'badmask_{i}.png',
                bbox_inches='tight', pad_inches=0.05, dpi=120)
    plt.close('all')
    plt.clf()
    gc.collect()


# TODO: Q: The plot shouldnt change between observations in the same orbit (same direct image)
def trace(d, meta, visnum, orbnum, i):
    """Plots the spectrum together with the trace."""
    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry141 as geo

    trace = geo.trace(meta.refpix[:, 1], meta.refpix[:, 2])                #determines trace coefficients
    trace_i = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_i + meta.LTV1)          #start of trace
    trace_f = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_f + meta.LTV1)          #end of trace
    tracex = np.linspace(trace_i, trace_f, 100)                    #x values over which to compute trace
    tracey = meta.refpix[orbnum, 1] + meta.LTV2 + trace[0][orbnum] + \
            trace[1][orbnum]*(tracex - tracex[0])                    #y values of trace
    plt.imshow(d[1].data, origin='lower', vmin=0, vmax=4000)            #plots raw image
    plt.colorbar()
    plt.plot(tracex, tracey, color='yellow', linewidth=2)                #plots trace on raw image frame
    plt.title(f'Spectrum w/ trace (line), visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    if meta.save_trace_plot:
        s20_trace_dir = meta.workdir + 'figs' / 's20_trace'
        if not s20_trace_dir.exists():
            s20_trace_dir.mkdir(parents=True)
        plt.savefig(s20_trace_dir / f'trace_{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def bkg_hist(fullframe_diff, skymedian, meta, i, ii):
    """Plot saving a histogram of the fluxes in the up-the-ramp sample.
    Showing the user decided background threshold and the median flux below the threshold.
    """
    histo = fullframe_diff.flatten()
    fig, ax = plt.subplots(3, 1, figsize=(10, 8))
    ax[0].hist(histo, int(len(histo)/2000), facecolor='k', alpha=0.2)
    plt.suptitle(f'UpTheRamp {i}-{ii}, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    ax[0].axvline(meta.background_thld, lw=2, ls='-', c='b', label='background_thld')
    ax[0].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Flux')
    ax[0].set_ylabel('log. Frequ.')
    ax[0].legend(loc=1)
    ax[1].hist(histo, int(len(histo) / 10000), facecolor='k', alpha=0.2)
    ax[1].axvline(meta.background_thld, lw=2, ls='-', c='b', label='background_thld')
    ax[1].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian')
    ax[1].set_xlabel('Flux')
    ax[1].set_ylabel('lin. Frequ.')
    ax[1].legend(loc=1)
    ax[2].axvline(skymedian, lw=2, ls='--', c='b', label=f'skymedian {skymedian:.3g}')
    ax[2].set_xlabel('Flux')
    ax[2].set_ylabel('lin. Frequ.')
    zoom = 150
    # yyy, xxx, _ = ax[1].hist(histo, 51, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    ax[2].hist(histo, 41, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    # xmax = xxx[np.argsort(yyy)[::-1]][0]
    ax[2].legend(loc=1)
    plt.tight_layout()
    # data_new = histo[((skymedian - zoom) < histo) & (histo < (skymedian + zoom))]
    # var = util.median_abs_dev(data_new)
    # from scipy.stats import norm
    # Fit a normal distribution to the data:
    # mu, std = norm.fit(data_new)
    # x_new = np.linspace(skymedian - zoom, skymedian + zoom, 100)
    # p = norm.pdf(x_new, mu, std)
    # ax[1].plot(x_new, p*max(yyy)*10, 'k', linewidth=2)
    # ax[1].axvline(x_new[np.argmax(p)], lw=2, ls='--', c='g', label='gaussian {0:.3g}'.format(x_new[np.argmax(p)]))

    if meta.save_bkg_hist_plot:
        s20_bkg_hist_dir = meta.workdir / 'figs' / 's20_bkg_hist'
        if not s20_bkg_hist_dir.exists():
            s20_bkg_hist_dir.mkdir(parents=True)
        plt.savefig(s20_bkg_hist_dir / f'bkg_hist{i}-{ii}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=100)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def utr(diff, meta, i, ii, orbnum, rowmedian, rowmedian_absder, peaks):
    """
    Saves a plot of up-the-ramp sample, the row by row sum and the derivate of the latter. It furthermore shows the aperture used for the analysis.
    """
    cmin = int(meta.refpix[
                   orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_i + meta.LTV1  # determines left column for extraction (beginning of the trace)
    cmax = min(int(meta.refpix[orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_f + meta.LTV1,
               meta.subarray_size-5)  # right column (end of trace, or edge of detector)

    fig, ax = plt.subplots(1,3, figsize=(10,8), sharey=True)

    im = ax[0].imshow(diff, vmin=0, vmax=300, origin='lower')
    # ax[2].axhline(idx, c = 'b', ls='--')
    # ax[2].plot([cmin-cmin, cmin-cmin, cmax-cmin, cmax-cmin, cmin-cmin], [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.85)
    ax[0].set_xlim(cmin-cmin, cmax-cmin)
    ax[0].set_ylim(meta.rmin, meta.rmax)
    ax[0].axhline(min(peaks), c='r', ls='--', lw=2)
    ax[0].axhline(max(peaks), c='r', ls='--', lw=2)
    ax[0].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3)
    ax[0].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3)
    ax[0].title.set_text('spectrum, cmin to cmax')

    p1 = np.arange(len(diff))
    ax[1].axhline(min(peaks), c='r', ls='--', lw=2)
    ax[1].axhline(max(peaks), c='r', ls='--', lw=2)
    ax[1].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3)
    ax[1].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3)
    rowmm = np.median(rowmedian)
    peaks_mid = int(np.mean(peaks))
    # print(peaks_mid)
    ax[1].plot(rowmedian, p1)
    ax[1].axvline(rowmm, ls='-', c='k', alpha=0.4, lw=3)

    def peak_to_median_ratio(peak_val, median_val, percentage):
        return (peak_val-median_val)*(1-percentage)+median_val

    ax[1].axvline(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.95),
                  ls='--', c='k', alpha=0.4, lw=3)
    ax[1].axvline(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.99),
                  ls='--', c='k', alpha=0.4, lw=3)
    ax[1].axvline(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.995),
                  ls='--', c='k', alpha=0.4, lw=3)
    ax[1].axvline(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.999),
                  ls='--', c='k', alpha=0.4, lw=3)
    axtf = ax[1].get_xaxis_transform()
    ax[1].text(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.95),
               .20, '95%', transform=axtf)
    ax[1].text(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.99),
               .15, '99%', transform=axtf)
    ax[1].text(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.995),
               .1, '99.5%', transform=axtf)
    ax[1].text(peak_to_median_ratio(rowmedian[peaks_mid], rowmm, 0.999),
               .05, '99.9%', transform=axtf)
    ax[1].set_xscale('log')
    ax[1].set_ylim(meta.rmin, meta.rmax)
    ax[1].set_xlabel('log Flux')
    ax[1].title.set_text('median row Flux')

    p2 = (p1[1:] + p1[:-1]) / 2
    ax[2].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3)
    ax[2].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3)
    ax[2].scatter(rowmedian_absder[peaks], peaks, marker='x', c='r')
    ax[2].plot(rowmedian_absder, p2)
    ax[2].set_xscale('log')
    ax[2].set_ylim(meta.rmin, meta.rmax)
    ax[2].set_xlabel('log Flux')
    ax[2].title.set_text('Derivative')

    plt.colorbar(im)
    plt.tight_layout()
    plt.suptitle(f'UpTheRamp {i}-{ii}, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}', y=1)
    if meta.save_utr_plot:
        s20_utr_dir = meta.workdir / 'figs' / 's20_utr'
        if not s20_utr_dir.exists():
            s20_utr_dir.mkdir(parents=True)
        plt.savefig(s20_utr_dir / f'utr{i}-{ii}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=90)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def sp1d(template_waves, spec_box, meta, i, spec_opt=False):
    """Plots the resulting spectrum. If the user did optimal extraction,
    a comparison between optextr and box sum will be shown.
    """
    plt.clf()
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 6.4))
    plt.plot(template_waves, spec_box, c='r', label='box sum')
    if meta.opt_extract:
        plt.plot(template_waves, spec_opt, c='k', label='optimal extraction')
    plt.legend()
    plt.title(f'sp1d, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    plt.tight_layout()
    if meta.save_sp1d_plot:
        s29_sp1d_dir = meta.workdir / 'figs' / 's29_sp1d'
        if not s29_sp1d_dir.exists():
            s29_sp1d_dir.mkdir(parents=True)
        plt.savefig(s29_sp1d_dir / f'sp1d_{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def bkg_evo(bkg_evo, meta):
    """Plot of the background flux as a function of up the ramp sample."""
    plt.clf()
    fig, ax = plt.subplots(1, 1)
    ax.plot(range(len(bkg_evo)), bkg_evo)
    ax.set_xlabel('# diff image')
    ax.set_ylabel('diff flux')
    # plt.legend()
    plt.tight_layout()
    if meta.save_bkg_evo_plot:
        s20_bkg_evo_dir = meta.workdir / 'figs' / 's20_bkg_evo'
        if not s20_bkg_evo_dir.exists():
            s20_bkg_evo_dir.mkdir(parents=True)
        plt.savefig(s20_bkg_evo_dir / f'bkg_evo.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def sp1d_diff(sp1d_all_diff, meta, wvl_hires):
    """Difference of 1D spectrum between two consecutive exposures."""
    s20_sp1d_diff_dir = meta.workdir / 'figs' / 's20_sp1d_diff'
    if not s20_sp1d_diff_dir.exists():
        s20_sp1d_diff_dir.mkdir(parents=True)

    ylimmin = np.nanmin(sp1d_all_diff)
    ylimmax = np.nanmax(sp1d_all_diff)
    for iiii in range(len(sp1d_all_diff)):
        plt.clf()
        fig, ax = plt.subplots(1, 1, figsize=(6.4, 6.4*1.5/2))
        ax.plot(wvl_hires, sp1d_all_diff[iiii])
        if meta.grism == 'G141':
            ax.set_xlim(10000, 17800)
        elif meta.grism == 'G102':
            ax.set_xlim(5000, 16000)
        ax.set_ylim(ylimmin*1.2, ylimmax*1.2)
        fig.suptitle(f'Background Diff, from v{meta.ivisit_sp[iiii]} o{meta.iorbit_sp[iiii]} to v{meta.ivisit_sp[iiii+1]} o{meta.iorbit_sp[iiii+1]}',
                     fontsize=12, y=0.99)
        plt.tight_layout()
        if meta.save_sp1d_diff_plot:
            plt.savefig(s20_sp1d_diff_dir / f'sp1d_diff_{iiii}.png',
                        bbox_inches='tight', pad_inches=0.05, dpi=120)
            plt.close('all')
            plt.clf()
            gc.collect()
        else:
            plt.show()
            plt.close('all')
            plt.clf()
            gc.collect()


def utr_aper_evo(peaks_all, meta):
    """Plot of the evolution in aperture size."""
    peaks = np.array(peaks_all)
    peaks = np.sort(peaks)
    # nsamps_mask = []
    # for i in nsamps:
    #    nsamps_mask.append(np.arange(i, dtype=int))
    # nsamps_mask = np.array(nsamps_mask)
    # nsamps_mask[::2] = nsamps_mask[::-1][::2]
    # nsamps_mask = np.array(nsamps_mask).flatten()

    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), sharex=True)

    # for ii in range(max(nsamps)):
    #    ax[0].plot(range(len(peaks[:,0][nsamps_mask==ii])), peaks[:,0][nsamps_mask==ii])
    #    ax[0].plot(range(len(peaks[:,1][nsamps_mask==ii])), peaks[:,1][nsamps_mask==ii])
    #    ax[1].plot(range(len(peaks[:,0][nsamps_mask==ii])), np.diff(peaks[nsamps_mask==ii]))

    # ax[0].plot(range(len(peaks[:,0])), peaks[:,0])
    # ax[0].plot(range(len(peaks[:,1])), peaks[:,1])
    ax.plot(range(len(peaks[:, 0])), np.diff(peaks))

    ax.set_xlabel('# diff image')
    # ax[0].set_ylabel('peak position')
    ax.set_ylabel('peak distance')

    # plt.legend()
    plt.tight_layout()
    if meta.save_utr_aper_evo_plot:
        s20_utr_aper_evo_dir = meta.workdir / 'figs' / 's20_utr_aper_evo'
        if not s20_utr_aper_evo_dir.exists():
            s20_utr_aper_evo_dir.mkdir(parents=True)
        plt.savefig(s20_utr_aper_evo_dir / f'utr_aper_evo.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def refspec_fit(modelx, modely, p0, datax, datay, leastsq_res, meta, i):
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    plt.suptitle(f'refspec_fit {i}, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    ax.plot(modelx, modely, label='refspec')
    # ax.plot((p0[0] + p0[1] * datax), datay * p0[2], label='initial guess')
    ax.plot((leastsq_res[0] + datax * leastsq_res[1]), datay * leastsq_res[2],
             label=f'spectrum fit, wvl = ({leastsq_res[0]:.5g}+{leastsq_res[1]:.5g}*x), {leastsq_res[2]:.5g}')
    ax.plot(datax, datay/max(datay), label='spectrum before fit')
    if meta.grism == 'G141':
        ax.set_xlim(9800, 18000)
    elif meta.grism == 'G102':
        ax.set_xlim(7000, 12500)
    # ax.set_xscale('log')
    ax.set_xlabel('Wavelength (angstrom)')
    ax.set_xlabel('rel. Flux')
    plt.legend()
    plt.tight_layout()
    if meta.save_refspec_fit_plot:
        s20_refspec_fit_dir = meta.workdir / 'figs' / 's20_refspec_fit'
        if not s20_refspec_fit_dir.exists():
            s20_refspec_fit_dir.mkdir(parents=True)
        plt.savefig(s20_refspec_fit_dir / f'/refspec_fit_{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def refspec_fit_lin(modelx, modely, p0, datax, datay, leastsq_res, meta, i):
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    plt.suptitle(f'refspec_fit {i}, visit {meta.ivisit_sp[i]}, orbit {meta.iorbit_sp[i]}')
    ax.plot(modelx, modely, label='refspec')
    # ax.plot((p0[0] + p0[1] * datax), datay * p0[2], label='initial guess')
    ax.plot((leastsq_res[0] + datax), datay * leastsq_res[1],
             label=f'spectrum fit, wvl = ({leastsq_res[0]:.5g}+x), {leastsq_res[1]:.5g}')
    ax.plot(datax, datay/max(datay), label='spectrum before fit')
    if meta.grism == 'G141':
        ax.set_xlim(9000, 20000)
    elif meta.grism == 'G102':
        ax.set_xlim(5000, 16000)
    ax.set_xscale('log')
    ax.set_xlabel('Wavelength (angstrom)')
    ax.set_xlabel('rel. Flux')
    plt.legend()
    plt.tight_layout()
    if meta.save_refspec_fit_plot:
        s20_refspec_fit_dir = meta.workdir / 'figs' / 's20_refspec_fit'
        if not s20_refspec_fit_dir.exists():
            s20_refspec_fit_dir.mkdir(parents=True)
        plt.savefig(s20_refspec_fit_dir / f'refspec_fit_{i}.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def drift(leastsq_res_all, meta):
    a = np.array([i[0] for i in leastsq_res_all])
    b = np.array([i[1] for i in leastsq_res_all])
    c = np.array([i[2] for i in leastsq_res_all])

    fig, ax = plt.subplots(3, 1, figsize=(8,8), sharex=True)
    ax[0].set_title('Wavelength calibration: Function: c*(a+b*x)')
    ax[0].plot(range(len(a)), a, label='lin. wvl shift')
    ax[1].plot(range(len(b)), b, label='wvl scaling')
    ax[2].plot(range(len(c)), c, label='flux scaling')

    [ax[0].axvline(i, ls='--', c='k', lw=2.2, alpha=0.55) for i in meta.new_visit_idx_sp]
    [ax[1].axvline(i, ls='--', c='k', lw=2.2, alpha=0.55) for i in meta.new_visit_idx_sp]
    [ax[2].axvline(i, ls='--', c='k', lw=2.2, alpha=0.55) for i in meta.new_visit_idx_sp]

    [ax[0].axvline(i, ls='--', c='k', lw=2, alpha=0.35) for i in meta.new_orbit_idx_sp]
    [ax[1].axvline(i, ls='--', c='k', lw=2, alpha=0.35) for i in meta.new_orbit_idx_sp]
    [ax[2].axvline(i, ls='--', c='k', lw=2, alpha=0.35) for i in meta.new_orbit_idx_sp]

    ax[0].set_ylabel('a')
    ax[1].set_ylabel('b')
    ax[2].set_ylabel('c')

    plt.subplots_adjust(hspace=0.05)

    if meta.save_drift_plot:
        s20_drift_dir = meta.workdir / 'figs' / 's20_drift'
        if not s20_drift_dir.exists():
            s20_drift_dir.mkdir(parents=True)
        plt.savefig(s20_drift_dir / 'drift.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


def drift_lin(leastsq_res_all, meta):
    a = np.array([i[0] for i in leastsq_res_all])
    c = np.array([i[1] for i in leastsq_res_all])

    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    ax[0].set_title('Wavelength calibration: Function: c*(a+b*x)')
    ax[0].plot(range(len(a)), a, label='lin. wvl shift')
    ax[1].plot(range(len(c)), c, label='flux scaling')

    [ax[0].axvline(i, ls='--', c='k', lw=2.2, alpha=0.55) for i in meta.new_visit_idx_sp]
    [ax[1].axvline(i, ls='--', c='k', lw=2.2, alpha=0.55) for i in meta.new_visit_idx_sp]

    [ax[0].axvline(i, ls='--', c='k', lw=2, alpha=0.35) for i in meta.new_orbit_idx_sp]
    [ax[1].axvline(i, ls='--', c='k', lw=2, alpha=0.35) for i in meta.new_orbit_idx_sp]

    ax[0].set_ylabel('a')
    ax[1].set_ylabel('c')

    plt.subplots_adjust(hspace=0.05)

    if meta.save_drift_plot:
        s20_drift_dir = meta.workdir / 'figs' / 's20_drift'
        if not s20_drift_dir.exists():
            s20_drift_dir.mkdir(parents=True)
        plt.savefig(s20_drift_dir / 'drift.png',
                    bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
        plt.clf()
        gc.collect()
    else:
        plt.show()
        plt.close('all')
        plt.clf()
        gc.collect()


# 21
def plot_wvl_bins(w_hires, f_interp, wave_bins, wvl_bins, dirname):
    """Plot of a 1D spectrum and the bins."""
    plt.plot(w_hires, f_interp)
    for wave in wave_bins:
        plt.axvline(wave, color='0.5')
    plt.ylabel("Photelectrons")
    plt.xlabel("Wavelength (angstroms)")
    plt.tight_layout()
    plt.savefig(dirname / f'spec_bins{wvl_bins}.png')
    # plt.show()


# 30
def plot_raw(data, meta):
    """Saves a plot with the raw light curve (which includes the systematics)."""
    # palette = sns.color_palette("husl", data.nvisit)
    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)
    fig, ax = plt.subplots(data.nvisit, 1, figsize=(6, 3.6*data.nvisit), sharex=True)
    if data.nvisit > 1:
        for i in range(data.nvisit):
            ind = data.vis_num == i
            # plt.subplot((data.nvisit) * 100 + 10 + i + 1)
            ax[i].plot(data.t_vis[ind] / 60., data.flux[ind], marker='o',
                       markersize=3.0, linestyle="none",
                       label="Visit {0}".format(i), color=palette[i])
            ax[i].set_xlim(((data.t_vis.min() - 0.02) / 60, (data.t_vis.max() + 0.05) / 60))
            ax[i].set_ylim((0.998 * data.flux.min(), 1.002 * data.flux.max()))
            ax[i].legend(loc=1)
            ax[i].set_ylabel("Flux (e-)")
    else:
        # plt.subplot((data.nvisit) * 100 + 10 + i + 1)
        ax.plot(data.t_vis / 60., data.flux, marker='o',
                markersize=3.0, linestyle="none",
                label="Visit 0", color=palette[0])
        ax.set_xlim(((data.t_vis.min() - 0.02) / 60, (data.t_vis.max() + 0.05) / 60))
        ax.set_ylim((0.998 * data.flux.min(), 1.002 * data.flux.max()))
        ax.legend(loc=1)
        ax.set_ylabel("Flux (e-)")
    plt.xlabel("Time after visit start (hours)")
    fig.suptitle('wvl = {0:0.3f} micron'.format(meta.wavelength), fontsize=15, y=0.998)
    plt.tight_layout()

    raw_lc_dir = meta.workdir / meta.fitdir / 'raw_lc'
    if not raw_lc_dir.exists():
        raw_lc_dir.mkdir(parents=True)
    plt.savefig(raw_lc_dir / f'raw_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    plt.close('all')
    plt.clf()
    gc.collect()


def save_plot_raw_data(data, meta):
    """Saves the data used for the raw light curve plot."""
    table = Table()
    if data.nvisit > 1:
        for i in range(data.nvisit):
            ind = data.vis_num == i
            table[f't_vis_v{i}'] = np.array(data.t_vis[ind], dtype=np.float64)
            table[f'flux_v{i}'] = np.array(data.flux[ind], dtype=np.float64)
    else:
        table['t_vis'] = np.array(data.t_vis, dtype=np.float64)
        table['flux'] = np.array(data.flux, dtype=np.float64)

    raw_lc_dir = meta.workdir / meta.fitdir / 'raw_lc'
    if not raw_lc_dir.exists():
        raw_lc_dir.mkdir(parents=True)
    ascii.write(table, raw_lc_dir / f'raw_lc_data_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt',
                format='rst', overwrite=True)


def rmsplot(model, data, meta, fitter=None):
    """Plot RMS vs. bin size looking for time-correlated noise
    Taken from POET: https://github.com/kevin218/POET/blob/master/code/lib/plots.py.
    """
    time = data.time
    residuals = model.norm_resid
    residuals = residuals[np.argsort(time)]

    rms, stderr, binsz = util.computeRMS(residuals, binstep=1)
    normfactor = 1e-6
    plt.rcParams.update({'legend.fontsize': 11})
    plt.figure(1111, figsize=(8, 6))
    plt.clf()
    plt.suptitle(f'Correlated Noise at wvl = {meta.wavelength:0.3f} micron', size=16)
    plt.loglog(binsz, rms/normfactor, color='black', lw=1.5, label='Fit RMS', zorder=3)    # our noise
    plt.loglog(binsz, stderr/normfactor, color='red', ls='-', lw=2, label='Std. Err.', zorder=1)  # expected noise
    plt.xlim(0.95, binsz[-1]*2)
    plt.ylim(stderr[-1]/normfactor/2., stderr[0]/normfactor*2.)
    plt.xlabel("Bin Size", fontsize=14)
    plt.ylabel("RMS (ppm)", fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.legend()

    dname = Path(f'{fitter}_res')
    fname = dname / f'corr_plot_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png'
    plt.savefig(meta.workdir / meta.fitdir / fname)
    plt.close('all')
    util.save_allandata(binsz, rms, stderr, meta, fitter=fitter)
    plt.clf()
    gc.collect()


def plot_fit_lc2(data, fit, meta, mcmc=False, nested=False):
    """Plots phase folded fit."""
    plt.clf()
    fig, ax = plt.subplots(2, 1)

    p = FormatParams(fit.params, data)  # FIXME
    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)

    # ind = model.phase > 0.5
    # model.phase[ind] -= 1.
    # calculate a range of times at higher resolution to make model look nice
    phase_hr = np.linspace(fit.phase.min() - 0.05, fit.phase.max() + 0.05, 1000)
    t_hr = phase_hr * p.per[0] + p.t0[0] + data.toffset

    # plot data
    # plot best fit model from first visit
    ax[0].plot(phase_hr, calc_astro(t_hr, fit.params, data, fit.myfuncs, 0))

    # plot systematics removed data
    for i in range(data.nvisit):
        ind = data.vis_num == i
        ax[0].plot(fit.phase[ind], fit.data_nosys[ind], color=palette[i], marker='o', markersize=3, linestyle="none")

    # add labels/set axes
    # xlo, xhi = np.min(model.phase)*0.9, np.max(model.phase)*1.1
    xlo, xhi = -0.1, 0.1
    ax[0].set_xlim(xlo, xhi)
    ax[0].set_ylabel("Relative Flux")

    # annotate plot with fit diagnostics
    # ax = plt.gca()
    ax[0].text(0.85, 0.29,
               r'$\chi^2_{\nu}$:    ' + '{0:0.2f}'.format(fit.chi2red) + '\n'
               + 'obs. rms:  ' + '{0:0.1f}'.format(fit.rms) + '\n'
               + 'exp. rms:  ' + '{0:0.1f}'.format(fit.rms_predicted),
               verticalalignment='top', horizontalalignment='left',
               transform=ax[0].transAxes, fontsize=12)

    # plot fit residuals
    ax[1].axhline(0, zorder=1, color='0.2', linestyle='dashed')

    for i in range(data.nvisit):
        ind = data.vis_num == i
        ax[1].plot(fit.phase[ind], 1.0e6 * fit.norm_resid[ind],
                   color=palette[i], marker='o', markersize=3,
                   linestyle="none")

    # add labels/set axes
    ax[1].set_xlim(xlo, xhi)
    ax[1].set_ylabel("Residuals (ppm)")
    ax[1].set_xlabel("Orbital phase")

    if mcmc:
        fig.suptitle(f'MCMC, {meta.wavelength:0.3f} micron', fontsize=15, y=0.998)
    elif nested:
        fig.suptitle(f'Nested Sampling, {meta.wavelength:0.3f} micron', fontsize=15, y=0.998)
    else:
        fig.suptitle(f'LSQ, {meta.wavelength:0.3f} micron', fontsize=15, y=0.998)

    plt.tight_layout()

    if mcmc:
        plt.savefig(meta.workdir / meta.fitdir / 'fit_lc' / f'mcmc_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    elif nested:
        plt.savefig(meta.workdir / meta.fitdir / 'fit_lc' / f'nested_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    else:
        plt.savefig(meta.workdir / meta.fitdir / 'fit_lc' / f'lsq_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    plt.close('all')
    plt.clf()
    gc.collect()


def plot_fit_lc3(data, fit, meta, mcmc=False):
    """Plots light curve without systematics model."""
    plt.clf()
    fig, ax = plt.subplots(1, 1)

    p = FormatParams(fit.params, data)  # FIXME
    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)

    time_model = np.linspace(data.time.min() - 0.05, data.time.max() + 0.05, 1000)
    flux_model = calc_astro(time_model, fit.params, data, fit.myfuncs, 0)
    ax.plot(time_model, flux_model)

    # plot systematics removed data
    ax.plot(data.time, fit.data_nosys, marker='o', markersize=3, linestyle="none")

    # add labels/set axes
    # xlo, xhi = np.min(model.phase)*0.9, np.max(model.phase)*1.1
    # xlo, xhi = -0.1, 0.1
    # ax[0].set_xlim(xlo, xhi)
    ax.set_ylabel("Relative Flux")

    # annotate plot with fit diagnostics
    # ax = plt.gca()
    ax.text(0.85, 0.29,
            r'$\chi^2_{\nu}$:    ' + '{0:0.2f}'.format(fit.chi2red) + '\n'
            + 'obs. rms:  ' + '{0:0d}'.format(int(fit.rms)) + '\n'
            + 'exp. rms:  ' + '{0:0d}'.format(int(fit.rms_predicted)),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=12)

    # add labels/set axes
    ax.set_xlabel("Time")

    if mcmc:
        fig.suptitle(f'MCMC, {meta.wavelength:0.3f} micron', fontsize=15, y=0.998)
    else:
        fig.suptitle(f'LSQ, {meta.wavelength:0.3f} micron', fontsize=15, y=0.998)

    plt.tight_layout()
    # plt.show()
    fit_lc_dir = meta.workdir + meta.fitdir + 'fit_lc'
    if not os.path.isdir(fit_lc_dir):
        os.makedirs(fit_lc_dir)
    if mcmc:
        plt.savefig(fit_lc_dir / f'newmcmc_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    else:
        plt.savefig(fit_lc_dir / f'newfit_lc_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png')
    # plt.waitforbuttonpress(0) # this will wait for indefinite time
    plt.close('all')
    plt.clf()
    gc.collect()


def save_astrolc_data(data, fit, meta):
    """Saves the data used to plot the astrophysical model (without the systematics)
    and the data (without the systematics) not phase folded.
    """
    table_model = Table()
    table_nosys = Table()

    p = FormatParams(fit.params, data)

    time_model = np.linspace(data.time.min() - p.per[0]/2, data.time.max() + p.per[0]/2, 1000)
    flux_model = calc_astro(time_model, fit.params, data, fit.myfuncs, 0)

    table_model['time_model'] = np.array(time_model, dtype=np.float64)
    table_model['flux_model'] = np.array(flux_model, dtype=np.float64)

    table_nosys['time_nosys'] = np.array(data.time, dtype=np.float64)
    table_nosys['flux_nosys'] = np.array(fit.data_nosys, dtype=np.float64)

    fit_lc_dir = meta.workdir + meta.fitdir + 'fit_lc'
    if not os.path.isdir(fit_lc_dir):
        os.makedirs(fit_lc_dir)
    ascii.write(table_model, fit_lc_dir / f'fit_lc_data_model_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt',
                format='rst', overwrite=True)
    ascii.write(table_nosys, fit_lc_dir / f'fit_lc_data_nosys_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt',
                format='rst', overwrite=True)


# def plot_fit_lc4(data, fit, meta, mcmc=False):
#     #datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
#     plt.clf()
#     fig, ax = plt.subplots(2,1)
#     #print(fit.params)
#     p = FormatParams(fit.params, data)  # FIXME
#     sns.set_palette("muted")
#     palette = sns.color_palette("muted", data.nvisit)
#
#     time_days = np.linspace(data.time.min() - 0.05, data.time.max() + 0.05, 1000)
#
#
#     # colors = ['blue', 'red', 'orange', 'gray']
#
#     # plot data
#     # plot best fit model from first visit
#
#     ax[0].plot(data.time, calc_sys(data.time, fit.params, data, fit.myfuncs, 0))
#     ax[0].set_xlim(data.time[0], data.time[29])
#
#     # plot systematics removed data
#     #ax[0].plot(data.time, fit.data_nosys, marker='o', markersize=3, linestyle="none")
#
#
#     # add labels/set axes
#     # xlo, xhi = np.min(model.phase)*0.9, np.max(model.phase)*1.1
#     #xlo, xhi = -0.1, 0.1
#     #ax[0].set_xlim(xlo, xhi)
#     ax[0].set_ylabel("Relative Flux")
#
#     # annotate plot with fit diagnostics
#     #ax = plt.gca()
#     ax[0].text(0.85, 0.29,
#             '$\chi^2_\\nu$:    ' + '{0:0.2f}'.format(fit.chi2red) + '\n'
#             + 'obs. rms:  ' + '{0:0d}'.format(int(fit.rms)) + '\n'
#             + 'exp. rms:  ' + '{0:0d}'.format(int(fit.rms_predicted)),
#             verticalalignment='top', horizontalalignment='left',
#             transform=ax[0].transAxes, fontsize=12
#             )
#
#     # plot fit residuals
#     ax[1].axhline(0, zorder=1, color='0.2', linestyle='dashed')
#
#     # for i in range(data.nvisit):
#     #     ind = data.vis_num == i
#     #     ax[1].plot(fit.phase[ind], 1.0e6 * fit.norm_resid[ind], color=palette[i], marker='o', markersize=3,
#     #              linestyle="none")
#
#     # add labels/set axes
#     #ax[1].set_xlim(xlo, xhi)
#     ax[1].set_ylabel("Residuals (ppm)")
#     ax[1].set_xlabel("Time")
#
#     fig.suptitle('Filename: {0}'.format(meta.run_file.split('/')[-1].split('.txt')[0]), fontsize=15, y=0.998)
#
#     plt.tight_layout()
#     #plt.show()
#     if not os.path.isdir(meta.workdir + meta.fitdir + '/fit_lc'):
#         os.makedirs(meta.workdir + meta.fitdir + '/fit_lc')
#     if mcmc:
#         plt.savefig(meta.workdir + meta.fitdir + '/fit_lc' + "/new2mcmc_lc_{0}.png".format(meta.s30_file_counter))
#     else:
#         plt.savefig(meta.workdir + meta.fitdir + '/fit_lc' + "/new2fit_lc_{0}.png".format(meta.s30_file_counter))
#     # plt.waitforbuttonpress(0) # this will wait for indefinite time
#     #plt.close()


def params_vs_wvl(vals, errs, idxs, meta):
    """Plots every fitted parameter as a function of bin. It is able to show how astrophysical
    & systematical parameters change over wavelength.
    """
    labels = meta.labels
    fig, ax = plt.subplots(len(idxs[0]), 1, figsize=(6.4, 25), sharex=True)

    for i in range(len(idxs[0])):
        if len(idxs[0]) == 1:
            ax.errorbar(range(len(idxs)), [vals[ii][idxs[0][i]] for ii in range(len(vals))],
                        yerr=[errs[ii][idxs[0][i]] for ii in range(len(errs))], fmt='.')
            ax.set_ylabel(labels[i])
        else:
            ax[i].errorbar(range(len(idxs)), [vals[ii][idxs[0][i]] for ii in range(len(vals))],
                           yerr=[errs[ii][idxs[0][i]] for ii in range(len(errs))], fmt='.')
            ax[i].set_ylabel(labels[i])
    plt.subplots_adjust(hspace=0.01)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(meta.workdir / meta.fitdir / 'lsq_res' / 'lsq_params_vs_wvl.png',
                dpi=450, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()


def params_vs_wvl_mcmc(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta):
    """Plots every fitted parameter as a function of bin. It is able to show
    how astrophysical & systematical parameters change over wavelength.
    """
    vals_mcmc = np.array(vals_mcmc)
    errs_lower_mcmc = np.array(errs_lower_mcmc)
    errs_upper_mcmc = np.array(errs_upper_mcmc)

    labels = meta.labels

    fig, ax = plt.subplots(len(vals_mcmc[0]), 1, figsize=(6.4, 20), sharex=True)
    for i in range(len(vals_mcmc.T)):
        if len(vals_mcmc[0]) == 1:
            ax.errorbar(meta.wavelength_list, vals_mcmc.T[i],
                        yerr=(errs_lower_mcmc.T[i], errs_upper_mcmc.T[i]), fmt='.')
            ax.set_ylabel(labels[i])
        else:
            ax[i].errorbar(meta.wavelength_list, vals_mcmc.T[i],
                           yerr=(errs_lower_mcmc.T[i], errs_upper_mcmc.T[i]), fmt='.')
            ax[i].set_ylabel(labels[i])
    plt.subplots_adjust(hspace=0.01)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(meta.workdir / meta.fitdir / 'mcmc_res' / 'mcmc_params_vs_wvl.png',
                dpi=500, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()


def params_vs_wvl_nested(vals_nested, errs_lower_nested, errs_upper_nested, meta):
    """Plots every fitted parameter as a function of bin. It is able to show
    how astrophysical & systematical parameters change over wavelength.
    """
    vals_nested = np.array(vals_nested)
    errs_lower_nested = np.array(errs_lower_nested)
    errs_upper_nested = np.array(errs_upper_nested)

    labels = meta.labels

    fig, ax = plt.subplots(len(vals_nested[0]), 1, figsize=(6.4, 20), sharex=True)
    for i in range(len(vals_nested.T)):
        if len(vals_nested[0]) == 1:
            ax.errorbar(meta.wavelength_list, vals_nested.T[i],
                        yerr=(errs_lower_nested.T[i], errs_upper_nested.T[i]), fmt='.')
            ax.set_ylabel(labels[i])
        else:
            ax[i].errorbar(meta.wavelength_list, vals_nested.T[i],
                           yerr=(errs_lower_nested.T[i], errs_upper_nested.T[i]), fmt='.')
            ax[i].set_ylabel(labels[i])
    plt.subplots_adjust(hspace=0.01)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(meta.workdir / meta.fitdir / 'nested_res' / 'nested_params_vs_wvl.png',
                dpi=500, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()


def lsq_rprs(vals, errs, idxs, meta):
    """Plots the spectrum (rprs vs wvl) as fitted by the least square routine."""
    rp_idx = np.where(np.array(meta.labels) == 'rp')[0][0]
    rprs_vals_lsq = [vals[ii][idxs[0][rp_idx]] for ii in range(len(vals))]
    rprs_errs_lsq = [errs[ii][idxs[0][rp_idx]] for ii in range(len(errs))]
    plt.errorbar(meta.wavelength_list, rprs_vals_lsq, yerr=rprs_errs_lsq, fmt='.', c='red')
    plt.xlabel('Wavelength (micron)')
    plt.ylabel('Transit Depth (ppm)')
    plt.savefig(meta.workdir / meta.fitdir / 'lsq_res' / 'lsq_rprs.png',
                dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()


def mcmc_chains(ndim, sampler, nburn, labels, meta):
    """Plots the temporal evolution of the MCMC chain."""
    plt.clf()
    fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, ndim))
    if ndim > 1:
        for i in range(0, ndim):
            axes[i].plot(sampler.chain[:, nburn:, i].T, alpha=0.4)
            # axes[i].yaxis.set_major_locator(MaxNLocator(5))
            axes[i].set_ylabel(labels[i])
    elif ndim == 1:
        axes.plot(sampler.chain[:, nburn:, 0].T, alpha=0.4)
        # axes.yaxis.set_major_locator(MaxNLocator(5))
        axes.set_ylabel(labels)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    if nburn == 0:
        fig.savefig(meta.workdir / meta.fitdir / 'mcmc_res' / f"mcmc_chains_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png")
    else:
        fig.savefig(meta.workdir / meta.fitdir / 'mcmc_res' / f"mcmc_chains_noburn_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png")
    plt.close('all')
    plt.clf()
    gc.collect()


# FIXME: make sure this works for cases when nvisit > 1
def mcmc_pairs(samples, params, meta, fit_par, data):
    """Plots a pairs plot of the MCMC."""
    labels = meta.labels
    fig = corner.corner(samples, labels=labels, show_titles=True,quantiles=[0.16, 0.5, 0.84],title_fmt='.4')
    figname = meta.workdir / meta.fitdir / 'mcmc_res' /\
            f"mcmc_pairs_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png"
    fig.savefig(figname)
    plt.close('all')
    plt.clf()
    gc.collect()

# FIXME: make sure this works for cases when nvisit > 1
def nested_pairs(samples, params, meta, fit_par, data):
    """Plots a pairs plot of the nested sampling."""
    labels = meta.labels
    fig = corner.corner(samples, labels=labels, show_titles=True,
                        quantiles=[0.16, 0.5, 0.84], title_fmt='.4')
    figname = meta.workdir / meta.fitdir / 'nested_res' /\
            f"nested_pairs_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png"
    fig.savefig(figname)
    plt.close('all')
    plt.clf()
    gc.collect()


def dyplot_runplot(results, meta):
    """Plot a summary of the run."""
    rfig, raxes = dyplot.runplot(results)
    plt.savefig(meta.workdir / meta.fitdir / 'nested_res/' /
            f"dyplot_runplot_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png")
    plt.close('all')
    plt.clf()
    gc.collect()


def dyplot_traceplot(results, meta):
    """Plot traces and 1-D marginalized posteriors."""
    tfig, taxes = dyplot.traceplot(results, labels=meta.labels)
    plt.savefig(meta.workdir / meta.fitdir / 'nested_res' /
            f"dyplot_traceplot_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png")
    plt.close('all')
    plt.clf()
    gc.collect()


def dyplot_cornerplot(results, meta):
    # Plot the 2-D marginalized posteriors.
    cfig, caxes = dyplot.cornerplot(results, show_titles=True, title_fmt='.4',labels=meta.labels, color='blue', hist_kwargs=dict(facecolor='blue', edgecolor='blue'))
    plt.savefig(meta.workdir / meta.fitdir / 'nested_res' /
                f"dyplot_cornerplot_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.png")
    plt.close('all')
    plt.clf()
    gc.collect()


def mcmc_rprs(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta):
    """Plots the spectrum (rprs vs wvl) as resulting from the MCMC."""
    vals_mcmc = np.array(vals_mcmc)
    errs_lower_mcmc = np.array(errs_lower_mcmc)
    errs_upper_mcmc = np.array(errs_upper_mcmc)

    rp_idx = np.where(np.array(meta.labels) == 'rp')[0][0]
    plt.errorbar(meta.wavelength_list, vals_mcmc.T[rp_idx],
                 yerr=(errs_lower_mcmc.T[rp_idx], errs_upper_mcmc.T[rp_idx]),
                 fmt='.', c='darkblue', alpha=0.9)

    plt.xlabel('Wavelength (micron)')
    plt.ylabel('Transit Depth (ppm)')
    plt.savefig(meta.workdir / meta.fitdir / 'mcmc_res' / 'mcmc_rprs.png',
                dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()


def nested_rprs(vals_nested, errs_lower_nested, errs_upper_nested, meta):
    """Plots the spectrum (rprs vs wvl) as resulting from the nested sampling."""
    vals_nested = np.array(vals_nested)
    errs_lower_nested = np.array(errs_lower_nested)
    errs_upper_nested = np.array(errs_upper_nested)

    rp_idx = np.where(np.array(meta.labels) == 'rp')[0][0]

    plt.errorbar(meta.wavelength_list, vals_nested.T[rp_idx],
                 yerr=(errs_lower_nested.T[rp_idx], errs_upper_nested.T[rp_idx]),
                 fmt='.', c='darkblue', alpha=0.9)

    plt.xlabel('Wavelength (micron)')
    plt.ylabel('Transit Depth (ppm)')
    plt.savefig(meta.workdir / meta.fitdir / 'nested_res' / 'nested_rprs.png',
                dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()
    gc.collect()
