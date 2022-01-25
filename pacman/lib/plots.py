import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import os
from astropy.time import Time
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from ..lib import util
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator


##00
def mjd_to_utc(time):
    """
    Converts a list of MJDs to a list of dates in years.
    """
    res = []
    for i in time:
        res.append(Time(i, format='mjd').jyear)
    return np.array(res)


def utc_to_mjd(time):
    """
    Converts a list of dates in years to a list of MJDs.
    """
    res = []
    for i in time:
        res.append(Time(i, format='jyear').mjd,)
    return np.array(res)


def obs_times(meta, times, ivisits, iorbits, updated=False):
    """
    Plot of the visit index as a function of observed time for the observations. Includes a table with the number of orbits in each visit and a zoom into each visit.

    Parameters
    ---------------
    updated:bool
        If the user decided to not use all visits but set some "which_visits" in the pcf, this bool is need to save a plot for all files in the data directory and a plot for the onces defined with "which_visits". It prevents that when the function is being called again, the previous plot isnt overwritten.

    """
    fig = plt.figure(figsize=(12,5))

    Nvisits = np.max(ivisits)+1 #Number of visits

    gs = GridSpec(nrows=Nvisits, ncols=2)

    ax1 = fig.add_subplot(gs[:, 0])

    for i in range(Nvisits):
        axn = fig.add_subplot(gs[i, 1])
        if i == 0:
            axn.set_title('each visit')
        axn.text(0.95, 0.41, 'v{0}'.format(i), horizontalalignment='center', verticalalignment='center', transform=axn.transAxes)
        visit_mask = ivisits == i
        times_i = times[visit_mask]
        axn.scatter(times_i, np.ones(len(times_i)), marker='s', c='r', s=5, alpha=0.5)
        axn.get_xaxis().set_ticks([])
        axn.get_yaxis().set_ticks([])
        axn.axes.xaxis.set_visible(False)
        axn.axes.yaxis.set_visible(False)

    ax1.scatter(times, ivisits, marker='s', c='r', s=5, alpha=0.5)
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_ylabel('visit number')
    ax1.set_xlabel('time (MJD)')
    ax1.set_ylim(np.max(ivisits)+0.9, np.min(ivisits)-0.9) #turn around y axis
    ax1sec = ax1.secondary_xaxis('top', functions=(mjd_to_utc, utc_to_mjd))
    ax1sec.set_xlabel('time (year)')

    #It checks the orbit index right before it goes back to 0 (because a new visit stars).
    #We also have to append the very last orbit number
    norbits = iorbits[np.where(np.diff(iorbits) < 0)]
    norbits = np.append(norbits, iorbits[-1:])+1
    table_data = np.array([np.arange(len(norbits), dtype=int), norbits])
    col_names = ['ivisit', '#orbits']
    table = ax1.table(cellText=table_data.T,cellLoc='center',
              colLabels=col_names,
              loc='right',bbox=[-0.45, 0., 0.3, 1])
    table.scale(0.3, 3.0)

    plt.subplots_adjust(wspace=0.04)

    if meta.save_obs_times_plot:
        if not os.path.isdir(meta.workdir + '/figs/s00_obs_dates/'):
            os.makedirs(meta.workdir + '/figs/s00_obs_dates/')
        if updated==False:
            plt.savefig(meta.workdir + '/figs/s00_obs_dates/obs_dates_all.png', bbox_inches='tight', pad_inches=0.05, dpi=120)
        else:
            plt.savefig(meta.workdir + '/figs/s00_obs_dates/obs_dates.png', bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


## 02
def barycorr(x,y,z,time, obsx, obsy, obsz, coordtable, meta):
    """
    This function plots the vectorfile positions of HST and where the observations where taken

    Parameters
    ---------------
    x: array
        X position from vectorfile
    y: array
        Y position from vectorfile
    z: array
        Z position from vectorfile
    time: array
        times from the vectorfile
    obsx: array
        X position of observations
    obsy: array
        Y position of observations
    obsz: array
        Z position of observations
    coordtable
        a list of files containing the vector information of HST downloaded in s01
    meta
        the name of the metadata file

    Returns
    ----------
    Saves and/or Shows a plot

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """
    plt.rcParams["figure.figsize"] = (8, 6)
    fig = plt.figure(1001)
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(x, y, z, c=time, s=10, cmap=cm.inferno)
    ax.scatter(x[0], y[0], z[0], s=200, marker='x', c='r', label='horizons start')
    ax.scatter(x[-1], y[-1], z[-1], s=200, marker='x', c='b', label='horizons end')
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
        if not os.path.isdir(meta.workdir + '/figs/s02_barycorr/'):
            os.makedirs(meta.workdir + '/figs/s02_barycorr/')
        plt.savefig(meta.workdir + '/figs/s02_barycorr/bjdcorr_{0}.png'.format(coordtable.split('/')[-1].split('.')[0]), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


## 03
def smooth(meta, x, y, y_smoothed):
    plt.plot(x, y, label='raw spectrum')
    plt.plot(x, y_smoothed, label='smoothed spectrum')
    plt.legend(loc=1)
    plt.xscale('log')
    plt.xlabel('wavelength (m)')
    plt.ylabel('normalized intensity')
    if meta.save_smooth_plot:
        if not os.path.isdir(meta.workdir + '/figs/s03_smooth/'):
            os.makedirs(meta.workdir + '/figs/s03_smooth/')
        plt.savefig(meta.workdir + '/figs/s03_smooth/smooth.png', bbox_inches='tight',
                    pad_inches=0.05,
                    dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def refspec(bp_wvl, bp_val, sm_wvl, sm_flux, ref_wvl, ref_flux, meta):
    """
    Plots the bandpass, the stellar spectrum and the product of the both
    """
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.figure(1003)
    plt.clf()
    plt.plot(bp_wvl, bp_val, label='bandpass')
    plt.plot(sm_wvl, sm_flux, label='stellar spectrum ({0})'.format(meta.sm))
    plt.plot(ref_wvl, ref_flux, label='stellar spectrum * bandpass')
    plt.xscale('log')
    if meta.grism == 'G102': plt.xlim(0.6e-6, 1.4e-6)
    elif meta.grism == 'G141': plt.xlim(0.8e-6, 2e-6)
    plt.xlabel('wavelength (m)')
    plt.ylabel('norm. intensity')
    plt.legend(loc=4)
    plt.tight_layout()
    if meta.save_refspec_plot:
        if not os.path.isdir(meta.workdir + '/figs/s03_refspec/'):
            os.makedirs(meta.workdir + '/figs/s03_refspec/')
        plt.savefig(meta.workdir + '/figs/s03_refspec/refspec.png', bbox_inches='tight',
                    pad_inches=0.05,
                    dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


## 10
def image_quick(ima, i, meta):
    """
    This plots the full direct image.
    """
    plt.figure(10044)
    plt.clf()
    fig, ax = plt.subplots(1,1, figsize=(8,8))

    plt.rcParams['image.cmap'] = 'viridis'

    nrow = len(ima[1].data[:, 0])
    ncol = len(ima[1].data[0, :])

    plt.suptitle("Direct image #{0}, visit #{1}, orbit #{2}".format(i, meta.ivisit_di[i], meta.iorbit_di[i]), y=0.94)

    ax.title.set_text('Full Direct Image')
    im = ax.imshow(ima[1].data * ima[0].header['exptime'], origin='lower', vmin=0, vmax=500)

    cmin = meta.di_cmin
    cmax = meta.di_cmax
    rmin = meta.di_rmin
    rmax = meta.di_rmax
    ax.plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.8)

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
        if not os.path.isdir(meta.workdir + '/figs/s10_images/'):
            os.makedirs(meta.workdir + '/figs/s10_images/')
        plt.savefig(meta.workdir + '/figs/s10_images/quick_di{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=180)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def image(dat, ima, results, i, meta):
    """
    This plots the full direct image with the guess of the target (defined using di_rmin, etc.) marked as a red box.
    It also plots a zoom into the guess position of the target with the gaussian fit solution marked with a cross.
    """
    plt.figure(1004)
    plt.clf()
    fig, ax = plt.subplots(1,2, figsize=(10,6))

    plt.rcParams['image.cmap'] = 'viridis'

    nrow = len(ima[1].data[:, 0])
    ncol = len(ima[1].data[0, :])

    plt.suptitle("Direct image #{0}, visit #{1}, orbit #{2}".format(i, meta.ivisit_di[i], meta.iorbit_di[i]), y=0.885)

    ax[0].title.set_text('Full Direct Image')
    ax[0].imshow(ima[1].data * ima[0].header['exptime'], origin='lower', vmin=0, vmax=500)
    cmin = meta.di_cmin
    cmax = meta.di_cmax
    rmin = meta.di_rmin
    rmax = meta.di_rmax
    ax[0].plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.8)
    ax[0].set_xlabel('columns')
    ax[0].set_ylabel('rows')

    ax[1].title.set_text('Cutout')
    image_cutout = dat * ima[0].header['exptime']
    vmax = np.max(image_cutout)
    im = ax[1].imshow(image_cutout, origin='lower', vmin=0, vmax=vmax/10)
    ax[1].plot(results[2], results[3], marker='x', color='orange', markeredgewidth=3., ms=10, label='centroid',
             linestyle="none")
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
        if not os.path.isdir(meta.workdir + '/figs/s10_images/'):
            os.makedirs(meta.workdir + '/figs/s10_images/')
        plt.savefig(meta.workdir + '/figs/s10_images/di_{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')



## 20
def sp2d(d, meta, i):
    """
    Plot the spectrum with a low vmax to make the background better visible
    """
    plt.imshow(d[1].data, origin = 'lower',  vmin=0, vmax=500)
    plt.colorbar()
    plt.tight_layout()
    plt.title('Spectrum w/ low vmax, visit {0}, orbit {1}'.format(meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.save_sp2d_plot:
        if not os.path.isdir(meta.workdir + '/figs/sp2d'):
            os.makedirs(meta.workdir + '/figs/sp2d')
        plt.savefig(meta.workdir + '/figs/sp2d/sp2d_{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def trace(d, meta, visnum, orbnum, i):
    """
    Plots the spectrum together with the trace
    #TODO: Q: The plot shouldnt change between observations in the same orbit (same direct image)
    """
    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry141 as geo

    trace = geo.trace(meta.refpix[:,1], meta.refpix[:,2])                #determines trace coefficients
    trace_i = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_i + meta.LTV1)          #start of trace
    trace_f = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_f + meta.LTV1)          #end of trace
    tracex = np.linspace(trace_i, trace_f,100)                    #x values over which to compute trace
    tracey = meta.refpix[orbnum,1] + meta.LTV2 + trace[0][orbnum] + \
    trace[1][orbnum]*(tracex - tracex[0])                    #y values of trace
    plt.imshow(d[1].data, origin = 'lower', vmin=0, vmax=4000)            #plots raw image
    plt.colorbar()
    plt.plot(tracex, tracey, color='yellow', linewidth=2)                #plots trace on raw image frame
    plt.title('Spectrum w/ trace (line), visit {0}, orbit {1}'.format(meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.save_trace_plot:
        if not os.path.isdir(meta.workdir + '/figs/trace/'):
            os.makedirs(meta.workdir + '/figs/trace/')
        plt.savefig(meta.workdir + '/figs/trace/trace_{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close()
    else:
        plt.show()
        plt.close()


def bkg_hist(fullframe_diff, skymedian, meta, i, ii):
    """
    Plot saving a histogram of the fluxes in the up-the-ramp sample. Showing the user decided background threshold and the median flux below the threshold.
    """
    histo = fullframe_diff.flatten()
    fig, ax = plt.subplots(2,1)
    yy, xx, _ = ax[0].hist(histo, int(len(histo)/2000), facecolor='k', alpha=0.2)
    plt.suptitle('UpTheRamp {0}-{1}, visit {2}, orbit {3}'.format(i, ii, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    ax[0].axvline(meta.background_thld, lw=2, ls='-', c='b', label='background_thld')
    ax[0].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Flux')
    ax[0].set_ylabel('log. Frequ.')
    ax[0].legend(loc=1)
    ax[1].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian {0:.3g}'.format(skymedian))
    ax[1].set_xlabel('Flux')
    ax[1].set_ylabel('lin. Frequ.')
    zoom = 150
    #yyy, xxx, _ = ax[1].hist(histo, 51, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    ax[1].hist(histo, 41, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    #xmax = xxx[np.argsort(yyy)[::-1]][0]
    ax[1].legend(loc=1)
    plt.tight_layout()
    #data_new = histo[((skymedian - zoom) < histo) & (histo < (skymedian + zoom))]
    #var = util.median_abs_dev(data_new)
    #from scipy.stats import norm
    # Fit a normal distribution to the data:
    #mu, std = norm.fit(data_new)
    #x_new = np.linspace(skymedian - zoom, skymedian + zoom, 100)
    #p = norm.pdf(x_new, mu, std)
    #ax[1].plot(x_new, p*max(yyy)*10, 'k', linewidth=2)
    #ax[1].axvline(x_new[np.argmax(p)], lw=2, ls='--', c='g', label='gaussian {0:.3g}'.format(x_new[np.argmax(p)]))

    if meta.save_bkg_hist_plot:
        if not os.path.isdir(meta.workdir + '/figs/bkg_hist/'):
            os.makedirs(meta.workdir + '/figs/bkg_hist/')
        plt.savefig(meta.workdir + '/figs/bkg_hist/bkg_hist{0}-{1}.png'.format(i,ii), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close()


def utr(diff, meta, i, ii, orbnum, rowsum, rowsum_absder, peaks):
    """
    Saves a plot of up-the-ramp sample, the row by row sum and the derivate of the latter. It furthermore shows the aperture used for the analysis.
    """
    cmin = int(meta.refpix[
                   orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_i + meta.LTV1 + meta.offset  # determines left column for extraction (beginning of the trace)
    cmax = min(int(meta.refpix[orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_f + meta.LTV1 - meta.offset,
               meta.subarray_size)  # right column (end of trace, or edge of detector)

    fig, ax = plt.subplots(1,3, figsize=(10,6))

    im = ax[0].imshow(diff, vmin=0, vmax=300, origin='lower')
    #ax[2].axhline(idx, c = 'b', ls='--')
    #ax[2].plot([cmin-cmin, cmin-cmin, cmax-cmin, cmax-cmin, cmin-cmin], [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.85)
    ax[0].set_xlim(cmin-cmin, cmax-cmin)
    ax[0].set_ylim(meta.rmin, meta.rmax)
    ax[0].axhline(min(peaks), c='r', ls='--', lw=2)
    ax[0].axhline(max(peaks), c='r', ls='--', lw=2)
    ax[0].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3)
    ax[0].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3)
    ax[0].title.set_text('spectrum, cmin to cmax')

    p1 = np.arange(len(diff))
    ax[1].axhline(min(peaks), c='r', ls='--', lw=2 )
    ax[1].axhline(max(peaks), c='r', ls='--', lw=2 )
    ax[1].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3 )
    ax[1].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3 )
    ax[1].plot(rowsum, p1)
    ax[1].set_ylim(meta.rmin, meta.rmax)
    ax[1].set_xlabel('lin Flux')
    ax[1].title.set_text('sum row Flux')

    p2 = (p1[1:] + p1[:-1]) / 2
    ax[2].axhline(min(peaks)-meta.window, c='r', ls='-', lw=3 )
    ax[2].axhline(max(peaks)+meta.window, c='r', ls='-', lw=3 )
    ax[2].scatter(rowsum_absder[peaks], peaks, marker = 'x', c='r')
    ax[2].plot(rowsum_absder, p2)
    ax[2].set_ylim(meta.rmin, meta.rmax)
    ax[2].set_xlabel('lin Flux')
    ax[2].title.set_text('Derivative')

    plt.colorbar(im)
    plt.tight_layout()
    plt.suptitle('UpTheRamp {0}-{1}, visit {2}, orbit {3}'.format(i, ii, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.save_utr_plot:
        if not os.path.isdir(meta.workdir + '/figs/utr/'):
            os.makedirs(meta.workdir + '/figs/utr/')
        plt.savefig(meta.workdir + '/figs/utr/utr{0}-{1}.png'.format(i, ii), bbox_inches='tight', pad_inches=0.05, dpi=100)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def sp1d(template_waves, spec_box, meta, i, spec_opt=False):
    """
    Plots the resulting spectrum. If the user did optimal extraction, a comparison between optextr and box sum will be shown.
    """
    plt.clf()
    #fig, ax = plt.subplots(1, 1, figsize=(6.4, 6.4))
    plt.plot(template_waves, spec_box, c='r', label='box sum')
    if meta.opt_extract:
        plt.plot(template_waves, spec_opt, c='k', label='optimal extraction')
    plt.legend()
    plt.title('sp1d, visit {0}, orbit {1}'.format(meta.ivisit_sp[i], meta.iorbit_sp[i]))
    plt.tight_layout()
    if meta.save_sp1d_plot:
        if not os.path.isdir(meta.workdir + '/figs/sp1d/'):
            os.makedirs(meta.workdir + '/figs/sp1d/')
        plt.savefig(meta.workdir + '/figs/sp1d/sp1d_{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def bkg_evo(bkg_evo, meta):
    """
    Plot of the background flux as a function of up the ramp sample
    """
    plt.clf()
    fig, ax = plt.subplots(1,1)
    ax.plot(range(len(bkg_evo)), bkg_evo)
    ax.set_xlabel('# diff image')
    ax.set_ylabel('diff flux')
    #plt.legend()
    plt.tight_layout()
    if meta.save_bkg_evo_plot:
        if not os.path.isdir(meta.workdir + '/figs/bkg_evo/'):
            os.makedirs(meta.workdir + '/figs/bkg_evo/')
        plt.savefig(meta.workdir + '/figs/bkg_evo/bkg_evo.png', bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def sp1d_diff(sp1d_all_diff, meta, wvl_hires):
    """
    Difference of 1D spectrum between two consecutive exposures
    """
    if not os.path.isdir(meta.workdir + '/figs/sp1d_diff/'):
        os.makedirs(meta.workdir + '/figs/sp1d_diff/')
    ylimmin = np.min(sp1d_all_diff)
    ylimmax = np.max(sp1d_all_diff)
    for iiii in range(len(sp1d_all_diff)):
        plt.clf()
        fig, ax = plt.subplots(1, 1, figsize=(6.4, 6.4*1.5/2))
        ax.plot(wvl_hires, sp1d_all_diff[iiii])
        if meta.grism == 'G141':
            ax.set_xlim(10000, 17800)
        elif meta.grism == 'G102':
            ax.set_xlim(5000, 16000)
        ax.set_ylim(ylimmin*1.2, ylimmax*1.2)
        fig.suptitle('Background Diff, from v{0} o{1} to v{2} o{3}'.format(meta.ivisit_sp[iiii], meta.iorbit_sp[iiii],meta.ivisit_sp[iiii+1], meta.iorbit_sp[iiii+1]),fontsize=12, y=0.99)
        plt.tight_layout()
        if meta.save_sp1d_diff_plot:
            plt.savefig(meta.workdir + '/figs/sp1d_diff/sp1d_diff_{0}.png'.format(iiii), bbox_inches='tight', pad_inches=0.05, dpi=120)
            plt.close('all')
        else:
            plt.show()
            plt.close('all')


def utr_aper_evo(peaks_all, meta):
    """
    Plot of the evolution in aperture size
    """
    peaks = np.array(peaks_all)
    peaks = np.sort(peaks)
    #nsamps_mask = []
    #for i in nsamps:
    #    nsamps_mask.append(np.arange(i, dtype=int))
    #nsamps_mask = np.array(nsamps_mask)
    #nsamps_mask[::2] = nsamps_mask[::-1][::2]
    #nsamps_mask = np.array(nsamps_mask).flatten()

    plt.clf()
    fig, ax = plt.subplots(1,1, figsize=(8,4), sharex=True)

    #for ii in range(max(nsamps)):
    #    ax[0].plot(range(len(peaks[:,0][nsamps_mask==ii])), peaks[:,0][nsamps_mask==ii])
    #    ax[0].plot(range(len(peaks[:,1][nsamps_mask==ii])), peaks[:,1][nsamps_mask==ii])
    #    ax[1].plot(range(len(peaks[:,0][nsamps_mask==ii])), np.diff(peaks[nsamps_mask==ii]))

    #ax[0].plot(range(len(peaks[:,0])), peaks[:,0])
    #ax[0].plot(range(len(peaks[:,1])), peaks[:,1])
    ax.plot(range(len(peaks[:,0])), np.diff(peaks))

    ax.set_xlabel('# diff image')
    #ax[0].set_ylabel('peak position')
    ax.set_ylabel('peak distance')

    #plt.legend()
    plt.tight_layout()
    if meta.save_utr_aper_evo_plot:
        if not os.path.isdir(meta.workdir + '/figs/utr_aper_evo/'):
            os.makedirs(meta.workdir + '/figs/utr_aper_evo/')
        plt.savefig(meta.workdir + '/figs/utr_aper_evo/utr_aper_evo.png', bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def refspec_comp(modelx, modely, p0, datax, datay, leastsq_res, meta, i):
    fig, ax = plt.subplots(1,1, figsize=(9,6))
    plt.suptitle('refspec_comp {0}, visit {1}, orbit {2}'.format(i, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    ax.plot(modelx, modely, label='refspec')
    #ax.plot((p0[0] + p0[1] * datax), datay * p0[2], label='initial guess')
    ax.plot((leastsq_res[0] + datax * leastsq_res[1]), datay * leastsq_res[2],
             label='spectrum fit, wvl = ({0:.5g}+{1:.5g}*x), {2:.5g}'.format(leastsq_res[0], leastsq_res[1], leastsq_res[2]))
    if meta.grism == 'G141':
        ax.set_xlim(9000, 20000)
    elif meta.grism == 'G102':
        ax.set_xlim(5000, 16000)
    ax.set_xscale('log')
    ax.set_xlabel('Wavelength (angstrom)')
    ax.set_xlabel('rel. Flux')
    plt.legend()
    plt.tight_layout()
    if meta.save_refspec_comp_plot:
        if not os.path.isdir(meta.workdir + '/figs/drift/'):
            os.makedirs(meta.workdir + '/figs/drift/')
        plt.savefig(meta.workdir + '/figs/drift/drift{0}.png'.format(i), bbox_inches='tight', pad_inches=0.05, dpi=120)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')



##21
def plot_wvl_bins(w_hires, f_interp, wave_bins, wvl_bins, dirname):
    """
    Plot of a 1D spectrum and the bins
    """
    plt.plot(w_hires, f_interp)
    for wave in wave_bins: plt.axvline(wave, color='0.5')
    plt.ylabel("Photelectrons")
    plt.xlabel("Wavelength (angstroms)")
    plt.tight_layout()
    plt.savefig(dirname + '/spec_bins{0}.png'.format(wvl_bins))
    #plt.show()










#30

from matplotlib import rc
import matplotlib
import seaborn as sns
from .formatter import FormatParams
from .model import Model, calc_sys, calc_astro

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction": "in", "ytick.direction": "in"})
matplotlib.rcParams.update({'lines.markeredgewidth': 0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset': False})


def plot_raw(data, meta):
    #palette = sns.color_palette("husl", data.nvisit)
    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)
    fig, ax = plt.subplots(data.nvisit, 1, figsize=(6.4, 2.4*data.nvisit), sharex=True)
    if data.nvisit>1:
        for i in range(data.nvisit):
            ind = data.vis_num == i
            #plt.subplot((data.nvisit) * 100 + 10 + i + 1)
            ax[i].plot(data.t_vis[ind] * 24., data.flux[ind], marker='o', \
                     markersize=3.0, linestyle="none", \
                     label="Visit {0}".format(i), color=palette[i])
            ax[i].set_xlim(((data.t_vis.min() - 0.02) * 24., (data.t_vis.max() + 0.05) * 24.))
            ax[i].set_ylim((0.998 * data.flux.min(), 1.002 * data.flux.max()))
            ax[i].legend(loc=1)
            ax[i].set_ylabel("Flux (e-)")
    else:

        #plt.subplot((data.nvisit) * 100 + 10 + i + 1)
        ax.plot(data.t_vis * 24., data.flux, marker='o', \
                 markersize=3.0, linestyle="none", \
                 label="Visit 0", color=palette[0])
        ax.set_xlim(((data.t_vis.min() - 0.02) * 24., (data.t_vis.max() + 0.05) * 24.))
        ax.set_ylim((0.998 * data.flux.min(), 1.002 * data.flux.max()))
        ax.legend(loc=1)
        ax.set_ylabel("Flux (e-)")
    plt.xlabel("Time after visit start (hours)")
    fig.suptitle('Filename: {0}'.format(meta.run_file.split('/')[-1]), fontsize=15, y=0.998)
    plt.tight_layout()
    plt.savefig(meta.workdir + meta.fitdir + "/raw_lc_{0}_{1}.png".format(meta.run_file.split('/')[-1], meta.fittime))
    plt.close()


# COMPUTE ROOT-MEAN-SQUARE AND STANDARD ERROR OF DATA FOR VARIOUS BIN SIZES
def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
    # data    = fit.normresiduals
    # maxnbin = maximum # of bins
    # binstep = Bin step size

    # bin data into multiple bin sizes
    npts = data.size
    if maxnbins is None:
        maxnbins = npts / 10.
    binsz = np.arange(1, maxnbins + binstep, step=binstep, dtype=int)
    nbins = np.zeros(binsz.size, dtype=int)
    rms = np.zeros(binsz.size)
    rmserr = np.zeros(binsz.size)
    for i in range(binsz.size):
        nbins[i] = int(np.floor(data.size / binsz[i]))
        bindata = np.zeros(nbins[i], dtype=float)
        # bin data
        # ADDED INTEGER CONVERSION, mh 01/21/12
        for j in range(nbins[i]):
            bindata[j] = data[j * binsz[i]:(j + 1) * binsz[i]].mean()
        # get rms
        rms[i] = np.sqrt(np.mean(bindata ** 2))
        rmserr[i] = rms[i] / np.sqrt(2. * nbins[i])
    # expected for white noise (WINN 2008, PONT 2006)
    stderr = (data.std() / np.sqrt(binsz)) * np.sqrt(nbins / (nbins - 1.))
    if isrmserr is True:
        return rms, stderr, binsz, rmserr
    else:
        return rms, stderr, binsz




# Plot RMS vs. bin size looking for time-correlated noise
def rmsplot(model, data, meta):

    time = data.time
    residuals = model.norm_resid
    residuals = residuals[np.argsort(time)]

    rms, stderr, binsz = computeRMS(residuals, binstep=1)
    normfactor = 1e-6
    plt.rcParams.update({'legend.fontsize':11})
    plt.figure(1111, figsize=(8,6))
    plt.clf()
    plt.suptitle(' Correlated Noise', size=16)
    plt.loglog(binsz, rms/normfactor, color='black', lw=1.5, label='Fit RMS', zorder=3)    # our noise
    plt.loglog(binsz, stderr/normfactor, color='red', ls='-', lw=2, label='Std. Err.', zorder=1) # expected noise
    plt.xlim(0.95, binsz[-1]*2)
    plt.ylim(stderr[-1]/normfactor/2., stderr[0]/normfactor*2.)
    plt.xlabel("Bin Size", fontsize=14)
    plt.ylabel("RMS (ppm)", fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.legend()
    #if savefile != None:
    plt.savefig(meta.workdir + meta.fitdir + '/corr_lc_{0}_len{1}.png'.format(meta.fittime, len(time)))
    plt.close()

from astropy.stats import sigma_clip
import time


def plot_fit_lc(data, fit, meta, mcmc=False):
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    plt.clf()
    fig, ax = plt.subplots(2,1)
    p = FormatParams(fit.params, data)  # FIXME

    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)

    # ind = model.phase > 0.5
    # model.phase[ind] -= 1.

    # calculate a range of times at higher resolution to make model look nice
    phase_hr = np.linspace(fit.phase.min() - 0.05, fit.phase.max() + 0.05, 1000)
    t_hr = phase_hr * p.per[0] + p.t0[0] + data.toffset

    # colors = ['blue', 'red', 'orange', 'gray']

    # plot data
    # plot best fit model from first visit

    ax[0].plot(phase_hr, calc_astro(t_hr, fit.params, data, fit.myfuncs, 0))


    # plot systematics removed data
    for i in range(data.nvisit):
        ind = data.vis_num == i
        ax[0].plot(fit.phase[ind], fit.data_nosys[ind], color=palette[i], marker='o', markersize=3, linestyle="none")

    clip_mask = np.ma.getmask(sigma_clip(fit.resid, sigma=meta.run_clipsigma, maxiters=1))
    ax[0].scatter(fit.phase[clip_mask], fit.data_nosys[clip_mask], marker='x', s=100, c='r')

    ax[1].axhline(np.median(1.0e6 * fit.norm_resid), c='r')
    clip_line = np.std(1.0e6 * fit.norm_resid)
    ax[1].axhline(np.median(1.0e6 * fit.norm_resid)+clip_line*meta.run_clipsigma, c='r', ls='--')
    ax[1].axhline(np.median(1.0e6 * fit.norm_resid)-clip_line*meta.run_clipsigma, c='r', ls='--')

    # add labels/set axes
    # xlo, xhi = np.min(model.phase)*0.9, np.max(model.phase)*1.1
    xlo, xhi = -0.1, 0.1
    ax[0].set_xlim(xlo, xhi)
    ax[0].set_ylabel("Relative Flux")

    # annotate plot with fit diagnostics
    #ax = plt.gca()
    ax[0].text(0.85, 0.29,
            '$\chi^2_\\nu$:    ' + '{0:0.2f}'.format(fit.chi2red) + '\n'
            + 'obs. rms:  ' + '{0:0d}'.format(int(fit.rms)) + '\n'
            + 'exp. rms:  ' + '{0:0d}'.format(int(fit.rms_predicted)),
            verticalalignment='top', horizontalalignment='left',
            transform=ax[0].transAxes, fontsize=12
            )

    # plot fit residuals
    ax[1].axhline(0, zorder=1, color='0.2', linestyle='dashed')

    for i in range(data.nvisit):
        ind = data.vis_num == i
        ax[1].plot(fit.phase[ind], 1.0e6 * fit.norm_resid[ind], color=palette[i], marker='o', markersize=3,
                 linestyle="none")

    # add labels/set axes
    ax[1].set_xlim(xlo, xhi)
    ax[1].set_ylabel("Residuals (ppm)")
    ax[1].set_xlabel("Orbital phase")

    plt.tight_layout()
    if mcmc:
        plt.savefig(meta.workdir + meta.fitdir + "/mcmc_lc_{0}_len{1}_{2}.png".format(meta.fittime, len(fit.phase), datetime))
    else:
        plt.savefig(meta.workdir + meta.fitdir + "/white_lc_{0}_len{1}_{2}.png".format(meta.fittime, len(fit.phase), datetime))
    # plt.waitforbuttonpress(0) # this will wait for indefinite time
    plt.close()


def plot_fit_lc2(data, fit, meta, mcmc=False):
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    plt.clf()
    fig, ax = plt.subplots(2,1)
    print(fit.params)
    p = FormatParams(fit.params, data)  # FIXME
    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)

    # ind = model.phase > 0.5
    # model.phase[ind] -= 1.
    # calculate a range of times at higher resolution to make model look nice
    phase_hr = np.linspace(fit.phase.min() - 0.05, fit.phase.max() + 0.05, 1000)
    #print(phase_hr)
    print(p.per)
    print(p.per[0])
    print(p.t0)
    print(p.t0[0])
    print(data.toffset)
    t_hr = phase_hr * p.per[0] + p.t0[0] + data.toffset

    # colors = ['blue', 'red', 'orange', 'gray']

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
    #ax = plt.gca()
    ax[0].text(0.85, 0.29,
            '$\chi^2_\\nu$:    ' + '{0:0.2f}'.format(fit.chi2red) + '\n'
            + 'obs. rms:  ' + '{0:0d}'.format(int(fit.rms)) + '\n'
            + 'exp. rms:  ' + '{0:0d}'.format(int(fit.rms_predicted)),
            verticalalignment='top', horizontalalignment='left',
            transform=ax[0].transAxes, fontsize=12
            )

    # plot fit residuals
    ax[1].axhline(0, zorder=1, color='0.2', linestyle='dashed')

    for i in range(data.nvisit):
        ind = data.vis_num == i
        ax[1].plot(fit.phase[ind], 1.0e6 * fit.norm_resid[ind], color=palette[i], marker='o', markersize=3,
                 linestyle="none")

    # add labels/set axes
    ax[1].set_xlim(xlo, xhi)
    ax[1].set_ylabel("Residuals (ppm)")
    ax[1].set_xlabel("Orbital phase")

    plt.tight_layout()
    if mcmc:
        plt.savefig(meta.workdir + meta.fitdir + "/mcmc_lc_{0}_len{1}_{2}.png".format(meta.fittime, len(fit.phase), datetime))
    else:
        plt.savefig(meta.workdir + meta.fitdir + "/white_lc_{0}_len{1}_{2}.png".format(meta.fittime, len(fit.phase), datetime))
    # plt.waitforbuttonpress(0) # this will wait for indefinite time
    plt.close()


def plot_chains(ndim, sampler, nburn, labels, meta):
    plt.clf()
    print(ndim)
    fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, ndim))
    for i in range(0, ndim):
        axes[i].plot(sampler.chain[:, nburn:, i].T, alpha=0.4)
        # axes[i].yaxis.set_major_locator(MaxNLocator(5))
        axes[i].set_ylabel(labels[i])
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    if nburn == 0:
        fig.savefig(meta.workdir + meta.fitdir + "/mcmc_chains_{0}_{1}.png".format(meta.run_file.split('/')[-1], meta.fittime))
    else:
        fig.savefig(meta.workdir + meta.fitdir + "/mcmc_chains_noburn_{0}_{1}.png".format(meta.run_file.split('/')[-1], meta.fittime))
    plt.close()