import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from ..lib import util


def barycorr(x,y,z,time, obsx, obsy, obsz, coordtable, meta):
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
        plt.savefig(meta.workdir + '/ancil/horizons/bjdcorr_{0}.png'.format(coordtable.split('/')[-1].split('.')[0]))
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def bandpass(wvl, bp_val, grism, i, meta):
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.figure(1002)
    plt.clf()
    plt.plot(wvl, bp_val, c='C0')
    plt.xlabel('Angstrom')
    plt.ylabel('Throughput')
    plt.title('{0}_v{1}'.format(grism, i))
    plt.tight_layout()
    if meta.save_bandpass_plot:
        plt.savefig(meta.workdir + '/ancil/bandpass/bandpass_v{0}.png'.format(i))
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def refspec(x, y, wvl_g, flux_g, f, vi, meta):
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.figure(1003)
    plt.clf()
    plt.plot(x, y, label='stellar spectrum')
    plt.plot(wvl_g, flux_g, label='bandpass')
    plt.plot(wvl_g, f(wvl_g)*flux_g, label='stellar spectrum * bandpass')
    plt.xlim(0.7, 2)
    plt.legend()
    plt.tight_layout()
    if meta.save_refspec_plot:
        plt.savefig(meta.workdir + '/ancil/bandpass/refspec_v{0}.png'.format(vi))
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def image(dat, ima, results, i, meta):
    plt.figure(1004)
    plt.clf()
    fig, ax = plt.subplots(1,2, figsize=(9,5.3))

    plt.rcParams['image.cmap'] = 'viridis'

    nrow = len(ima[1].data[:, 0])
    ncol = len(ima[1].data[0, :])

    plt.suptitle("Direct image #{0}, visit #{1}, orbit #{2}".format(i, meta.ivisit_di[i], meta.iorbit_di[i]))

    ax[0].title.set_text('Full Direct Image')
    ax[0].imshow(ima[1].data * ima[0].header['exptime'], origin='lower', vmin=0, vmax=500)
    cmin = meta.di_cmin
    cmax = meta.di_cmax
    rmin = meta.di_rmin
    rmax = meta.di_rmax
    ax[0].plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], lw=1, c='r', alpha=0.85)
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
        if not os.path.isdir(meta.workdir + '/figs/images/'):
            os.makedirs(meta.workdir + '/figs/images/')
        plt.savefig(meta.workdir + '/figs/images/{0}.png'.format(i))
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def bkg(d, meta, i):
    plt.imshow(d[1].data, origin = 'lower',  vmin=0, vmax=300)
    plt.colorbar()
    plt.tight_layout()
    plt.title('Background, visit {0}, orbit {1}'.format(meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.all_image_output:
        if not os.path.isdir(meta.workdir + '/figs/background'):
            os.makedirs(meta.workdir + '/figs/background')
        plt.savefig(meta.workdir + '/figs/background/background{0}.png'.format(i))
        plt.close('all')
    if meta.all_image_show:
        plt.show()
        plt.close('all')


def plot_trace(d, meta, visnum, orbnum, i):
    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry as geo
    else:
        print('Error: GRISM in obs_par.cf is neither G102 nor G141!')

    trace = geo.trace(meta.refpix[:,1], meta.refpix[:,2])                #determines trace coefficients

    trace_i = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_i + meta.LTV1)          #start of trace
    trace_f = meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale + float(meta.BEAMA_f + meta.LTV1)          #end of trace
    tracex = np.linspace(trace_i, trace_f,100)                    #x values over which to compute trace
    tracey = meta.refpix[orbnum,1] + meta.LTV2 + trace[0][orbnum] + \
    trace[1][orbnum]*(tracex - tracex[0])                    #y values of trace
    plt.imshow(d[1].data, origin = 'lower', vmin=0, vmax=40000)            #plots raw image
    plt.colorbar()
    plt.plot(tracex, tracey, color='yellow', linewidth=2)                #plots trace on raw image frame
    plt.title('Trace marked with a line, visit {0}, orbit {1}'.format(meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.save_trace_plot:
        if not os.path.isdir(meta.workdir + '/figs/trace/'):
            os.makedirs(meta.workdir + '/figs/trace/')
        plt.savefig(meta.workdir + '/figs/trace/{0}.png'.format(i))
        plt.close()
    else:
        plt.show()
        plt.close()
    return [trace_i, trace_f]


def bkg_histogram(fullframe_diff, skymedian, meta, i, ii):

    histo = fullframe_diff.flatten()
    fig, ax = plt.subplots(2,1)
    yy, xx, _ = ax[0].hist(histo, int(len(histo)/500), facecolor='k', alpha=0.2)
    plt.suptitle('UpTheRamp {0}-{1}, visit {2}, orbit {3}'.format(i, ii, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    ax[0].axvline(meta.background_thld, lw=2, ls='-', c='b', label='background_thld')
    ax[0].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Flux')
    ax[0].set_ylabel('log. Frequ.')
    ax[0].legend()
    ax[1].axvline(skymedian, lw=2, ls='--', c='b', label='skymedian {0:.3g}'.format(skymedian))
    ax[1].set_xlabel('Flux')
    ax[1].set_ylabel('lin. Frequ.')
    zoom = 200
    #yyy, xxx, _ = ax[1].hist(histo, 51, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    ax[1].hist(histo, 51, range=(skymedian - zoom, skymedian + zoom), facecolor='k', alpha=0.2)
    #xmax = xxx[np.argsort(yyy)[::-1]][0]
    ax[1].legend()
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


    if not os.path.isdir(meta.workdir + '/figs/bkg_histogram/'):
        os.makedirs(meta.workdir + '/figs/bkg_histogram/')
    plt.savefig(meta.workdir + '/figs/bkg_histogram/bkg_histogram{0}-{1}.png'.format(i,ii))
    #plt.show()
    plt.close('all')



def uptheramp(diff, meta, i, ii, idx, orbnum):
    rmin, rmax = max(idx - meta.window, 0), min(idx + meta.window, meta.rmax)
    cmin = int(meta.refpix[
                   orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_i + meta.LTV1 + meta.offset  # determines left column for extraction (beginning of the trace)
    cmax = min(int(meta.refpix[orbnum, 2] + meta.POSTARG1 / meta.platescale) + meta.BEAMA_f + meta.LTV1 - meta.offset,
               meta.subarray_size)  # right column (end of trace, or edge of detector)

    fig, ax = plt.subplots(1,2)


    im = ax[1].imshow(diff, vmin=0, vmax=300, origin='lower')
    ax[1].axhline(idx, c = 'b', ls='--')
    ax[1].plot([cmin-cmin, cmin-cmin, cmax-cmin, cmax-cmin, cmin-cmin], [rmin, rmax, rmax, rmin, rmin], lw=2, c='r', alpha=0.85)
    ax[1].set_xlim(cmin-cmin, cmax-cmin)
    ax[1].set_ylim(0, 512)
    ax[0].plot(np.sum(diff, axis = 1), range(len(diff)))
    ax[0].set_xscale('log')
    ax[0].set_xlim(1e4, 1e7)
    ax[0].set_ylim(0,512)
    ax[0].axhline(idx, c = 'b', ls='--')
    ax[0].set_xlabel('log Flux')

    plt.colorbar(im)
    plt.tight_layout()
    plt.suptitle('UpTheRamp {0}-{1}, visit {2}, orbit {3}'.format(i, ii, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    if meta.all_image_output:
        if not os.path.isdir(meta.workdir + '/figs/uptheramp/'):
            os.makedirs(meta.workdir + '/figs/uptheramp/')
        plt.savefig(meta.workdir + '/figs/uptheramp/uptheramp{0}-{1}.png'.format(i, ii))
        plt.close('all')
    if meta.all_image_show:
        plt.show()
        plt.close('all')


def refspec_comp(x_vals, y_vals, modelx, modely, p0, datax, datay, leastsq_res, meta, i):

    fig, ax = plt.subplots(1,1, figsize=(9,6))
    plt.suptitle('refspec_comp {0}, visit {1}, orbit {2}'.format(i, meta.ivisit_sp[i], meta.iorbit_sp[i]))
    ax.plot(x_vals, y_vals, label='refspec')
    ax.plot(modelx, modely, label='refspec smoothed')
    ax.plot((p0[0] + p0[1] * datax), datay * p0[2], label='initial guess')
    ax.plot((leastsq_res[0] + datax * leastsq_res[1]), datay * leastsq_res[2],
             label='spectrum fit, wvl = ({0:.5g}+{1:.5g}*x), {2:.5g}'.format(leastsq_res[0], leastsq_res[1], leastsq_res[2]))
    ax.set_xlim(min(x_vals)/1.5, max(x_vals)*1.5)
    ax.set_xscale('log')
    ax.set_xlabel('Wavelength (angstrom)')
    ax.set_xlabel('rel. Flux')
    plt.legend()
    plt.tight_layout()
    if not os.path.isdir(meta.workdir + '/figs/drift/'):
        os.makedirs(meta.workdir + '/figs/drift/')
    plt.savefig(meta.workdir + '/figs/drift/drift{0}.png'.format(i))
    plt.close('all')

