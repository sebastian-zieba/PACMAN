import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import math
from importlib import reload
import multiprocessing as mp
import os
from tqdm import tqdm
import numpy.ma as ma
from scipy.interpolate import interp1d
from ..lib import plots
from scipy.optimize import leastsq
from .sort_nicely import sort_nicely as sn


#s00

def readfiles(meta):
    """
    Reads in the files saved in datadir and saves them into a list

    Parameters
    -----------
    meta
        metadata object

    Returns
    ----------
    meta
        metadata object but adds segment_list to metadata containing the sorted data fits files

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021

    """
    meta.segment_list = []
    for fname in os.listdir(str(meta.datadir)):
        if fname.endswith(meta.suffix + '.fits'):
            meta.segment_list.append(str(meta.datadir) +'/'+ fname)
    return meta



#s02

def ancil(meta, s10=False, s20=False):
    """
    This function loads in a lot of useful arrays and values into meta

    The following additional information are being loading into meta:

    - **norbit:**   number of orbits
    - **nvisit:**   number of visits
    - **files_sp:** all spectra files
    - **files_di:** all direct image files
    - **ra:** RA of the target in radians (from the header) (Note: this data is taken from the first spectrum file)
    - **dec:** DEC of the target in radians (from the header) (Note: this data is taken from the first spectrum file)
    - **coordtable:** a list of files containing the vector information of HST downloaded in s01


    Parameters
    -----------
    meta
        metadata object
    s10: bool
        Is set to True when s10 is being performed
    s20: bool
        Is set to True when s20 is being performed

    Returns
    ----------
    meta
        metadata object

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """

    filelist = ascii.read(meta.workdir + '/filelist.txt')

    aa = filelist['iorbit'].data
    bb = np.diff(aa)
    meta.norbit, meta.nvisit = len(np.insert(aa[1:][np.where(bb != 0)], 0, aa[0])), len(set(filelist['ivisit'].data))

    meta.mask_sp = np.array([i[0] for i in filelist['instr']]) == 'G'
    meta.mask_di = np.array([i[0] for i in filelist['instr']]) == 'F'
    meta.files_sp = [meta.datadir + '/' + i for i in filelist['filenames'][meta.mask_sp].data]
    meta.files_di = [meta.datadir + '/' + i for i in filelist['filenames'][meta.mask_di].data]

    f = fits.open(meta.files_sp[0])

    meta.ra = f[0].header['ra_targ'] * math.pi / 180.0  # stores right ascension
    meta.dec = f[0].header['dec_targ'] * math.pi / 180.0  # stores declination

    meta.coordtable = []  # table of spacecraft coordinates
    for i in range(max(filelist['ivisit']) + 1): meta.coordtable.append(
        meta.workdir + '/ancil/horizons/' + "/horizons_results_v" + str(i) + ".txt")

    ###
    # 03
    # TODO: Q: Is it okay if I assume everything in datadir uses same Filter Grism pair?
    meta.filter = filelist['instr'][meta.mask_di][0]
    meta.grism = filelist['instr'][meta.mask_sp][0]

    meta.scans_sp = filelist['scan'][meta.mask_sp].data
    meta.iorbit_sp = filelist['iorbit'][meta.mask_sp].data
    meta.ivisit_sp = filelist['ivisit'][meta.mask_sp].data

    # the following list stores the cumulative orbit number
    meta.iorbit_sp_cumulative = np.zeros(len(meta.iorbit_sp), dtype=int)
    c = 0
    for i in range(len(meta.iorbit_sp) - 1):
        if meta.iorbit_sp[i + 1] != meta.iorbit_sp[i]:
            c = c + 1
        meta.iorbit_sp_cumulative[i + 1] = c

    # the following lists store the indices of spectral files where the orbit/visit increases
    meta.new_orbit_idx_sp = np.concatenate(([0], np.where(np.diff(meta.iorbit_sp)!=0)[0]+1))
    meta.new_visit_idx_sp = np.concatenate(([0], np.where(np.diff(meta.ivisit_sp)!=0)[0]+1))

    meta.iorbit_di = filelist['iorbit'][meta.mask_di].data
    meta.ivisit_di = filelist['ivisit'][meta.mask_di].data

    # the following list stores the cumulative orbit number for the direct images
    meta.iorbit_di_cumulative = np.zeros(len(meta.iorbit_di), dtype=int)
    c = 0
    for i in range(len(meta.iorbit_di) - 1):
        if meta.iorbit_di[i + 1] != meta.iorbit_di[i]:
            c = c + 1
        meta.iorbit_di_cumulative[i + 1] = c

    meta.t_mjd_sp = filelist['t_mjd'][meta.mask_sp].data

    meta.t_orbit_sp = filelist['t_orbit'][meta.mask_sp].data
    meta.t_visit_sp = filelist['t_visit'][meta.mask_sp].data

    meta.platescale = 0.13  # IR detector has plate scale of 0.13 arcsec/pixel

    meta.POSTARG1 = f[0].header['POSTARG1']  # x-coordinate of the observer requested target offset
    meta.POSTARG2 = f[0].header['POSTARG2']  # y-coordinate of the observer requested target offset
    meta.LTV1 = int(f[1].header['LTV1'])     # X offset to get into physical pixels
    meta.LTV2 = int(f[1].header['LTV2'])     # Y offset to get into physical pixels

    meta.subarray_size = f[1].header['SIZAXIS1']  # size of subarray

    if meta.grism == 'G102':
        meta.BEAMA_i = 41
        meta.BEAMA_f = 248
    elif meta.grism == 'G141':
        meta.BEAMA_i = 15
        meta.BEAMA_f = 196


    if s10:
        if 't_bjd' in filelist.keys():
            meta.t_bjd_sp = filelist['t_bjd'][meta.mask_sp].data
        if 't_bjd' in filelist.keys():
            meta.t_bjd_di = filelist['t_bjd'][meta.mask_di].data

    if s20:
        meta.rdnoise = 22.0 #read noise
        if meta.grism == 'G102':
            meta.flat = meta.pacmandir + '/ancil/flats/WFC3.IR.G102.flat.2.fits'
        elif meta.grism == 'G141':
            meta.flat = meta.pacmandir + '/ancil/flats/WFC3.IR.G141.flat.2.fits'

    return meta

#03
def gaussian_kernel(meta, x, y):
    """
    https://matthew-brett.github.io/teaching/smoothing_intro.html
    """

    sigma = meta.smooth_sigma * 1e-10

    y = y / max(y)
    y_smoothed = np.zeros(y.shape)

    for x_idx, x_val in enumerate(x):
        kernel = np.exp(-(x - x_val) ** 2 / (2 * sigma ** 2))
        kernel = kernel / sum(kernel)
        y_smoothed[x_idx] = sum(y * kernel)
    y_smoothed = y_smoothed / max(y_smoothed)

    if meta.save_smooth_plot:
        plots.smooth(meta, x, y, y_smoothed)

    return (x, y_smoothed)


#s10

def di_reformat(meta):
    """
    This function was introduced because some observations have several DIs per orbit. The user can set in the pcf how they want to determine the DI target position in this case.
    """
    iorbit_max = max(meta.iorbit_sp_cumulative)
    control_array = np.arange(iorbit_max + 1)

    reffile = ascii.read(meta.workdir + '/xrefyref.txt')
    #f = open(meta.workdir + "/xrefyref.txt", 'w')

    meta.refpix = np.zeros((iorbit_max+1, 3))

    # First case: Every orbit has just one DI
    if np.array_equal(reffile['iorbit_cumul'], control_array):
        print('There is one DI per orbit.')
        for i in range(iorbit_max + 1):
            meta.refpix[i] = [reffile['t_bjd'][i], reffile['pos1'][i], reffile['pos2'][i]]
            #print(reffile['t_bjd'][i], reffile['pos1'][i], reffile['pos2'][i], file=f)
        #f.close()
    # Second case: Every orbit contains at least one DI. But there is at least one orbit with more than one DI.
    elif set(control_array) == set(reffile['iorbit_cumul']):
        print('There is at least one orbit with at least more than one DI 1')
        for i in range(iorbit_max + 1):
            mask_i = reffile['iorbit'] == i
            if meta.di_multi == 'median':
                meta.refpix[i] = [np.median(reffile['t_bjd'][mask_i]), np.median(reffile['pos1'][mask_i]), np.median(reffile['pos2'][mask_i])]
                #print(np.median(reffile['t_bjd'][mask_i]), np.median(reffile['pos1'][mask_i]), np.median(reffile['pos2'][mask_i]), file=f)
            elif meta.di_multi == 'latest':
                meta.refpix[i] = [reffile['t_bjd'][mask_i][-1], reffile['pos1'][mask_i][-1], reffile['pos2'][mask_i][-1]]
                #print(reffile['t_bjd'][mask_i][-1], reffile['pos1'][mask_i][-1], reffile['pos2'][mask_i][-1], file=f)
        #f.close()

    # TODO this here
    # Third case. Not every orbit has a DI.
    if set(control_array) != set(reffile['iorbit_cumul']) and len(set(control_array)) < len(
            set(reffile['iorbit_cumul'])):
        print('There is at least one orbit with at least more than one DI 2')
        #f.close()

#s20

def get_wave_grid(meta):
    """
    Gets grid of wavelength solutions for each orbit and row.
    """
    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry141 as geo

    wave_grid = np.empty((meta.norbit*meta.nvisit, meta.subarray_size, meta.subarray_size))

    #calculates wavelength solution row by row for each orbit
    for i in range(meta.norbit):
        for j in range(meta.subarray_size):
            disp_solution = geo.dispersion(meta.refpix[i,1], -meta.LTV2+j)
            delx = 0.5 + np.arange(meta.subarray_size) - (meta.refpix[i,2] + meta.LTV1 + meta.POSTARG1/meta.platescale)
            wave_grid[i, j, :] = disp_solution[0] + delx*disp_solution[1]

    return wave_grid


def get_flatfield(meta):                    #function that flatfields a data array D, which starts at [minr, cmin] of hdu[1].data
    """
    Opens the flat file and uses it for bad pixel masking.
    """
    flat = fits.open(meta.flat)                #reads in flatfield cube
    WMIN = flat[0].header['WMIN']                #constants for computing flatfield coefficients
    WMAX = flat[0].header['WMAX']
    #print('flatfield:', WMIN, WMAX)
    a0 = flat[0].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]
    a1 = flat[1].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]
    a2 = flat[2].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]

    flatfield = []
    for i in range(meta.norbit*meta.nvisit):
        wave = meta.wave_grid[i, :]
        x = (wave - WMIN)/(WMAX-WMIN)
        flatfield.append(a0+a1*x+a2*x**2)
        flatfield[i][flatfield[i] < 0.5] = -1.        #sets flatfield pixels below 0.5 to -1 so they can be masked
    return flatfield


def median_abs_dev(vec):
    """
    Used to determine the variance for the background count estimate
    """
    med = ma.median(vec)
    return ma.median(abs(vec - med))


def read_refspec(meta):
    #https://matthew-brett.github.io/teaching/smoothing_intro.html
    refspec = np.loadtxt(meta.workdir + '/ancil/refspec/refspec.txt').T
    x_refspec, y_refspec = refspec[0]*1e10, refspec[1]/max(refspec[1])

    return (x_refspec, y_refspec)


def residuals2(params, x1, y1, x2, y2):
    """
    calculate residuals for leastsq.
    """
    a, b, c = params
    x1=np.array(x1)
    x2=np.array(x2)
    y1=np.array(y1)
    y2=np.array(y2)

    f = interp1d(x1, y1, kind='cubic')
    fit = f(a+b*x2)*c

    return fit - y2


def correct_wave_shift_fct_0(meta, orbnum, cmin, cmax, spec_opt, i):
    template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]

    g102mask = template_waves > 8200  # we dont use the spectrum below 8200 angstrom for the interpolation as the reference bandpass cuts out below this wavelength

    x_refspec, y_refspec = read_refspec(meta)

    # TODO: This is so bad
    x_refspec_new = np.concatenate((np.linspace(-5000, min(x_refspec), 10, endpoint=False),
                                    x_refspec,
                                    np.linspace(max(x_refspec) + 350, 30000, 10, endpoint=False)))
    y_refspec_new = np.concatenate((np.zeros(10),
                                    y_refspec,
                                    np.zeros(10)))

    # TODO: will break if optimal extractions isnt used!
    x_data = template_waves[g102mask]
    y_data = (spec_opt / max(spec_opt))[g102mask]

    p0 = [0, 1, 1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2, p0, args=(x_refspec_new, y_refspec_new, x_data, y_data))[0]

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit(x_refspec_new, y_refspec_new, p0, x_data, y_data, leastsq_res, meta, i)

    # for all other but first exposure in visit exposures
    x_data_firstexpvisit = leastsq_res[0] + template_waves * leastsq_res[1]
    y_data_firstexpvisit = np.copy(y_data)

    # np.savetxt('testing_exp0.txt', list(zip(x_data_firstexpvisit, y_data_firstexpvisit)))
    # np.savetxt('testing_firstexp.txt', list(zip(template_waves, spec_opt/max(spec_opt))))


    return x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res


def correct_wave_shift_fct_1(meta, orbnum, cmin, cmax, spec_opt, x_data_firstexpvisit, y_data_firstexpvisit, i):
    # TODO: So bad too
    x_model = np.concatenate((np.linspace(-5000, min(x_data_firstexpvisit), 10, endpoint=False),
                              x_data_firstexpvisit,
                              np.linspace(max(x_data_firstexpvisit) + 350, 30000, 10, endpoint=False)))
    y_model = np.concatenate((np.zeros(10),
                              y_data_firstexpvisit,
                              np.zeros(10)))

    x_data = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]
    y_data = spec_opt / max(spec_opt)

    p0 = [0, 1, 1]
    leastsq_res = leastsq(residuals2, p0, args=(x_model, y_model, x_data, y_data))[0]

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit(x_model, y_model, p0, x_data, y_data, leastsq_res, meta, i)

    wvls = leastsq_res[0] + x_data * leastsq_res[1]

    return wvls, leastsq_res


def make_lsq_rprs_txt(vals, errs, idxs, meta):
    """
    Saves the rprs vs wvl as a txt file as resulting from the lsq.
    """
    if not os.path.isdir(meta.workdir + meta.fitdir + '/lsq_res'):
        os.makedirs(meta.workdir + meta.fitdir + '/lsq_res')
    f_lsq = open(meta.workdir + meta.fitdir + '/lsq_res' + "/lsq_rprs.txt", 'w')
    rprs_vals_lsq = [vals[ii][idxs[0][1]] for ii in range(len(vals))]
    rprs_errs_lsq = [errs[ii][idxs[0][1]] for ii in range(len(errs))]
    file_header = ['wavelength (micron)', 'rprs', 'rprs_err', 'chi2red']
    print("#{: <24} {: <25} {: <25} {: <25}".format(*file_header), file=f_lsq)
    for row in zip(meta.wavelength_list, rprs_vals_lsq, rprs_errs_lsq, meta.chi2red_list):
        print("{: <25} {: <25} {: <25} {: <25}".format(*row), file=f_lsq)
    f_lsq.close()


def make_mcmc_rprs_txt(meta):
    """
    Saves the rprs vs wvl as a txt file as resulting from the MCMC.
    """
    files_mcmc_res = glob.glob(os.path.join(meta.workdir + meta.fitdir + '/mcmc_res', "mcmc_out_*.p"))
    files_mcmc_res = sn(files_mcmc_res)

    medians = []
    errors_lower = []
    errors_upper = []

    for f in files_mcmc_res:
        handle = open(f, 'rb')
        data, params, chain = pickle.load(handle)
        ndim = chain.shape[-1]
        samples = chain[:, meta.run_nburn:, :].reshape((-1, ndim))
        # TODO samples[:, 1] has to be fixeD!!!!!
        q = quantile(samples[:, 1], [0.16, 0.5, 0.84])
        medians.append(q[1])
        errors_lower.append(abs(q[1] - q[0]))
        errors_upper.append(abs(q[2] - q[1]))

    medians = np.array(medians)
    errors_lower = np.array(errors_lower)
    errors_upper = np.array(errors_upper)

    if not os.path.isdir(meta.workdir + meta.fitdir + '/mcmc_res'):
        os.makedirs(meta.workdir + meta.fitdir + '/mcmc_res')
    f_mcmc = open(meta.workdir + meta.fitdir + '/mcmc_res' + "/mcmc_rprs.txt", 'w')
    file_header = ['wavelength (micron)', 'rprs', 'rprs_err_lower', 'rprs_err_upper']
    print("#{: <24} {: <25} {: <25} {: <25}".format(*file_header), file=f_lsq)
    for row in zip(meta.wavelength_list, medians, errors_lower, errors_upper):
        print("{: <25} {: <25} {: <25} {: <25}".format(*row), file=f_lsq)
    f_lsq.close()


def quantile(x, q):
    return np.percentile(x, [100. * qi for qi in q])


def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
    ind = err != 0.0
    weights = 1.0/err[ind]**2
    mu = np.sum(data[ind]*weights)/np.sum(weights)
    var = 1.0/np.sum(weights)
    return [mu, np.sqrt(var)]                




def residuals(params, template_waves, template, spectrum, error):
    shift, scale = params
    fit = scale*np.interp(template_waves, template_waves-shift, spectrum)
    x = (template-fit)/error
    return (template-fit)/error




def interpolate_spectrum(spectrum, error, template, template_waves):
    p0 = [1., 1.0]                                        #initial guess for parameters shift and scale
    plsq, success  = leastsq(residuals, p0, args=(template_waves, template, spectrum, error))
    shift, scale = plsq
    interp_spectrum = np.interp(template_waves, template_waves-shift, spectrum)
    interp_error = np.interp(template_waves, template_waves-shift, error)
    return [interp_spectrum, interp_error, shift]

# def read_dict(filename):
#     #with open(filename,"r") as text:
#     #    return dict(line.strip().split() for line in text)
#     dict = {}
#     with open(filename,"r") as text:
#         for line in text:
#             key, value = line.split()
#             value = str_to_num(value)
#             if value == 'True': value = True
#             if value == 'False': value = False
#             dict[key] = value
#     return dict



#def getphase(t):
#    phase = (t - t0)/period
#    return phase - int(phase)









