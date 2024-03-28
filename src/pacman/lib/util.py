import math
from pathlib import Path
from typing import Optional

import numpy as np
import numpy.ma as ma
from astropy.io import ascii, fits
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from scipy.signal import find_peaks

from .options import OPTIONS
from .sort_nicely import sort_nicely as sn
from ..lib import plots


# s00
def readfiles(meta):
    """Reads in the files saved in datadir and saves them into a list.

    Parameters
    ----------
    meta :
        Metadata object'

    Returns
    -------
    meta :
        Metadata object but adds segment_list to metadata containing the
        sorted data fits files.

    Notes
    -----
    History:
        Written by Sebastian Zieba      December 2021
    """
    meta.segment_list = []
    for fname in meta.datadir.iterdir():
        if fname.name.endswith(f'{meta.suffix}.fits'):
            meta.segment_list.append(meta.datadir / fname)
    return meta


# s02
def ancil(meta, s10: Optional[bool] = False, s20: Optional[bool] = False):
    """This function loads in a lot of useful arrays and values into meta.

    The following additional information are being loading into meta:

    - **norbit:**   number of orbits
    - **nvisit:**   number of visits
    - **files_sp:** all spectra files
    - **files_di:** all direct image files
    - **ra:** RA of the target in radians (from the header) (Note: this data is taken from the first spectrum file)
    - **dec:** DEC of the target in radians (from the header) (Note: this data is taken from the first spectrum file)
    - **coordtable:** a list of files containing the vector information of HST downloaded in s01


    Parameters
    ----------
    meta
        metadata object
    s10: bool
        Is set to True when s10 is being performed
    s20: bool
        Is set to True when s20 is being performed

    Returns
    -------
    meta
        metadata object

    Notes
    -----
    History:
        Written by Sebastian Zieba      December 2021
    """
    filelist = ascii.read(meta.workdir / 'filelist.txt')

    aa = filelist['iorbit'].data
    bb = np.diff(aa)
    meta.norbit, meta.nvisit = len(np.insert(aa[1:][np.where(bb != 0)], 0, aa[0])), len(set(filelist['ivisit'].data))

    meta.mask_sp = np.array([i[0] for i in filelist['instr']]) == 'G'
    meta.mask_di = np.array([i[0] for i in filelist['instr']]) == 'F'
    meta.files_sp = [Path(meta.datadir) / i for i in filelist['filenames'][meta.mask_sp].data]
    meta.files_di = [Path(meta.datadir) / i for i in filelist['filenames'][meta.mask_di].data]

    f = fits.open(meta.files_sp[0])
    meta.ra = f[0].header['ra_targ'] * math.pi / 180.0  # stores right ascension
    meta.dec = f[0].header['dec_targ'] * math.pi / 180.0  # stores declination

    meta.coordtable = []  # table of spacecraft coordinates
    for i in range(max(filelist['ivisit']) + 1): meta.coordtable.append(
        meta.workdir / 'ancil' / 'horizons' / f'horizons_results_v{i}.txt')

    # 03
    # TODO: Q: Is it okay if I assume everything in datadir uses same Filter Grism pair?
    meta.filter = filelist['instr'][meta.mask_di][0]
    meta.grism = filelist['instr'][meta.mask_sp][0]

    meta.scans_sp = filelist['scan'][meta.mask_sp].data
    meta.iorbit_sp = filelist['iorbit'][meta.mask_sp].data
    meta.ivisit_sp = filelist['ivisit'][meta.mask_sp].data
    meta.iexp_orb_sp = filelist['iexp_orb'][meta.mask_sp].data

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

    # TODO: this is a problem if there is one DI per visit. It will give us meta.iorbit_di_cumulative == [0,0,0,...]
    # TODO: with as many zeros as visits
    # TODO: I implemented a quick fix (29th march 2023):
    # TODO: i check if there are as many DIs as visits and if the DIs have been taken same time apart.
    t_bjd_di_tmp = filelist['t_mjd'][meta.mask_di].data
    if (len(meta.iorbit_di) == (max(meta.ivisit_sp)+1)) and (np.all(np.diff(t_bjd_di_tmp)>1/24)):
        meta.iorbit_di_cumulative = np.arange(len(meta.iorbit_di), dtype=int)
    meta.t_mjd_sp = filelist['t_mjd'][meta.mask_sp].data

    meta.t_orbit_sp = filelist['t_orbit'][meta.mask_sp].data
    meta.t_visit_sp = filelist['t_visit'][meta.mask_sp].data
    meta.platescale = 0.13  # IR detector has plate scale of 0.13 arcsec/pixel

    f_sp0 = fits.open(meta.files_sp[0])
    meta.POSTARG1_sp = f_sp0[0].header['POSTARG1']  # x-coordinate of the observer requested target offset
    meta.POSTARG2_sp = f_sp0[0].header['POSTARG2']  # y-coordinate of the observer requested target offset
    meta.LTV1 = int(f_sp0[1].header['LTV1'])     # X offset to get into physical pixels
    meta.LTV2 = int(f_sp0[1].header['LTV2'])     # Y offset to get into physical pixels
    meta.subarray_size = f_sp0[1].header['SIZAXIS1']  # size of subarray
    f_sp0.close()

    f_di0 = fits.open(meta.files_di[0])
    meta.POSTARG1_di = f_di0[0].header['POSTARG1']  # x-coordinate of the observer requested target offset
    meta.POSTARG2_di = f_di0[0].header['POSTARG2']  # y-coordinate of the observer requested target offset

    meta.POSTARG1 = meta.POSTARG1_sp - meta.POSTARG1_di
    meta.POSTARG2 = meta.POSTARG2_sp - meta.POSTARG2_di

    if meta.grism == 'G102':
        meta.BEAMA_i = 41
        meta.BEAMA_f = 248
    elif meta.grism == 'G141':
        meta.BEAMA_i = 15
        meta.BEAMA_f = 196
    f_di0.close()

    if s10:
        if 't_bjd' in filelist.keys():
            meta.t_bjd_sp = filelist['t_bjd'][meta.mask_sp].data
        if 't_bjd' in filelist.keys():
            meta.t_bjd_di = filelist['t_bjd'][meta.mask_di].data

    if s20:
        meta.rdnoise = 22.0  # read noise
        if meta.grism == 'G102':
            meta.flat = meta.pacmandir / 'data' / 'flats' / 'WFC3.IR.G102.flat.2.fits'
        elif meta.grism == 'G141':
            meta.flat = meta.pacmandir / 'data' / 'flats' / 'WFC3.IR.G141.flat.2.fits'
    return meta

# 03
def gaussian_kernel(meta, x, y):
    """Performs a gaussian kernel over an array. Used in smoothing of the stellar spectrum.
    Taken from: https://stackoverflow.com/questions/24143320/gaussian-sum-filter-for-irregular-spaced-points.
    """
    y = y / max(y)
    x.sort()

    resolution = min(len(x), 3000)

    x_eval = np.linspace(min(x), max(x), resolution)
    sigma = meta.smooth_sigma * 1e-10 / 5

    delta_x = x_eval[:, None] - x
    weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
    weights /= np.sum(weights, axis=1, keepdims=True)
    y_eval = np.dot(weights, y)

    y_eval = y_eval / max(y_eval)

    if meta.save_smooth_plot:
        plots.smooth(meta, x, y, x_eval, y_eval)
    return x_eval, y_eval


# def gaussian_kernel_old(meta, x, y):
#     """
#     https://matthew-brett.github.io/teaching/smoothing_intro.html
#     """
#
#     sigma = meta.smooth_sigma * 1e-10
#
#     y = y / max(y)
#     y_smoothed = np.zeros(y.shape)
#     for x_idx, x_val in tqdm(enumerate(x), desc='smoothing...', total=len(x)):
#         kernel = np.exp(-(x - x_val) ** 2 / (2 * sigma ** 2))
#         kernel = kernel / sum(kernel)
#         y_smoothed[x_idx] = sum(y * kernel)
#     y_smoothed = y_smoothed / max(y_smoothed)
#
#     if meta.save_smooth_plot:
#         plots.smooth(meta, x, y, x, y_smoothed)
#
#     return (x, y_smoothed)


# s10
def di_reformat(meta):
    """This function was introduced because some observations have several DIs per orbit.
    The user can set in the pcf how they want to determine the DI target position in this case.
    """
    iorbit_max = max(meta.iorbit_sp_cumulative)
    ivisit_max = max(meta.ivisit_sp)
    control_one_per_orbit = np.arange(iorbit_max + 1)
    control_one_per_visit = np.arange(ivisit_max + 1)
    reffile = ascii.read(meta.workdir / 'xrefyref.txt')

    meta.ivisits_new_orbit = meta.iorbit_sp[meta.new_visit_idx_sp]
    meta.nvisits_in_orbit = np.diff(np.append(meta.iorbit_sp[meta.new_visit_idx_sp], np.array([max(meta.iorbit_sp) + 1])))
    meta.norbits_in_visit = np.append(meta.iorbit_sp[meta.new_visit_idx_sp[1:]-1],meta.iorbit_sp[-1])+1

    # First case: Every orbit has just one DI
    if np.array_equal(reffile['iorbit_cumul'], control_one_per_orbit):
        print('There is one DI per orbit.')
        meta.refpix = np.zeros((iorbit_max + 1, 3))
        for i in range(iorbit_max + 1):
            meta.refpix[i] = [reffile['t_bjd'][i], reffile['pos1'][i], reffile['pos2'][i]]
            # print(reffile['t_bjd'][i], reffile['pos1'][i], reffile['pos2'][i], file=f)
        # f.close()

    # Second case: Every visit has just one DI
    elif np.array_equal(reffile['iorbit_cumul'], control_one_per_visit):
        print('There is one DI per visit.')
        meta.refpix = []
        counter = 0
        for i in range(len(meta.norbits_in_visit)):
            for j in range(meta.norbits_in_visit[counter]):
                meta.refpix.append([reffile['t_bjd'][counter], reffile['pos1'][counter], reffile['pos2'][counter]])
            counter = counter + 1
        meta.refpix = np.array(meta.refpix)

    # Third case: Every orbit contains at least one DI. But there is at least one orbit with more than one DI.
    elif set(control_one_per_orbit) == set(reffile['iorbit_cumul']):
        print('There is at least one orbit with at least more than one DI')
        meta.refpix = np.zeros((iorbit_max + 1, 3))
        for i in range(iorbit_max + 1):
            mask_i = reffile['iorbit'] == i
            if meta.di_multi == 'median':
                meta.refpix[i] = [np.median(reffile['t_bjd'][mask_i]), np.median(reffile['pos1'][mask_i]), np.median(reffile['pos2'][mask_i])]
                # print(np.median(reffile['t_bjd'][mask_i]), np.median(reffile['pos1'][mask_i]), np.median(reffile['pos2'][mask_i]), file=f)
            elif meta.di_multi == 'latest':
                meta.refpix[i] = [reffile['t_bjd'][mask_i][-1], reffile['pos1'][mask_i][-1], reffile['pos2'][mask_i][-1]]
                # print(reffile['t_bjd'][mask_i][-1], reffile['pos1'][mask_i][-1], reffile['pos2'][mask_i][-1], file=f)
        # f.close()

    # TODO Add this check
    # Third case. Not every orbit has a DI.
    if set(control_one_per_orbit) != set(reffile['iorbit_cumul']) and len(set(control_one_per_orbit)) < len(
            set(reffile['iorbit_cumul'])):
        print('Not every orbit has a DI.')
        # f.close()

# s20
def get_wave_grid(meta):
    """Gets grid of wavelength solutions for each orbit and row."""
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


def get_flatfield(meta):                    # function that flatfields a data array D, which starts at [minr, cmin] of hdu[1].data
    """Opens the flat file and uses it for bad pixel masking."""
    flat = fits.open(meta.flat)             # reads in flatfield cube
    WMIN = flat[0].header['WMIN']           # constants for computing flatfield coefficients
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


def peak_finder(array, i, ii, orbnum, meta):
    """Finding peaks in an array."""
    # determine the aperture cutout
    rowmedian = np.median(array, axis=1)  # median of every row
    rowmedian_absder = abs(
        rowmedian[1:] - rowmedian[:-1])  # absolute derivative in order to determine where the spectrum is
    # we use scipy.signal.find_peaks to determine in which rows the spectrum starts and ends
    # TODO: Think about better values for height and distance
    # TODO: Finding the row with the highest change in flux (compared to row above and below) isnt robust against outliers!
    peaks, _ = find_peaks(rowmedian_absder, height=max(rowmedian_absder * 0.3), distance=2)
    # peaks = peaks[:2] # only take the two biggest peaks (there should be only 2)
    peaks = np.array([min(peaks[:2]), max(peaks[:2]) + 1])

    if meta.save_utr_plot or meta.show_utr_plot:
        plots.utr(array, meta, i, ii, orbnum, rowmedian, rowmedian_absder, peaks)
    return peaks


def median_abs_dev(vec):
    """Used to determine the variance for the background count estimate."""
    med = ma.median(vec)
    return ma.median(abs(vec - med))


def zero_pad_x(array):
    # TODO: Try to make this look nicer
    array_new = np.concatenate((np.linspace(-5000, min(array)+1, 10, endpoint=True),
                                array, np.linspace(max(array) + 350, 30000, 10, endpoint=False)))
    return array_new


def zero_pad_y(array):
    return np.concatenate((np.zeros(10), array, np.zeros(10)))


def read_refspec(meta):
    """Reads in the reference spectrum."""
    refspec = np.loadtxt(meta.workdir / 'ancil' / 'refspec' / 'refspec.txt').T
    x_refspec, y_refspec = refspec[0]*1e10, refspec[1]/max(refspec[1])

    x_refspec = zero_pad_x(x_refspec)
    y_refspec = zero_pad_y(y_refspec)
    return x_refspec, y_refspec


def residuals2(params, x1, y1, x2, y2):
    """Calculate residuals for leastsq."""
    a, b, c = params
    x1 = np.array(x1)
    x2 = np.array(x2)
    y1 = np.array(y1)
    y2 = np.array(y2)

    f = interp1d(x1, y1, kind='linear')
    fit = f(a+b*x2)*c
    return fit - y2


def residuals2_lin(params, x1, y1, x2, y2):
    """Calculate residuals for leastsq."""
    a, c = params
    x1 = np.array(x1)
    x2 = np.array(x2)
    y1 = np.array(y1)
    y2 = np.array(y2)

    f = interp1d(x1, y1, kind='linear')
    fit = f(a+x2)*c
    return fit - y2

def get_boxspat(meta,i,d,rmin,rmax,cmin,cmax,imageIndex):#,x_ref, y_ref,x_data, y_data,i,direction='y'):
    fullframe_first = d[imageIndex].data[rmin:rmax,cmin:cmax] #first non-destructive exposure
                
    # NOTE: Background subtraction
    if not meta.background_box:
        below_threshold = fullframe_first < meta.background_thld # mask with all pixels below the user defined threshold
        skymedian_first = np.median(fullframe_first[below_threshold].flatten())  # estimates the background counts by taking the flux median of the pixels below the flux threshold
    else:
        bg_rmin, bg_rmax, bg_cmin, bg_cmax = int(meta.bg_rmin), int(meta.bg_rmax), int(meta.bg_cmin), int(meta.bg_cmax),
        skymedian_first = np.median(fullframe_first[bg_rmin:bg_rmax, bg_cmin:bg_cmax].flatten())

    spectrum_first = fullframe_first - skymedian_first
    spat_box_0_first = spectrum_first.sum(axis = 1)                            #box-extracted spatial profile, sums along dispersion direction (sums columns in a row)   
    
    # NOTE: norm profile and remove nans.
    r_array = np.linspace(rmin,rmax+1,rmax-rmin) #array with row indexes
    spat_box_0_first = spat_box_0_first/np.nanmax(spat_box_0_first)
    
    r_array = r_array[~np.isnan(spat_box_0_first)]
    spat_box_0_first = spat_box_0_first[~np.isnan(spat_box_0_first)]
    
    # NOTE: interpolate spectrum for fit
    r_array_interp = np.linspace(rmin,rmax+1,(rmax-rmin)*1000)
    spat_box_0_first_interp = np.interp(r_array_interp, r_array, spat_box_0_first)
    
    return r_array_interp, spat_box_0_first_interp


def calculate_rowshift(meta,x_ref, y_ref,x_data, y_data,i):
    p0 = [0,1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2_lin, p0, args=(x_ref, y_ref,x_data, y_data))[0]
    if meta.save_rowshift_plot or meta.show_rowshift_plot:
        plots.rowshift_fit(x_ref, y_ref, x_data, y_data, leastsq_res, meta,i)
    return leastsq_res[0]#float
    
def calculate_stretch(meta,x_ref, y_ref,x_data, y_data,i):
    #all np arrays f(a+b*x2)*c
    p0 = [0,1,1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2, p0, args=(x_ref, y_ref,x_data, y_data))[0]
    if meta.save_rowshift_stretch_plot or meta.show_rowshift_stretch_plot:
        plots.stretch_fit(x_ref, y_ref, x_data, y_data, leastsq_res, meta,i)
    return leastsq_res[1]#float

def correct_wave_shift_fct_0(meta, orbnum, cmin, cmax, spec_opt, i):
    """Use the reference spectrum for the wave cal."""
    template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]

    g102mask = template_waves > 820  # we dont use the spectrum below 8200 angstrom for the interpolation as the reference bandpass cuts out below this wavelength
    x_refspec, y_refspec = read_refspec(meta)

    # TODO: will break if optimal extractions isnt used!
    #x_data & y_data is the spectrum if no wavelength correction would be performed
    x_data = template_waves[g102mask]
    y_data = (spec_opt / max(spec_opt))[g102mask]

    p0 = [0, 1, 1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2, p0, args=(x_refspec, y_refspec, x_data, y_data))[0]
    #leastsq_res[0] = leastsq_res[0] + 60

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit(x_refspec, y_refspec, p0, x_data, y_data, leastsq_res, meta, i)

    # for all other but first exposure in visit exposures
    x_data_firstexpvisit = leastsq_res[0] + x_data * leastsq_res[1]
    y_data_firstexpvisit = np.copy(y_data)

    #np.savetxt('testing_template.txt', list(zip(x_refspec, y_refspec)))
    #np.savetxt('testing_exp0.txt', list(zip(x_data, y_data)))
    return x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res


def correct_wave_shift_fct_0_lin(meta, orbnum, cmin, cmax, spec_opt, i):
    """Use the reference spectrum for the wave cal."""
    template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]

    g102mask = template_waves > 820  # we dont use the spectrum below 8200 angstrom for the interpolation as the reference bandpass cuts out below this wavelength
    x_refspec, y_refspec = read_refspec(meta)

    # TODO: will break if optimal extractions isnt used!
    #x_data & y_data is the spectrum if no wavelength correction would be performed
    x_data = template_waves[g102mask]
    y_data = (spec_opt / max(spec_opt))[g102mask]

    p0 = [0, 1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2_lin, p0, args=(x_refspec, y_refspec, x_data, y_data))[0]
    # leastsq_res[0] = leastsq_res[0] + 60

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit_lin(x_refspec, y_refspec, p0, x_data, y_data, leastsq_res, meta, i)

    # for all other but first exposure in visit exposures
    x_data_firstexpvisit = leastsq_res[0] + x_data
    y_data_firstexpvisit = np.copy(y_data)

    # np.savetxt('/home/zieba/Desktop/Projects/Observations/Hubble/KELT11_15255/run_2022-05-25_16-57-51_new/spectra_drift/' + 'testing_template.txt', list(zip(x_refspec_new, y_refspec_new)))
    # np.savetxt('/home/zieba/Desktop/Projects/Observations/Hubble/KELT11_15255/run_2022-05-25_16-57-51_new/spectra_drift/' + 'testing_exp0.txt', list(zip(x_data, y_data)))
    return x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res


def correct_wave_shift_fct_00(meta, orbnum, cmin, cmax, spec_opt, i):
    """Use the first exposure in the visit for wave cal."""
    template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]

    x_refspec, y_refspec = template_waves, spec_opt / max(spec_opt)

    # TODO: will break if optimal extractions isnt used!
    #x_data & y_data is the spectrum if no wavelength correction would be performed
    x_data = template_waves
    y_data = (spec_opt / max(spec_opt))

    p0 = [0, 1, 1]  # initial guess for least squares
    leastsq_res = leastsq(residuals2, p0, args=(x_refspec, y_refspec, x_data, y_data))[0]

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit(x_refspec, y_refspec, p0, x_data, y_data, leastsq_res, meta, i)

    # for all other but first exposure in visit exposures
    x_data_firstexpvisit = leastsq_res[0] + x_data * leastsq_res[1]
    y_data_firstexpvisit = np.copy(y_data)

    # np.savetxt('testing_exp0.txt', list(zip(x_data_firstexpvisit, y_data_firstexpvisit)))
    # np.savetxt('testing_firstexp.txt', list(zip(template_waves, spec_opt/max(spec_opt))))
    return x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res


def correct_wave_shift_fct_1(meta, orbnum, cmin, cmax, spec_opt,
                             x_data_firstexpvisit, y_data_firstexpvisit, i):
    x_model = zero_pad_x(x_data_firstexpvisit)
    y_model = zero_pad_y(y_data_firstexpvisit)

    x_data = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]
    y_data = spec_opt / max(spec_opt)

    p0 = [0, 1, 1]
    leastsq_res = leastsq(residuals2, p0, args=(x_model, y_model, x_data, y_data))[0]

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit(x_model, y_model, p0, x_data, y_data, leastsq_res, meta, i)

    wvls = leastsq_res[0] + x_data * leastsq_res[1]
    #np.savetxt('testing_exp{0}.txt'.format(i), list(zip(x_data, y_data)))
    return wvls, leastsq_res


def correct_wave_shift_fct_1_lin(meta, orbnum, cmin, cmax, spec_opt,
                                 x_data_firstexpvisit, y_data_firstexpvisit, i):
    x_model = zero_pad_x(x_data_firstexpvisit)
    y_model = zero_pad_y(y_data_firstexpvisit)

    x_data = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]
    y_data = spec_opt / max(spec_opt)

    p0 = [0, 1]
    leastsq_res = leastsq(residuals2_lin, p0, args=(x_model, y_model, x_data, y_data))[0]

    if meta.save_refspec_fit_plot or meta.show_refspec_fit_plot:
        plots.refspec_fit_lin(x_model, y_model, p0, x_data, y_data, leastsq_res, meta, i)

    wvls = leastsq_res[0] + x_data
    # np.savetxt('/home/zieba/Desktop/Projects/Observations/Hubble/KELT11_15255/run_2022-05-25_16-57-51_new/spectra_drift/' + 'testing_exp{0}.txt'.format(i), list(zip(x_data, y_data)))
    return wvls, leastsq_res


# 30
def read_fitfiles(meta):
    """Read in the files (white or spectroscopic) which will be fitted."""
    if meta.s30_fit_white:
        print('White light curve fit will be performed')
        files = []
        if meta.s30_most_recent_s20:
            lst_dir = sn((meta.workdir / "extracted_lc").iterdir())
            # the following line makes sure that only directories starting with a "2" are considered
            # this was implemented after issue #10 was raised (see issue for more info)
            # this works because the dates will always start with a "2"
            lst_dir_new = [lst_dir_i for lst_dir_i in lst_dir if lst_dir_i.name.startswith("2")]
            white_dir = lst_dir_new[-1]
            files.append(meta.workdir / "extracted_lc" / white_dir / "lc_white.txt")
            print('using most recent s20 run: {0}'.format(white_dir))
        else:
            files.append(meta.s30_white_file_path)
            print('using user set white file: '.format(meta.s30_white_file_path))

        if meta.grism == 'G102':
            meta.wavelength_list = [1.0]
        elif meta.grism == 'G141':
            meta.wavelength_list = [1.4]

    elif meta.s30_fit_spec:
        print('Spectroscopic light curve fit(s) will be performed')
        if meta.s30_most_recent_s21:
            # find most recent bins directory
            lst_dir = np.array([path for path in (meta.workdir / "extracted_sp").iterdir()
                                if path.is_dir()])
            dirs_bool = np.array(['bins' in i.name for i in lst_dir])
            lst_dir = lst_dir[dirs_bool]
            dirs_times = [i.name[-19:] for i in lst_dir]  # There are 19 digits in the typical '%Y-%m-%d_%H-%M-%S' format
            # sort the times
            times_sorted = sn(dirs_times)
            # most recent time
            recent_time = times_sorted[-1]
            idx = 0
            for i in range(len(lst_dir)):
                if lst_dir[i].name[-19:] == recent_time:
                    idx = i
            spec_dir_full = lst_dir[idx]
            #spec_dir_full = meta.workdir + "/extracted_sp/" + spec_dir
            files = sn(spec_dir_full.glob("*.txt"))
            print(f'using most recent s21 run: {spec_dir_full}')
        else:
            spec_dir_full = meta.s30_spec_dir_path
            files = sn(spec_dir_full.glob("*.txt"))

            print('using user set spectroscopic directory: '.format(meta.s30_spec_dir_path))
        spec_dir_wvl_file = spec_dir_full / 'wvl_table.dat'
        meta.wavelength_list = ascii.read(spec_dir_wvl_file)['wavelengths']
    else:
        print('Neither s30_fit_white nor s30_fit_spec are True!')

    print('Identified file(s) for fitting:', files)
    meta.nfits = len(files)
    return files, meta


def return_free_array(nvisit, fixed_array, tied_array):
    """Reads in the fit_par.txt and determines which parameters are free."""
    free_array = []
    for i in range(len(fixed_array)):
        if fixed_array[i].lower() == 'true' and tied_array[i] == -1:
            for ii in range(nvisit):
                free_array.append(False)
        if fixed_array[i].lower() == 'true' and not tied_array[i] == -1:
            free_array.append(False)
        if fixed_array[i].lower() == 'false' and tied_array[i] == -1:
            free_array.append(True)
            for ii in range(nvisit-1):
                free_array.append(False)
        if fixed_array[i].lower() == 'false' and not tied_array[i] == -1:
            free_array.append(True)
    free_array = np.array(free_array)
    return free_array

def return_untied_array(nvisit, fixed_array, tied_array):
    """Reads in the fit_par.txt and determines which parameters are untied."""
    untied_array = []
    for i in range(len(fixed_array)):
        if fixed_array[i].lower() == 'true' and tied_array[i] == -1:
            for ii in range(nvisit):
                untied_array.append(False)
        if fixed_array[i].lower() == 'true' and not tied_array[i] == -1:
            untied_array.append(True)
        if fixed_array[i].lower() == 'false' and tied_array[i] == -1:
            untied_array.append(False)
            for ii in range(nvisit-1):
                untied_array.append(False)
        if fixed_array[i].lower() == 'false' and not tied_array[i] == -1:
            untied_array.append(True)
    untied_array = np.array(untied_array)
    return untied_array

def format_params_for_Model(theta, params, nvisit,
                            fixed_array, tied_array, free_array, untied_array):
    nvisit = int(nvisit)
    params_updated = []

    # #TODO: Oida?
    # def repeated(array, index, n_times):
    #     array_new = []
    #     index_index = 0
    #     for ii in range(len(array)):
    #         if ii in index:
    #             array_new.append([array[index[index_index]]] * n_times)
    #         else:
    #             array_new.append([array[ii]])
    #         index_index = index_index + 1
    #     array_new2 = np.array([item for sublist in array_new for item in sublist])
    #     return array_new2

    #new function that copies parameters of visit 0 in visits>0 if tied=-1.:
    theta_new = np.copy(theta)
    params_updated = np.copy(params)
    params_updated[free_array] = theta_new
    for i_p in range(0,len(untied_array), nvisit):
        if not untied_array[i_p]:
            for ii_p in range(1,nvisit):
                params_updated[i_p+ii_p] = params_updated[i_p]
    # print('free_array', free_array)
    # print('params_updated', params_updated)
    return np.array(params_updated)


def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
    """COMPUTE ROOT-MEAN-SQUARE AND STANDARD ERROR OF DATA FOR VARIOUS BIN SIZES
    Taken from POET: https://github.com/kevin218/POET/blob/master/code/lib/correlated_noise.py
    """
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


def format_params_for_sampling(params, meta, fit_par):
    # FIXME: make sure this works for cases when nvisit>1
    nvisit = int(meta.nvisit)

    fixed_array = np.array(fit_par['fixed'])
    tied_array = np.array(fit_par['tied'])
    free_array = []

    for i in range(len(fixed_array)):
        if fixed_array[i].lower() == 'true' and tied_array[i] == -1:
            for ii in range(nvisit):
                free_array.append(False)
        if fixed_array[i].lower() == 'true' and not tied_array[i] == -1:
            free_array.append(False)
        if fixed_array[i].lower() == 'false' and tied_array[i] == -1:
            free_array.append(True)
            for ii in range(nvisit-1):
                free_array.append(False)
        if fixed_array[i].lower() == 'false' and not tied_array[i] == -1:
            free_array.append(True)
    free_array = np.array(free_array)
    return np.array(params[free_array])


def make_lsq_rprs_txt(vals, errs, idxs, meta):
    """Saves the rprs vs wvl as a txt file as resulting from the lsq."""
    with (meta.workdir / meta.fitdir / 'lsq_res'/ "lsq_rprs.txt")\
            .open('w', encoding=OPTIONS["encoding"]) as f_lsq:
        rp_idx = np.where(np.array(meta.labels) == 'rp')[0][0]
        rprs_vals_lsq = [vals[ii][idxs[0][rp_idx]] for ii in range(len(vals))]
        rprs_errs_lsq = [errs[ii][idxs[0][rp_idx]] for ii in range(len(errs))]
        file_header = ['wavelength (micron)', 'rprs', 'rprs_err']
        print("#{: <24} {: <25} {: <25}".format(*file_header), file=f_lsq)
        for row in zip(meta.wavelength_list, rprs_vals_lsq, rprs_errs_lsq):
            print("{: <25} {: <25} {: <25}".format(*row), file=f_lsq)


def make_rprs_txt(vals, errs_lower, errs_upper, meta, fitter=None):
    """Saves the rprs vs wvl as a txt file as resulting from the sampler."""
    rp_idx = np.where(np.array(meta.labels) == 'rp')[0][0]
    medians = np.array(vals).T[rp_idx]
    errors_lower = np.array(errs_lower).T[rp_idx]
    errors_upper = np.array(errs_upper).T[rp_idx]
    if fitter == 'mcmc':
        with (meta.workdir / meta.fitdir / 'mcmc_res' / "mcmc_rprs.txt")\
                .open('w', encoding=OPTIONS["encoding"]) as f_mcmc:
            file_header = ['wavelength (micron)', 'rprs', 'rprs_err_lower', 'rprs_err_upper']
            print("#{: <24} {: <25} {: <25} {: <25}".format(*file_header), file=f_mcmc)
            for row in zip(meta.wavelength_list, medians, errors_lower, errors_upper):
                print("{: <25} {: <25} {: <25} {: <25}".format(*row), file=f_mcmc)
    if fitter == 'nested':
        with (meta.workdir / meta.fitdir / 'nested_res' / "nested_rprs.txt")\
                .open('w', encoding=OPTIONS["encoding"]) as f_nested:
            file_header = ['wavelength (micron)', 'rprs', 'rprs_err_lower', 'rprs_err_upper']
            print("#{: <24} {: <25} {: <25} {: <25}".format(*file_header), file=f_nested)
            for row in zip(meta.wavelength_list, medians, errors_lower, errors_upper):
                print("{: <25} {: <25} {: <25} {: <25}".format(*row), file=f_nested)


def quantile(x, q):
    return np.percentile(x, [100. * qi for qi in q])


def weighted_mean(data, err):
    """Calculates the weighted mean for data points data with std devs. err."""
    ind = err != 0.0
    weights = 1.0/err[ind]**2
    mu = np.sum(data[ind]*weights)/np.sum(weights)
    var = 1.0/np.sum(weights)
    return [mu, np.sqrt(var)]


def create_res_dir(meta):
    """Creates the result directory depending on which fitters were used."""
    #print('H11', meta.workdir)
    #print('H22', meta.fitdir)

    lsq_res_dir = meta.workdir / meta.fitdir / 'lsq_res'
    #print('H33', lsq_res_dir)
    #print('DOES IT EXIST', lsq_res_dir.exists())
    if not lsq_res_dir.exists():
        lsq_res_dir.mkdir(parents=True)
    #print('DOES THIS STILL RUN?')
    mcmc_res_dir = meta.workdir / meta.fitdir / 'mcmc_res'
    nested_res_dir = meta.workdir / meta.fitdir / 'nested_res'
    if meta.run_lsq:
        if not lsq_res_dir.exists():
            lsq_res_dir.mkdir(parents=True)
    if meta.run_mcmc:
        if not mcmc_res_dir.exists():
            mcmc_res_dir.mkdir(parents=True)
    if meta.run_nested:
        if not nested_res_dir.exists():
            nested_res_dir.mkdir(parents=True)


def log_run_setup(meta):
    """Prepares lists in meta where fit statistics will be saved into."""
    meta.rms_pred_list_lsq = []
    meta.rms_pred_list_mcmc = []
    meta.rms_pred_list_nested = []

    meta.rms_list_lsq = []
    meta.rms_list_mcmc = []
    meta.rms_list_nested = []

    meta.chi2_list_lsq = []
    meta.chi2_list_mcmc = []
    meta.chi2_list_nested = []

    meta.chi2red_list_lsq = []
    meta.chi2red_list_mcmc = []
    meta.chi2red_list_nested = []

    meta.bic_list_lsq = []
    meta.bic_list_mcmc = []
    meta.bic_list_nested = []

    # meta.bic_alt_list_lsq = []
    # meta.bic_alt_list_mcmc = []
    # meta.bic_alt_list_nested = []

    meta.ln_like_list_lsq = []
    meta.ln_like_list_mcmc = []
    meta.ln_like_list_nested = []

    if 'uncmulti' in meta.s30_myfuncs or meta.rescale_uncert:
        meta.chi2_notrescaled_list_lsq = []
        meta.chi2_notrescaled_list_mcmc = []
        meta.chi2_notrescaled_list_nested = []

        meta.chi2red_notrescaled_list_lsq = []
        meta.chi2red_notrescaled_list_mcmc = []
        meta.chi2red_notrescaled_list_nested = []

        meta.bic_notrescaled_list_lsq = []
        meta.bic_notrescaled_list_mcmc = []
        meta.bic_notrescaled_list_nested = []

        # meta.bic_alt_notrescaled_list_lsq = []
        # meta.bic_alt_notrescaled_list_mcmc = []
        # meta.bic_alt_notrescaled_list_nested = []

        meta.ln_like_notrescaled_list_lsq = []
        meta.ln_like_notrescaled_list_mcmc = []
        meta.ln_like_notrescaled_list_nested = []

        if 'uncmulti' in meta.s30_myfuncs:
            meta.uncmulti_mcmc = []
            meta.uncmulti_nested = []
    return meta


def append_fit_output(fit, meta, fitter=None, medians=None):
    """Appends fit statistics like rms or chi2 into meta lists."""
    precision = 3
    if fitter == 'lsq':
        meta.rms_pred_list_lsq.append(round(fit.rms_predicted, precision))
        meta.rms_list_lsq.append(round(fit.rms, precision))
        meta.chi2_list_lsq.append(round(fit.chi2, precision))
        meta.chi2red_list_lsq.append(round(fit.chi2red, precision))
        meta.bic_list_lsq.append(round(fit.bic, precision))
        # meta.bic_alt_list_lsq.append(round(fit.bic_alt, precision))
        meta.ln_like_list_lsq.append(round(fit.ln_like, precision))
    if fitter == 'mcmc':
        meta.rms_pred_list_mcmc.append(round(fit.rms_predicted, precision))
        meta.rms_list_mcmc.append(round(fit.rms, precision))
        meta.chi2_list_mcmc.append(round(fit.chi2, precision))
        meta.chi2red_list_mcmc.append(round(fit.chi2red, precision))
        meta.bic_list_mcmc.append(round(fit.bic, precision))
        # meta.bic_alt_list_mcmc.append(round(fit.bic_alt, precision))
        meta.ln_like_list_mcmc.append(round(fit.ln_like, precision))
    if fitter == 'nested':
        meta.rms_pred_list_nested.append(round(fit.rms_predicted, precision))
        meta.rms_list_nested.append(round(fit.rms, precision))
        meta.chi2_list_nested.append(round(fit.chi2, precision))
        meta.chi2red_list_nested.append(round(fit.chi2red, precision))
        meta.bic_list_nested.append(round(fit.bic, precision))
        # meta.bic_alt_list_nested.append(round(fit.bic_alt, precision))
        meta.ln_like_list_nested.append(round(fit.ln_like, precision))
    if ('uncmulti' in meta.s30_myfuncs) or (meta.rescale_uncert):
        if fitter == 'mcmc':
            meta.chi2_notrescaled_list_mcmc.append(round(fit.chi2_notrescaled, precision))
            meta.chi2red_notrescaled_list_mcmc.append(round(fit.chi2red_notrescaled, precision))
            meta.bic_notrescaled_list_mcmc.append(round(fit.bic_notrescaled, precision))
            # meta.bic_alt_notrescaled_list_mcmc.append(round(fit.bic_alt_notrescaled, precision))
            meta.ln_like_notrescaled_list_mcmc.append(round(fit.ln_like_notrescaled, precision))
            if ('uncmulti' in meta.s30_myfuncs):
                meta.uncmulti_mcmc.append(round(medians[-1], precision))
        if fitter == 'nested':
            meta.chi2_notrescaled_list_nested.append(round(fit.chi2_notrescaled, precision))
            meta.chi2red_notrescaled_list_nested.append(round(fit.chi2red_notrescaled, precision))
            meta.bic_notrescaled_list_nested.append(round(fit.bic_notrescaled, precision))
            # meta.bic_alt_notrescaled_list_nested.append(round(fit.bic_alt_notrescaled, precision))
            meta.ln_like_notrescaled_list_nested.append(round(fit.ln_like_notrescaled, precision))
            if ('uncmulti' in meta.s30_myfuncs):
                meta.uncmulti_nested.append(round(medians[-1], precision))


def save_fit_output(fit, data, meta):
    """Saves all the fit statistics like rms or chi2 into an astropy table."""
    t = Table()

    t['wave'] = meta.wavelength_list
    t['wave'].info.format = '.3f'
    if meta.run_mcmc:
        t['rms_pred'] = meta.rms_pred_list_mcmc
        t['rms'] = meta.rms_list_mcmc
        t['chi2'] = meta.chi2_list_mcmc
        t['chi2red'] = meta.chi2red_list_mcmc
        t['bic'] = meta.bic_list_mcmc
        # t['bic_alt']  = meta.bic_alt_list_mcmc
        t['ln_like'] = meta.ln_like_list_mcmc
        if ('uncmulti' in meta.s30_myfuncs) or (meta.rescale_uncert):
            t['chi2_notrescaled'] = meta.chi2_notrescaled_list_mcmc
            t['chi2red_notrescaled'] = meta.chi2red_notrescaled_list_mcmc
            t['bic_notrescaled'] = meta.bic_notrescaled_list_mcmc
            # t['bic_alt_notrescaled'] = meta.bic_alt_notrescaled_list_mcmc
            t['ln_like_notrescaled'] = meta.ln_like_notrescaled_list_mcmc
            if 'uncmulti' in meta.s30_myfuncs:
                t['uncmulti'] = meta.uncmulti_mcmc

    if meta.run_nested:
        t['rms_pred'] = meta.rms_pred_list_nested
        t['rms'] = meta.rms_list_nested
        t['chi2'] = meta.chi2_list_nested
        t['chi2red'] = meta.chi2red_list_nested
        t['bic'] = meta.bic_list_nested
        # t['bic_alt']  = meta.bic_alt_list_nested
        t['ln_like'] = meta.ln_like_list_nested
        if ('uncmulti' in meta.s30_myfuncs) or (meta.rescale_uncert):
            t['chi2_notrescaled'] = meta.chi2_notrescaled_list_nested
            t['chi2red_notrescaled'] = meta.chi2red_notrescaled_list_nested
            t['bic_notrescaled'] = meta.bic_notrescaled_list_nested
            # t['bic_alt_notrescaled'] = meta.bic_alt_notrescaled_list_nested
            t['ln_like_notrescaled'] = meta.ln_like_notrescaled_list_nested
            if 'uncmulti' in meta.s30_myfuncs:
                t['uncmulti'] = meta.uncmulti_nested

    if 'uncmulti' in meta.s30_myfuncs:
        t['uncmulti2'] = t['uncmulti'] ** 2
        t['uncmulti2'].info.format = '.3f'

    t['npoints'] = [data.npoints] * meta.nfits
    t['nfree_param'] = [data.nfree_param] * meta.nfits
    t['dof'] = [data.dof] * meta.nfits

    ascii.write(t, meta.workdir / meta.fitdir / 'fit_results.txt',
                format='rst', overwrite=True)
    print('Saved fit_results.txt file')


def save_allandata(binsz, rms, stderr, meta, fitter=None):
    """Saves the data used to create the Allan deviation plot."""
    t = Table()

    t['binsz'] = binsz
    t['rms'] = rms
    t['stderr'] = stderr

    fname = Path(f'{fitter}_res') /\
            f'corr_data_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt'
    ascii.write(t, meta.workdir / meta.fitdir / fname, format='rst', overwrite=True)


# def residuals(params, template_waves, template, spectrum, error):
#     shift, scale = params
#     fit = scale*np.interp(template_waves, template_waves-shift, spectrum)
#     x = (template-fit)/error
#     return (template-fit)/error


# def interpolate_spectrum(spectrum, error, template, template_waves):
#     p0 = [1., 1.0]                                        #initial guess for parameters shift and scale
#     plsq, success  = leastsq(residuals, p0, args=(template_waves, template, spectrum, error))
#     shift, scale = plsq
#     interp_spectrum = np.interp(template_waves, template_waves-shift, spectrum)
#     interp_error = np.interp(template_waves, template_waves-shift, error)
#     return [interp_spectrum, interp_error, shift]


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


# def getphase(t):
#    phase = (t - t0)/period
#    return phase - int(phase)
