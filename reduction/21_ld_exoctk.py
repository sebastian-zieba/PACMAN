import os
from bokeh.plotting import figure
from bokeh.io import output_notebook
from exoctk import modelgrid
from svo_filters import svo
from exoctk.limb_darkening import limb_darkening_fit as lf
import astropy.units as q

import yaml
obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)

if obs_par['GRISM'] == 'G141':
    grism = 'WFC3_IR.G141'
elif obs_par['GRISM'] == 'G102':
    grism = 'WFC3_IR.G102'

n_bins=obs_par['wvl_bins']
wave_min=obs_par['wvl_min']*q.um
wave_max=obs_par['wvl_max']*q.um

Teff, logg, MH = obs_par['Teff'], obs_par['logg'], obs_par['MH']

G141 = svo.Filter(grism, n_bins=n_bins, wave_min=wave_min, wave_max=wave_max)
G141.plot()
#
#The ATLAS9 includes models for abundances [M/H]=0.0, -0.5, -1.0, -1.5, -2.0, -2.5, +0.5, +0.2 and gravity range from log g= 0.0 to +5.0 in steps of +0.5. The range in effective temperature from 3500 K to 50000 K 
#https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas.html

model_grid = modelgrid.ACES(resolution=700)
model_grid.customize(Teff_rng=(3000,5000), logg_rng=(4,5), FeH_rng=(-0.5,0.5))
ld = lf.LDC(model_grid)
ld.calculate(Teff, logg, MH,  'quadratic', bandpass=G141)
print(ld.results)
ld.plot(show=True)
