import numpy as np
import yaml

obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)

if obs_par['GRISM'] == 'G141':
    grism = 'g141_throughput.txt'
elif obs_par['GRISM'] == 'G102':
    grism = 'g102_throughput.txt'

Teff, logg, MH = obs_par['Teff'], obs_par['logg'], obs_par['MH']

n_bins=obs_par['wvl_bins']
wave_min=obs_par['wvl_min']*1e4
wave_max=obs_par['wvl_max']*1e4

wvl_bins = np.linspace(wave_min, wave_max, n_bins)


f = open('./limb-darkening/input_files/ld_inputfile.txt', 'w')

for i in range(n_bins-1):
    params = [i, Teff, logg, MH, '-1', grism, 'P100', wvl_bins[i], wvl_bins[i+1]]
    params_bin = '\t'.join(map(str,params))
    print(params_bin)
    print(params_bin, file = f)
f.close()
