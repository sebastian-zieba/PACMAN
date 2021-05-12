from selenium import webdriver
from selenium.webdriver.support.select import Select 
import time

import yaml

obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)


sleeptime = 2

#https://sites.google.com/chromium.org/driver/
path = './util/chromedriver'

web = webdriver.Chrome(path)
web.get('https://exoctk.stsci.edu/limb_darkening')

time.sleep(sleeptime)


teff, logg, feh = str(obs_par['Teff']), str(obs_par['logg']), str(obs_par['MH'])
teff=3501 #FIXME temps lower than 3500 dont work online!

field_teff = web.find_element_by_id('teff')
field_logg = web.find_element_by_id('logg')
field_feh = web.find_element_by_id('feh')

field_teff.clear()
field_teff.send_keys(teff)
time.sleep(sleeptime)

field_logg.clear()
field_logg.send_keys(logg)
time.sleep(sleeptime)

field_feh.clear()
field_feh.send_keys(feh)
time.sleep(sleeptime)


#stellarmodel = web.find_element_by_id('aces').is_selected()
#if stellarmodel:
#    print('aces already selected')
#else:
#    web.find_element_by_id('aces').click()
#    print('aces now selected')



if obs_par['GRISM'] == 'G141':
    grism = 'WFC3_IR.G141'
elif obs_par['GRISM'] == 'G102':
    grism = 'WFC3_IR.G102'


select_bandpass = web.find_element_by_id("filterselect")
select_bandpass = Select(select_bandpass)
select_bandpass.select_by_value(grism)
time.sleep(sleeptime)




wave_min = str(obs_par['wvl_min'])
wave_max = str(obs_par['wvl_max'])
n_bins = str(obs_par['wvl_bins'])



field_wave_min = web.find_element_by_id('wave_min')
field_wave_max = web.find_element_by_id('wave_max')
field_n_bins = web.find_element_by_id('n_bins')

field_wave_min.clear()
field_wave_min.send_keys(wave_min)
time.sleep(sleeptime)

field_wave_max.clear()
field_wave_max.send_keys(wave_max)
time.sleep(sleeptime)

field_n_bins.clear()
field_n_bins.send_keys(n_bins)
time.sleep(sleeptime)

web.find_element_by_id('calculate_submit').click()

time.sleep(20)

r = web.find_elements_by_xpath ("/html/body/div[2]/div/table/tbody/tr")
c = web.find_elements_by_xpath ("/html/body/div[2]/div/table/tbody/tr[1]/td")

rc = len(r)
cc = len(c)

f = open('./ancil/LD_coefficients/ld_exoctk_online_output.txt', 'w')

header = ['lam_eff (microns)', 'lam_min (microns)', 'lam_max (microns)', 'c1', 'e1', 'c2', 'e2']
print('\t'.join(header), file = f)

for i in range(1, rc + 1) :
    curr_row = []
    for j in range(1, cc + 1) :        
        d = web.find_element_by_xpath("/html/body/div[2]/div/table/tbody/tr["+str(i)+"]/td["+str(j)+"]").text
        curr_row.append(d)
    print('\t'.join(curr_row), file = f)

f.close()




