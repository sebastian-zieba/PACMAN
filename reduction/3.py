# Based on a perl script found on https://renenyffenegger.ch/notes/Wissenschaft/Astronomie/Ephemeriden/JPL-Horizons
# Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi (see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi)
# Also helpful: https://github.com/kevin218/POET/blob/master/code/doc/spitzer_Horizons_README.txt

import os
import numpy as np
from astropy.io import ascii
#from astropy.time import Time
import urllib.request 

filelist_path = './config/filelist.txt'
data = ascii.read(filelist_path)

nvisit = data['nvisit']
times = data['times']


settings = [
	"COMMAND= -48", #Hubble
	"CENTER= 500@0", #Solar System Barycenter (SSB) [500@0]
	"MAKE_EPHEM= YES",
	"TABLE_TYPE= VECTORS",
	#"START_TIME= $ARGV[0]",
	#"STOP_TIME= $ARGV[1]",
	"STEP_SIZE= 5m", # 5 Minute interval 
	"OUT_UNITS= KM-S",
	"REF_PLANE= FRAME",
	"REF_SYSTEM= J2000",
	"VECT_CORR= NONE",
	"VEC_LABELS= YES",
	"VEC_DELTA_T= NO",
	"CSV_FORMAT= NO",
	"OBJ_DATA= YES",
	"VEC_TABLE= 3"]

for i, setting in enumerate(settings):
    settings[i] = settings[i].replace(" =", "=").replace("= ", "=")
    settings[i] = settings[i].replace(" ", "%20")
    settings[i] = settings[i].replace("&", "%26")
    settings[i] = settings[i].replace(";", "%3B")
    settings[i] = settings[i].replace("?", "%3F")

settings = '&'.join(settings)
settings = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&' + settings
#print(settings)

for i in range(max(nvisit)+1): 
    print('Retrieving Horizons file for visit {0}/{1}'.format(i, max(nvisit)))
    times_visit = times[np.where(nvisit == i)]
    t_start = min(times_visit) + 2400000.5 - 1/24 #Start of Horizons file one hour before first exposure in visit
    t_end = max(times_visit) + 2400000.5 + 1/24 #End of Horizons file one hour after last exposure in visit

    
    set_start = "START_TIME=JD{0}".format(t_start)
    set_end = "STOP_TIME=JD{0}".format(t_end)

    settings_new = settings + '&' + set_start + '&' + set_end

    #print(settings_new )

    #print(t_start)
    #print(t_end)

    #t_start = Time(t_start, format='jd', scale='utc')
    #t_end = Time(t_end, format='jd', scale='utc')
    #print(t_start.isot)
    
    #os.system('perl ./util/horizons.pl JD{0} JD{1} > ./ancil/bjd_conversion/horizons_results_v{2}.txt'.format(t_start, t_end, i))
    
    urllib.request.urlretrieve(settings_new, './ancil/bjd_conversion/horizons_results_v{0}.txt'.format(i))




