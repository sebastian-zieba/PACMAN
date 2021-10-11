import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload
from astropy.io import ascii

import wfc3_reduction.reduction.s30_run as s30

from wfc3_reduction.lib.reload_meta import reload_meta
from wfc3_reduction.lib import manageevent as me

eventlabel = 'L-98-59_Hubble15856'

#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-19_19-22-39_L-98-59_Hubble15856/' #[all], window=10, correct_shift=False, corrected refpix
workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-21_02-00-29_L-98-59_Hubble15856/'

reload_meta(workdir, eventlabel)
meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

nvisit = meta.nvisit


fit_par =   ascii.read(workdir + "/fit_par_new.txt", Reader=ascii.CommentedHeader)
print(fit_par)
#fit_par.write(workdir + "/fit_par2.txt", format='ascii.basic', overwrite=True, delimiter='\t')

f = open(workdir + "fit_par_new2.txt", 'w')

print("#{: <11} {: <8} {: <8} {: <11} {: <8} {: <11} {: <8} {: <8} {: <8} {: <8} {: <8} {: <8}".format(*fit_par.colnames), file=f)
for row in fit_par:
    if row[2] != '-1':
        print(row[2])
        for i in range(nvisit):
            row[2] = i
            print("{: <12} {: <8} {: <8} {: <11} {: <8} {: <11} {: <8} {: <8} {: <8} {: <8} {: <8} {: <8}".format(*row),
                  file=f)
    else:
        print("{: <12} {: <8} {: <8} {: <11} {: <8} {: <11} {: <8} {: <8} {: <8} {: <8} {: <8} {: <8}".format(*row), file=f)
f.close()

# print("#{: <11} {: <8} {: <8} {: <11} {: <8} {: <11} {: <8} {: <8} {: <8} {: <8} {: <8} {: <8}".format(*fit_par.colnames), file=f)
# for row in fit_par:
#     print("{: <12} {: <8} {: <8} {: <11} {: <8} {: <11} {: <8} {: <8} {: <8} {: <8} {: <8} {: <8}".format(*row), file=f)
# f.close()



fit_par2 = ascii.read(workdir + "fit_par_new2.txt", Reader=ascii.CommentedHeader)
print(fit_par2)

def remove_dupl(seq):
    #https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


print(fit_par2['parameter'])
print(remove_dupl(fit_par2['parameter']))