from . import read_pcf as rd
from . import manageevent as me


def update_meta(eventlabel, workdir):

    meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    pcffile = workdir + 'obs_par.pcf'
    pcf = rd.read_pcf(pcffile)
    rd.store_pcf(meta, pcf)

    meta.workdir = workdir
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Successfully reloaded meta file')

    return 0
