from . import readECF as rd
from . import manageevent as me


def reload_meta(workdir, eventlabel):

    meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    ecffile = workdir + 'obs_par.ecf'
    ecf = rd.read_ecf(ecffile)
    rd.store_ecf(meta, ecf)

    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Succsesfully reloaded meta file')

    return 0
