from pathlib import Path

import numpy as np

from .lib import manageevent as me
from .lib.options import OPTIONS


def run22(eventlabel: str, workdir: Path, meta=None):
    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    if meta.grism == 'G141':
        grism = 'g141_throughput.txt'
    elif meta.grism == 'G102':
        grism = 'g102_throughput.txt'

    Teff, logg, MH = meta.Teff, meta.logg, meta.MH

    for wvl_bins in meta.wvl_bins:
        # What bins do you want?
        # Wave_bins = np.linspace(1.125, 1.65, 22)*1e4
        wave_bins = np.linspace(meta.wvl_min, meta.wvl_max, wvl_bins)*1e4
        print(wave_bins)

        dirname = meta.workdir / "extracted_sp"
        if not dirname.exists():
            dirname.mkdir(parents=True)

        with (meta.workdir / "extracted_sp" /
              f'ld_inputfile_bins{wvl_bins}.txt').open('w', encoding=OPTIONS["encoding"]) as f:
            for i in range(wvl_bins-1):
                params = [i, Teff, logg, MH, '-1', grism, 'P100',
                          wave_bins[i], wave_bins[i+1]]
                params_bin = '\t'.join(map(str, params))
                # print(params_bin)
                print(params_bin, file=f)

    print('Finished s21 \n')
    return meta
