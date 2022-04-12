# ! /usr/bin/env python

import numpy as np
import itertools
"""ramp effect model
2 means two types of traps

original author: Daniel Apai

Version 0.3 fixing trapping parameters

Version 0.2.1 introduce two types of traps, slow traps and fast traps

Version 0.2: add extra keyword parameter to indicate scan or staring
mode observations for staring mode, the detector receive flux in the
same rate during overhead time as that during exposure
precise mathematics forms are included

Version 0.1: Adapted original IDL code to python by Yifan Zhou

"""


def ackBar2(
        cRates,
        tExp,
        exptime=180,
        trap_pop_s=0,
        trap_pop_f=0,
        dTrap_s=0,
        dTrap_f=0,
        dt0=0,
        lost=0,
        mode='scanning'
):
    """Hubble Space Telescope ramp effet model

    Parameters:
    cRates -- intrinsic count rate of each exposures, unit e/s
    tExp -- start time of every exposures
    expTime -- (default 180 seconds) exposure time of the time series
    trap_pop -- (default 0) number of occupied traps at the beginning of the observations
    dTrap -- (default [0])number of extra trap added between two orbits
    dt0 -- (default 0) possible exposures before very beginning, e.g.,
     guiding adjustment
    lost -- (default 0, no lost) proportion of trapped electrons that are not eventually detected
    (mode) -- (default scanning, scanning or staring, or others), for scanning mode
      observation , the pixel no longer receive photons during the overhead
      time, in staring mode, the pixel keps receiving elctrons
    """
    nTrap_s = 1525.38  # 1320.0
    eta_trap_s = 0.013318  # 0.01311
    tau_trap_s = 1.63e4
    nTrap_f = 162.38
    eta_trap_f = 0.008407
    tau_trap_f = 281.463
    try:
        dTrap_f = itertools.cycle(dTrap_f)
        dTrap_s = itertools.cycle(dTrap_s)
        dt0 = itertools.cycle(dt0)
    except TypeError:
        # if dTrap, dt0 provided in scala, convert them to list
        dTrap_f = itertools.cycle([dTrap_f])
        dTrap_s = itertools.cycle([dTrap_s])
        dt0 = itertools.cycle([dt0])
    obsCounts = np.zeros(len(tExp))
    # ensure initial values do not exceed the total trap numbers
    trap_pop_s = min(trap_pop_s, nTrap_s)
    trap_pop_f = min(trap_pop_f, nTrap_f)
    for i in range(len(tExp)):
        try:
            dt = tExp[i+1] - tExp[i]
        except IndexError:
            dt = exptime
        f_i = cRates[i]
        c1_s = eta_trap_s * f_i / nTrap_s + 1 / tau_trap_s  # a key factor
        c1_f = eta_trap_f * f_i / nTrap_f + 1 / tau_trap_f
        # number of trapped electron during one exposure
        dE1_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * (1 - np.exp(-c1_s * exptime))
        dE1_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * (1 - np.exp(-c1_f * exptime))
        dE1_s = min(trap_pop_s + dE1_s, nTrap_s) - trap_pop_s
        dE1_f = min(trap_pop_f + dE1_f, nTrap_f) - trap_pop_f
        trap_pop_s = min(trap_pop_s + dE1_s, nTrap_s)
        trap_pop_f = min(trap_pop_f + dE1_f, nTrap_f)
        obsCounts[i] = f_i * exptime - dE1_s - dE1_f
        if dt < 5 * exptime:  # whether next exposure is in next batch of exposures
            # same orbits
            if mode == 'scanning':
                # scanning mode, no incoming flux between exposures
                dE2_s = - trap_pop_s * (1 - np.exp(-(dt - exptime)/tau_trap_s))
                dE2_f = - trap_pop_f * (1 - np.exp(-(dt - exptime)/tau_trap_f))
            elif mode == 'staring':
                # for staring mode, there is flux between exposures
                dE2_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * (1 - np.exp(-c1_s * (dt - exptime)))
                dE2_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * (1 - np.exp(-c1_f * (dt - exptime)))
            else:
                # others, same as scanning
                dE2_s = - trap_pop_s * (1 - np.exp(-(dt - exptime)/tau_trap_s))
                dE2_f = - trap_pop_f * (1 - np.exp(-(dt - exptime)/tau_trap_f))
            trap_pop_s = min(trap_pop_s + dE2_s, nTrap_s)
            trap_pop_f = min(trap_pop_f + dE2_f, nTrap_f)
        elif dt < 1200:
            # considering in-orbit buffer download scenario
            if mode == 'staring':
                trap_pop_s = min(trap_pop_s * np.exp(-(dt-exptime)/tau_trap_s), nTrap_s)
                trap_pop_f = min(trap_pop_f * np.exp(-(dt-exptime)/tau_trap_f), nTrap_f)
        else:
            # switch orbit
            dt0_i = next(dt0)
            trap_pop_s = min(trap_pop_s * np.exp(-(dt-exptime-dt0_i)/tau_trap_s) + next(dTrap_s), nTrap_s)
            trap_pop_f = min(trap_pop_f * np.exp(-(dt-exptime-dt0_i)/tau_trap_f) + next(dTrap_f), nTrap_f)
            f_i = cRates[i + 1]
            c1_s = eta_trap_s * f_i / nTrap_s + 1 / tau_trap_s  # a key factor
            c1_f = eta_trap_f * f_i / nTrap_f + 1 / tau_trap_f
            dE3_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * (1 - np.exp(-c1_s * dt0_i))
            dE3_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * (1 - np.exp(-c1_f * dt0_i))
            dE3_s = min(trap_pop_s + dE3_s, nTrap_s) - trap_pop_s
            dE3_f = min(trap_pop_f + dE3_f, nTrap_f) - trap_pop_f
            trap_pop_s = min(trap_pop_s + dE3_s, nTrap_s)
            trap_pop_f = min(trap_pop_f + dE3_f, nTrap_f)
        trap_pop_s = max(trap_pop_s, 0)
        trap_pop_f = max(trap_pop_f, 0)

    return obsCounts


if __name__ == '__main__':
    pass
    # import matplotlib.pyplot as plt
    # t1 = np.linspace(0, 2700, 80)
    # t2 = np.linspace(5558, 8280, 80)
    # t = np.concatenate((t1, t2))
    # crate = 100
    # crates = crate * np.ones(len(t))
    # dataDIR = '/Users/ZhouYf/Documents/HST14241/alldata/2M0335/DATA/'
    # from os import path
    # import pandas as pd

    # info = pd.read_csv(
    #     path.expanduser('~/Documents/HST14241/alldata/2M0335/2M0335_fileInfo.csv'),
    #     parse_dates=True,
    #     index_col='Datetime')
    # info['Time'] = np.float32(info.index - info.index.values[0]) / 1e9
    # expTime = info['Exp Time'].values[0]
    # grismInfo = info[info['Filter'] == 'G141']
    # tExp = grismInfo['Time'].values
    # # cRates = np.ones(len(LC)) * LC.mean() * 1.002
    # cRates = np.ones(len(tExp)) * 100
    # obs = ackBar2(cRates, tExp, exptime=expTime, lost=0,
    #               mode='scanning')
    # plt.close('all')
    # # plt.plot(tExp, LC*expTime, 'o')
    # plt.plot(tExp, obs, '-')
    # # plt.ylim([crate * 0.95, crate * 1.02])
    # plt.show()
