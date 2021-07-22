import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import seaborn as sns
from .formatter import FormatParams
from .model import Model, calc_sys, calc_astro

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})
matplotlib.rcParams.update({'lines.markeredgewidth':0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset':False})

def plot_raw(data):
    palette = sns.color_palette("husl", data.nvisit)
    for i in range(data.nvisit):
        ind = data.vis_num==i
        plt.subplot((data.nvisit)*100+10+i+1)
        plt.plot(data.t_vis[ind]*24., data.flux[ind], marker='o', \
                markersize=3.0, linestyle="none", \
                label = "Visit {0}".format(i), color= palette[i])
        plt.xlim(((data.t_vis.min()-0.02)*24., (data.t_vis.max()+0.05)*24.))
        plt.ylim((0.998*data.flux.min(), 1.002*data.flux.max()))
        plt.legend()
    plt.xlabel("Time after visit start (hours)")
    plt.ylabel("Flux (e-)")
    plt.tight_layout()
    plt.savefig("raw_lc.pdf")
    plt.show()


def plot_fit(data, model):
    p = FormatParams(model.params, data)    #FIXME 

    sns.set_palette("muted")
    palette = sns.color_palette("muted", data.nvisit)

    #ind = model.phase > 0.5
    #model.phase[ind] -= 1.

    #calculate a range of times at higher resolution to make model look nice
    phase_hr = np.linspace(model.phase.min()-0.05, model.phase.max()+0.05, 1000)
    t_hr = phase_hr*p.per[0]+p.t0[0] + data.toffset

    #colors = ['blue', 'red', 'orange', 'gray']

    #plot data
    plt.subplot(211)
    #plot best fit model from first visit
    plt.plot(phase_hr, calc_astro(t_hr, model.params, data, model.myfuncs, 0))

    #plot systematics removed data
    for i in range(data.nvisit):
        ind = data.vis_num == i
        plt.plot(model.phase[ind], model.data_nosys[ind], color = palette[i], marker = 'o', markersize = 3, linestyle = "none")

    #add labels/set axes
    #xlo, xhi = np.min(model.phase)*0.9, np.max(model.phase)*1.1
    xlo, xhi = -0.1, 0.1
    plt.xlim(xlo,xhi)
    plt.ylabel("Relative Flux")

    #annotate plot with fit diagnostics
    ax = plt.gca()
    ax.text( 0.85, 0.29, 
             '$\chi^2_\\nu$:    ' + '{0:0.2f}'.format(model.chi2red) + '\n'
            + 'obs. rms:  ' + '{0:0d}'.format(int(model.rms)) + '\n' 
            + 'exp. rms:  ' + '{0:0d}'.format(int(model.rms_predicted)), 
            verticalalignment='top',horizontalalignment='left', 
            transform=ax.transAxes, fontsize = 12
    )
    
    #plot fit residuals
    plt.subplot(212)
    plt.axhline(0, zorder=1, color='0.2', linestyle='dashed')

    for i in range(data.nvisit):
        ind = data.vis_num == i
        plt.plot(model.phase[ind], 1.0e6*model.norm_resid[ind], color = palette[i], marker = 'o', markersize = 3, linestyle = "none")

    #add labels/set axes
    plt.xlim(xlo,xhi)
    plt.ylabel("Residuals (ppm)")
    plt.xlabel("Orbital phase")
    
    plt.tight_layout()
    plt.savefig("white_lc.pdf")
    plt.show()
    #plt.waitforbuttonpress(0) # this will wait for indefinite time
    plt.close()


