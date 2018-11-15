from __future__ import print_function
import sys, os, shutil
import matplotlib.pyplot as plt

#------------------------------------------------------------------------
class SimFemPlot(object):
    """Plot.
    """

#------------------------------------------------------------------------
    def setNames(self):
        if self.methods:
            self.methodsnames = {}
            for method in self.methods:
                self.methodsnames[method] = method
        self.paramname = self.param
#------------------------------------------------------------------------
    def __init__(self, methods=None, params=None, param='hmean'):
        self.methods = methods
        self.params = params
        self.param = param
        self.filename = None
        self.order = True
        self.setNames()
#---------------------------------------------------------#
    def ploterrors(self, datatoplot, keys=None, scale='loglog'):
        nrows = 1
        ncols = len(datatoplot)
        fig = plt.figure(figsize=(5*ncols,5*nrows))
        markers=['o', 'x', '+', '*', 'd', 's', '<', '>']
        colors=['r', 'g', 'b', 'y', 'm', 'c', 'k', 'gray']
        iplot = 1
        if not keys: keys = datatoplot.keys()
        for key in keys:
            data = datatoplot[key]
            ax = fig.add_subplot(nrows, ncols, iplot)
            iplot = iplot + 1
            for im,method in enumerate(self.methods):
                if scale =='loglog':
                    ax.loglog(self.params, data[method], marker=markers[im], color=colors[im], label=self.methodsnames[method], lw=2)
                else:
                    ax.plot(self.params, data[method], marker=markers[im], color=colors[im], label=self.methodsnames[method], lw=2)
            if self.param == "hmean" and self.order:
                    if key in ["H1", "E"]:
                        ax.loglog(self.params, self.params, linestyle="--", color='k', label='order 1')
                    if key in ["L1", "L2", "Linf"]:
                        paramsq = [p**2 for p in self.params]
                        ax.loglog(self.params, paramsq, linestyle="--", color='k', label='order 2')
            keyname = key
            if key == "L1": keyname = r"$L^1(\Omega)$"
            elif key == "L2": keyname = r"$L^2(\Omega)$"
            elif key == "H1": keyname = r"$H^1(\Omega)$"
            elif key == "Linf": keyname = r"$L^{\infty}(\Omega)$"
            ax.set_title("%s" %(keyname))
            ax.legend(loc='best')
            ax.grid(True)
            plt.xlabel(self.paramname)
            # plt.ylabel('e(%3s)' %(key))
        if self.filename is not None: plt.savefig(self.filename)
        plt.show()
