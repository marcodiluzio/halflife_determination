"""Script to visualize results of half-life_determination

Summary:
...

it contains functions:
????, ????

it contains classes:
Plotbirks

This module can be imported into another script with:
"import half-life_determination"                      #whole package
"from half-life_determination import visualization"   #single module
giving access to all its methods and classes

author:  Marco Di Luzio
email:   m.diluzio@inrim.it
date:    2026-01-22
version: 0.0.1
"""

#imports
import os
import datetime
#from itertools import product
import configparser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Plotbirks:
    """Class for visualization of the birks fits
    
    Attributes
    ----------
    self.fig : matpotlib.figure.Figure
        figure containing all the axes
    self.ax1 : matpotlib.axes.Axes
        axes diplaying the projection of the intercept of fits
    self.ax2 : matpotlib.axes.Axes
        axes diplaying normalized activity over TDCR ratios and corresponding linear fits
    self.ax4 : matpotlib.axes.Axes
        axes diplaying the residuals of the fits
    self.pkwargs : dict
        dictionary with datapoints plotting settings
    self.lkwargs : dict
        dictionary with lines plotting settings
    self.axmin : list
        list of the lowest abscissa for each series of data
    
    Methods
    -------
    add_dataset()
        add data to the plot
    show()
        display the figure
    """
    
    def __init__(self, figsize=(7,5), title=''):
        """Plotbirks class constructor
        
        Parameters
        ----------
        figsize : tuple
            tuple with width and height dimensions of plot in inches (default (7,5))
        title : str
            set the title of the entire figure (default '')
            
        Return
        ------
        None        
        """
        self.fig, ((self.ax1, self.ax2), (ax3, self.ax4)) = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, squeeze=True, width_ratios=(1,3), height_ratios=(3,1), figsize=figsize)
        self.pkwargs = {'marker':'o', 'markersize':3, 'linestyle':'', 'elinewidth':0.3, 'ecolor':'k'}
        self.lkwargs = {'marker':'', 'linestyle':'-', 'linewidth':0.75, 'color':'k'}
        self.title = title
        self.axmin = []
        
        ax3.remove()
        
    def add_dataset(self, x, y, uy, res, unc, residuals, label):
        """add dataset to show on self.axes
        
        Parameters
        ----------
        x : np.array
            x data
        y : np.array
            y data
        uy : res
            standard uncertainties of y data
        res : float
            intercept of fit
        unc : float
            uncertianties of res
        residuals : np.array
            residuals of fit on y
        label : str
            descr

        Return
        ------
        None
        """
        
        self.ax1.errorbar(0, res, yerr=2*unc, **self.pkwargs, label=label)
        self.ax1.plot([0,np.max(x)], [res, 1], **self.lkwargs)
        
        self.ax2.errorbar(x, y, yerr=2*uy, **self.pkwargs, label=label)
        self.ax2.plot([0,np.max(x)], [res, 1], **self.lkwargs)
        
        #residuals
        self.ax4.errorbar(x, residuals, yerr=2*uy, **self.pkwargs)

        self.axmin.append(np.min(x))
        
    def show(self):
        """adjust the axes, set the title and show the figure
        
        Parameters
        ----------
        None
        
        Return
        ------
        None
        """
        self.ax1.set_xlim(-0.1, 0.1)
        self.ax1.set_xticks([0.0], ['intercept'])
        self.ax1.set_ylabel(r'$A_\mathrm{norm}$ / 1')
        self.ax2.set_xlim(np.min(self.axmin)-0.01,1)
        self.ax2.set_xticklabels([])
        self.ax4.set_xlim(np.min(self.axmin)-0.01,1)
        self.ax4.set_xlabel('TDCR / 1')
        self.ax4.set_ylabel('residuals')
        self.ax4.axhline(y=0, xmin=0, xmax=1, **self.lkwargs)
        self.fig.suptitle(self.title, fontsize=12)
        
        self.ax2.legend(loc='best', ncols=2)
        
        self.fig.tight_layout()
        plt.show()
        

def _plot(filename):
    with pd.read_csv(filename) as _df:
        print(_df)