"""Script to visualize results of half-life_determination

Summary:
...

it contains functions:
_fit_plot, _fit_plot_M, _distribution_plot, plot_results, _recall_dataframe

it contains classes:
Plotbirks

This module can be imported into another script with:
"import halflife_determination"                      #whole package
"from halflife_determination import visualization"   #single module
giving access to all its methods and classes

author:  Marco Di Luzio
email:   m.diluzio@inrim.it
date:    2026-01-22
"""

#imports
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
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


def _fit_plot(x, y, fit_x, fit_y, residuals, uy=None, k=2, suptitle='', fitted_HL='0', return_plot=False):
    """Display fit (fit_x, fit_y) on the dataset (x, y)
    
    Parameters
    ----------
    x : numpy.array
        normalized time series
    y : numpy.array
        normalized activity
    fit_x : numpy.array
        values of x to fit
    fit_y : numpy.array
        fitted values of y
    residuals : numpy.array
        residuals of the fit
    uy : numpy.array
        normalized uncertainty of y (default None)
    k : float
        coverage factor for expanded uncertainty (default 2)
    suptitle : str
        title of the plot (default '')
    fitted_HL : str
        the fitted half-life value (default '0')
    return_plot : bool
        whether to return a tuple of fig and axes objects or None (default False)
    
    Return
    ------
    fig : plt.Figure    [if return_plot=True]
        figure containing Axes
    (ax1, ax2) : tuple  [if return_plot=True]
        tuple of figure Axes
    None                [if return_plot=False]
    """
    style_points = {'linestyle':'', 'marker':'o', 'markersize':2.5, 'markerfacecolor':'k', 'color':'k'}
    style_line = {'linestyle':'-', 'marker':'', 'color':'r'}

    if uy is None:
        uy = x * 0

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7,5), height_ratios=(2,1))
    ax1.errorbar(x, y, yerr=k*uy, **style_points, elinewidth=0.75)
    ax1.plot(fit_x, fit_y, **style_line)
    ax1.set_ylabel(r'ln$\left(A_\mathrm{n}\right)$ / 1')
    
    ax2.errorbar(x, residuals, yerr=k*uy, **style_points, elinewidth=0.75)
    ax2.axhline(y=0, xmin=0, xmax=1, **style_line)
    ax2.set_xlabel(r'$t_\mathrm{n}$ / $d$')
    ax2.set_ylabel('residuals / 1')
    
    ax1.text(0.75, 0.9, rf'$t_{{{"1/2"}}}$ = {fitted_HL} d', transform=ax1.transAxes)
    
    fig.suptitle(suptitle, fontsize=12)
    
    fig.tight_layout()
    if return_plot:
        return fig, (ax1, ax2)
    plt.show()
    
def _fit_plot_M(x, y, fit_x, fit_y, residuals, masque, uy=None, k=2, suptitle='', fitted_HL='0', return_plot=False):
    """Display fit (fit_x, fit_y) on the masqued dataset (x[masque], y[masque])
    
    Parameters
    ----------
    x : numpy.array
        normalized time series
    y : numpy.array
        normalized activity
    fit_x : numpy.array
        values of x to fit
    fit_y : numpy.array
        fitted values of y
    residuals : numpy.array
        residuals of the fit
    masque : list
        list of bool of data adopted for the elaboration
    uy : numpy.array
        relative uncertainty of y (default None)
    k : float
        coverage factor for expanded uncertainty (default 2)
    suptitle : str
        title of the plot (default '')
    fitted_HL : str
        the fitted half-life value (default '0')
    return_plot : bool
        whether to return a tuple of fig and axes objects or None (default False)
    
    Return
    ------
    fig : plt.Figure    [if return_plot=True]
        figure containing Axes
    (ax1, ax2) : tuple  [if return_plot=True]
        tuple of figure Axes
    None                [if return_plot=False]
    """
    style_points = {'linestyle':'', 'marker':'o', 'markersize':2.5, 'markerfacecolor':'k', 'color':'k', 'elinewidth':0.75}
    style_crosses = {'linestyle':'', 'marker':'x', 'markersize':3, 'markerfacecolor':'k', 'color':'k', 'elinewidth':0.75}
    style_line = {'linestyle':'-', 'marker':'', 'color':'r'}
    masque = np.array(masque)

    if uy is None:
        uy = x * 0
    
    MX = np.array([xi for xi,mi in zip(x,masque) if mi == True])
    MY = np.array([yi for yi,mi in zip(y,masque) if mi == True])
    MUY = np.array([uyi for uyi,mi in zip(uy,masque) if mi == True])
    
    XX = np.array([xi for xi,mi in zip(x,~masque) if mi == True])
    XY = np.array([yi for yi,mi in zip(y,~masque) if mi == True])
    XUY = np.array([uyi for uyi,mi in zip(uy,~masque) if mi == True])

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7,5), height_ratios=(2,1))
    ax1.errorbar(MX, MY, yerr=k*MUY, **style_points)
    ax1.errorbar(XX, XY, yerr=k*XUY, **style_crosses)
    ax1.plot(fit_x, fit_y, **style_line)
    ax1.set_ylabel(r'ln$\left(A_\mathrm{n}\right)$ / 1')
    
    ax2.errorbar(MX, residuals, yerr=k*MUY, **style_points)
    ax2.axhline(y=0, xmin=0, xmax=1, **style_line)
    ax2.set_xlabel(r'$t_\mathrm{n}$ / $d$')
    ax2.set_ylabel('residuals / 1')
    
    ax1.text(0.75, 0.9, rf'$t_{{{"1/2"}}}$ = {fitted_HL} d', transform=ax1.transAxes)
    
    fig.suptitle(suptitle, fontsize=12)
    
    fig.tight_layout()
    if return_plot:
        return fig, (ax1, ax2)
    plt.show()

def _distribution_plot(Y, suptitle='', fitted_HL='0', bins=20, return_plot=False):
    """Display distribution of estimated result of a MonteCarlo fit
    
    Parameters
    ----------
    Y : numpy.array
        array of results of the fit
    suptitle : str
        title of the plot (default '')
    fitted_HL : str
        the fitted half-life value (default '0')
    bins : int
        number of bins of the histogram visualization (default 20)
    return_plot : bool
        whether to return a tuple of fig and axes objects or None (default False)
    
    Return
    ------
    fig : plt.Figure    [if return_plot=True]
        figure containing Axes
    (ax1, ax2) : tuple  [if return_plot=True]
        tuple of figure Axes
    None                [if return_plot=False]
    """
    style_points = {'linestyle':'', 'marker':'o', 'markersize':2.5, 'markerfacecolor':'k', 'color':'k'}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7,5), width_ratios=(2.5,1), sharey=True)
    
    ax1.plot(Y, **style_points)
    ax1.text(0.75, 0.95, rf'$t_{{{"1/2"}}}$ = {fitted_HL} d', transform=ax1.transAxes)
    
    ax2.hist(Y, bins=bins, orientation='horizontal', color='r')
    ax2.set_xticks([])
    ax2.set_xlim(0, None)
    
    fig.suptitle(suptitle, fontsize=12)
    
    fig.tight_layout()
    if return_plot:
        return fig, (ax1, ax2)
    plt.show()

def plot_results(fitted_data, value_label, uncertainty_label, averages=None, k=2, bins=None, title='', return_plot=False):
    """Diplay the results for any independent series of measurements
    
    Parameters
    ----------
    fitted_data : pandas.DataFrame
        dataframe containing the fitted data
    value_label : str
        label of the dataframe to recall the values column
    uncertainty_label : str
        label of the dataframe to recall the uncertainties column
    averages : dict
        dict containing tuples with results of averaging procedure performed on these data (default None)
    k : float
        coverage factor for expanded uncertainty (default 2)
    bins : int
        number of bins for the histogram representation of results (default None)
    title : str
        title displaying on top of the plot (default '')
    return_plot : bool
        whether to return a tuple of fig and axes objects or None (default False)
    
    Return
    ------
    fig : plt.Figure    [if return_plot=True]
        figure containing Axes
    (ax1, ax2) : tuple  [if return_plot=True]
        tuple of figure Axes
    None                [if return_plot=False]
    """
    style_points = {'linestyle':'', 'marker':'o', 'markersize':2.5, 'markerfacecolor':'k', 'color':'k', 'elinewidth':0.75}
    
    sorted_data = fitted_data.sort_values(value_label, axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None)

    xlabels = [f'EDT: {_A}, COI: {_B}, CT: {_C}, ID: {_D}' for _A, _B, _C, _D in zip(sorted_data['ext_dead_time'], sorted_data['coincidence_w'], sorted_data['cocktail'], sorted_data['source'])]

    x_plot = np.arange(len(sorted_data[value_label]))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,7), width_ratios=(2.5,1), sharey=True)
    
    ax1.errorbar(x_plot, sorted_data[value_label], yerr=k*sorted_data[uncertainty_label]**0.5, **style_points)
    ax1.set_xticks(x_plot, xlabels, rotation=90)
    ax1.set_ylabel(r'$t_{1/2}$ / $d$')
    
    colors = ['r', '#00008B', '#FFEB2A', '#AA00FF', '#66FF00', '#1E2136', '#62866C', '#81B2D9', '#BBA6DD', '#64557B', '#1E2136']
    if averages is not None:
        for (_key, _value) in averages.items():
            _color = colors.pop(0)
            ax1.axhline(y=_value[0], xmin=0, xmax=1, linestyle='-', color=_color, linewidth=1.75, label=_key)
            ax1.axhline(y=_value[0]-k*_value[1], xmin=0, xmax=1, linestyle='--', color=_color, linewidth=1.0)
            ax1.axhline(y=_value[0]+k*_value[1], xmin=0, xmax=1, linestyle='--', color=_color, linewidth=1.0)    
    if bins is None:
        bins = int(len(sorted_data[value_label]) / 5)
    
    ax2.hist(sorted_data[value_label], bins=bins, orientation='horizontal', color='#00008B')
    ax2.set_xticks([])
    ax2.set_xlim(0, None)
    
    fig.suptitle(title)
    fig.legend(loc='lower right')

    fig.tight_layout()
    if return_plot:
        return fig, (ax1, ax2)
    plt.show()

def _recall_dataframe(filename, multiindex=True):
    """Diplay the results for any independent series of measurements
    
    Parameters
    ----------
    filename : str (or Path)
        filename for data to be returned as DataFrame
    multiindex : bool
        switches the index with a multiindex containing multiple columns (default True)
    
    Return
    ------
    _df : pandas.DataFrame
        corresponding data as a pandas.DataFrame
    """
    
    original_data = {'nuclide', 'source', 'meas_date', 'meas_time',
    'meas_dead_time', 'device', 'ext_dead_time', 'coincidence_w',
    'ndfilter', 'cocktail', 'volume', 'rates_DoublesRaw',
    'rates_DoublesRaw_u', 'rates_TriplesRaw', 'rates_TriplesRaw_u',
    'rates_Doubles', 'rates_Doubles_u', 'rates_Triples', 'rates_Triples_u',
    'TDCR', 'TDCR_u', 'MSefficiencyD', 'MSefficiencyD_u', 'C_Activity',
    'C_Activity_u', 'birks', 'filename'}
    
    elaborated_data = {'nuclide', 'meas_date', 'ext_dead_time', 'coincidence_w',
    'ndfilter', 'cocktail', 'volume', 'source', 'Act', 'Act_u', 'birks',
    'filename', 'norm_Act', 'norm_date'}
    
    fitted_constant = {'ext_dead_time', 'coincidence_w', 'cocktail', 'birks', 'source'}
    fitted_variable = ({'wl_HL', 'wl_uHL', 'wl_uHF', 'wl_uMF', 'wl_uLF', 'wl_combined_variance'},
    {'l_HL', 'l_uHL', 'l_uHF', 'l_uMF', 'l_uLF', 'l_combined_variance'},
    {'wexp_HL', 'wexp_uHL', 'wexp_uHF', 'wexp_uMF', 'wexp_uLF', 'wexp_combined_variance'},
    {'exp_HL', 'exp_uHL', 'exp_uHF', 'exp_uMF', 'exp_uLF', 'exp_combined_variance'},
    {'mcl_HL', 'mcl_uHL', 'mcl_uHF', 'mcl_uMF', 'mcl_uLF', 'mcl_combined_variance'},
    {'mcexp_HL', 'mcexp_uHL', 'mcexp_uHF', 'mcexp_uMF', 'mcexp_uLF', 'mcexp_combined_variance'})
    
    results_variable = ({'montecarlo exponential fit_value / d', 'montecarlo exponential fit_uncertainty / d'},
    {'exponential fit_value / d', 'exponential fit_uncertainty / d'},
    {'weighted exponential fit_value / d', 'weighted exponential fit_uncertainty / d'},
    {'linear fit_value / d', 'linear fit_uncertainty / d'},
    {'weigthed linear fit_value / d', 'weigthed linear fit_uncertainty / d'},
    {'montecarlo linear fit_value / d', 'montecarlo linear fit_uncertainty / d'})
    
    _df = pd.read_csv(filename)
    
    fitted_check = [subset.issubset(_df.columns) for subset in fitted_variable]
    result_check = [subset.issubset(_df.columns) for subset in results_variable]

    if original_data.issubset(_df.columns):
        #original_data
        if multiindex:
            _df.set_index(['nuclide', 'device', 'ext_dead_time', 'coincidence_w', 'cocktail', 'source', 'birks'], inplace=True)
        return _df
    
    elif elaborated_data.issubset(_df.columns):
        #elaborated_data
        if multiindex:
            _df.set_index(['nuclide', 'ext_dead_time', 'coincidence_w', 'cocktail', 'source', 'birks'], inplace=True)
        return _df
        
    elif fitted_constant.issubset(_df.columns) and np.sum(fitted_check) > 0:
        if multiindex:
            _df.set_index(['ext_dead_time', 'coincidence_w', 'cocktail', 'source', 'birks'], inplace=True)
        return _df
        
    elif np.sum(result_check) > 0:
        _df.set_index('Unnamed: 0', inplace=True)
        return _df

    return _df