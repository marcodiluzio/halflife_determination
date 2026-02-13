"""Script to perform half-life measurement on TDCR data

Summary:
half-life_determination is a module for half-life evaluation of radionuclides which activity is monitored collecting information on beta decay through a series of TDCR (Triple-to-Douple Coincidence Ratio) liquid scintillation measurements.
It is suitable for long-lived radionuclides having half-lives in the order of weeks to years.
Uncertainty evaluation is performed by adopting GUM procedures.

it contains functions:
open_result_file, _get_time_value, _get_activity, _get_category_value, _linear_fitting_procedure_birks, _get_filenames, _get_birks_best_value, get_HL_data_from_dir, renormalize_data, _exp, _exponential_fitting_procedure, _montecarlo_fitting_procedure, _linear_fitting_procedure, _linear_fitting_procedure_M, fit_data, _get_autocorrelation, PMM_method, BirgeAdjust, DerSimonianLairdp, CoxProcedureA, CoxProcedureB, iterative_procedure, get_result, elaboration, read_info, load_config

This module can be imported into another script with:
"from halflife_determination import hl_elaboration"   #single module
giving access to all their corrsponding methods and classes

Or can be directly used in the command line with:
"python -m halflife_determination {argument_1} {optional argument_2}"
in this case it takes 2 arguments
{argument_1} = name or path of the folder to search
{optional argument_2} = name of the configuration file

author:  Marco Di Luzio
email:   m.diluzio@inrim.it
"""

#imports
import os
import datetime
from itertools import product
import configparser
import pickle
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2
import consensusgen as csg
try:
    from halflife_determination import visualization
except ImportError:
    import visualization


def _get_time_value(line):
    """Return a time value (converted in s) from the corresponding dataframe line
    
    Parameters
    ----------
    line : dataframe.row
        row of a dataframe containing labels: 'Value' 'Unit'
    
    Return
    ------
    value : float
        time value in seconds
    """
    
    conversion = {'ns':1E-9, 'us':1E-6, 'ms':1E-3, 'µs':1E-6, '%':1E-2}

    value = float(line.Value)    
    return value * conversion.get(line.Unit, 1)
    
def _get_category_value(line, unit=True):
    """Return a str from the corresponding dataframe line
    
    Parameters
    ----------
    line : dataframe.row
        row of a dataframe containing labels: 'Value' 'Unit'
    unit : bool
        whether the unit has to be appended to the value (default True)
    
    Return
    ------
    category : str
        str combination of Value (and Unit)
    """
    
    if unit:
        return f'{line.Value} {line.Unit}'
    return  f'{line.Value}'

def _get_activity(df, key, secondary_key='Activity'):
    """Return an activity value and uncertainty (converted in Bq) from the corresponding dataframe line (could be also used to retrieve efficiency data)
    
    Parameters
    ----------
    df : pandas.DataFrame
        dataframe containing multiindex and labels: 'Value' 'Unit'
    key : str
        1st key of the multi index of the dataframe
    secondary_key : str
        2nd key (default 'Activity')
    
    Return
    ------
    value : float
        activity value in Bq (or np.nan)
    unc : float
        uncertainty of activity value in Bq (or np.nan)
    """
    
    conversion = {'kBq':1E3}
    
    try:
        line = df.loc[(key, secondary_key)]
        value = float(line.Value)
        unc = float(line.Uncertainty)
    except (KeyError, ValueError, TypeError):
        return np.nan, np.nan
    return value * conversion.get(line.Unit, 1), unc * conversion.get(line.Unit, 1)
    
def _linear_fitting_procedure_birks(X, Y, UY=None):
    """Perform linear fit on the dataset (x, y) to obtain estimate for intercept of the straight line
    
    Parameters
    ----------
    X : numpy.array
        array containing TDCR series
    Y : numpy.array
        array containing the normalized activity at different TDCR values
    UY : numpy.array
        array containing the normalized activity uncertainties at different TDCR values (default None)
    
    Return
    ------
    relevant_parameter : float
        intercept of the fitting straight line
    relevant_parameter_uncertainty : float
        uncertainty of the intercept of the fitting straight line
    """
    #creation of independent variables matrix
    W = X[:, np.newaxis]**[0,1]
    
    if UY is not None:
        V = np.linalg.inv(np.diag(np.power(UY, 2)))#weighted on variances
    else:
        V = np.identity(W.shape[0])

    #parameters estimation
    parameters = np.linalg.inv(W.T@V@W)@(W.T@V@Y)
    
    #residuals
    residuals = Y - parameters@W.T
    
    #covariance matrix
    n, k = Y.shape[0], W.shape[1]
    mcov = np.true_divide(1, n-k) * (residuals.T@V@residuals) * np.linalg.inv(W.T@V@W)
    
    #intercept (index 0) is taken as representative result
    
    return parameters[0], np.sqrt(np.diag(mcov))[0]

def open_result_file(filename):
    """Retrieve relevant half-life information from a csv result file using pandas module
    
    Parameters
    ----------
    filename : str (Path)
        file path of the file from which the information are recalled
    
    Return
    ------
    row : pandas.Series (or None)
        Series containing all the relevant information of the open csv file
    """
    # in case filename is not absolute
    filename = os.path.abspath(filename)

    try:
        df = pd.read_csv(filename, sep=';', header=0, index_col=(0,1), dtype=object)
    except (PermissionError, IndexError, pd.errors.EmptyDataError):
        return None
    
    try:
        #nuclide
        nuclide = df.loc[('General', 'Nuclide')]
        nuclide = nuclide.Value
    except KeyError:
        return None
        
    try:
        #source lab (actually the sample ID from the filename)
        manipualte_name = os.path.basename(filename)
        source_split = manipualte_name.split('_')
        if len(source_split) == 12:
            source_string = manipualte_name.split('_')[6]
        else:
            source_string = manipualte_name.split('_')[5]
            
        #default        
        sourcelab = '9999'
        for _nn, char in enumerate(source_string):
            if char.isnumeric():
                sourcelab = source_string[_nn:]
    except KeyError:
        return None
    
    try:
        #measurement date
        meas_date = df.loc[('General', 'Measurement Date')]
        #datetime string: 2025-06-30 09:05:00
        format_string = '%Y-%m-%d %H:%M:%S'
        meas_date = datetime.datetime.strptime(meas_date.Value, format_string)
    except (KeyError, ValueError):
        return None
    
    try:
        #measurement time
        meas_time = df.loc[('DAQ', 'Measurement Time')]
        meas_time = float(meas_time.Value)
        #print(meas_time.Uncertainty)
    except (KeyError, ValueError):
        return None
        
    try:
        #device
        device = df.loc[('DAQ', 'Device')]
        device = _get_category_value(device, unit=False)
    except KeyError:
        return None
    
    try:
        #dead time
        dead_time = df.loc[('DAQ', 'Dead Time (Extended)')]
        dead_time = _get_category_value(dead_time)
    except KeyError:
        return None
        
    try:
        #coincidence window
        coincidence_w = df.loc[('DAQ', 'Coincidence Window')]
        coincidence_w = _get_category_value(coincidence_w)
    except KeyError:
        return None

    try:
        #cocktail
        cocktail = df.loc[('Source', 'LS Cocktail')]
        cocktail = cocktail.Value
    except KeyError:
        return None
        
    try:
        #volume
        volume = df.loc[('Source', 'Volume')]
        volumeF = float(volume.Value)
        volume = volume.Value
    except (KeyError, TypeError):
        return None        
    
    try:
        #filter
        ndfilter = df.loc[('Source', 'Grey Filter ND')]
        ndfilter = ndfilter.Value
    except KeyError:
        return None
    
    #only double and triple coincidence should matter for us
    #rates (all in s-1)
    #rates_A = df.loc[('Rates', 'Singles Raw (A)')]
    #rates_A_u = float(rates_A.Uncertainty)
    #rates_A = float(rates_A.Value)
    
    #rates_B = df.loc[('Rates', 'Singles Raw (B)')]
    #rates_B_u = float(rates_B.Uncertainty)
    #rates_B = float(rates_B.Value)
    
    #rates_C = df.loc[('Rates', 'Singles Raw (C)')]
    #rates_C_u = float(rates_C.Uncertainty)
    #rates_C = float(rates_C.Value)
    
    #doubles
    #rates_AB = df.loc[('Rates', 'Doubles Raw (AB)')]
    #rates_AB_u = float(rates_AB.Uncertainty)
    #rates_AB = float(rates_AB.Value)
    
    #rates_BC = df.loc[('Rates', 'Doubles Raw (BC)')]
    #rates_BC_u = float(rates_BC.Uncertainty)
    #rates_BC = float(rates_BC.Value)
    
    #rates_AC = df.loc[('Rates', 'Doubles Raw (AC)')]
    #rates_AC_u = float(rates_AC.Uncertainty)
    #rates_AC = float(rates_AC.Value)
    
    try:
        rates_DoublesRaw = df.loc[('Rates', 'Doubles Raw (D)')]
        rates_DoublesRaw_u = float(rates_DoublesRaw.Uncertainty)
        rates_DoublesRaw = float(rates_DoublesRaw.Value)
    except (KeyError, TypeError):
        rates_DoublesRaw, rates_DoublesRaw_u = np.nan, np.nan
    try:
        rates_TriplesRaw = df.loc[('Rates', 'Triples Raw (T)')]
        rates_TriplesRaw_u = float(rates_TriplesRaw.Uncertainty)
        rates_TriplesRaw = float(rates_TriplesRaw.Value)
    except (KeyError, TypeError):
        rates_TriplesRaw, rates_TriplesRaw_u = np.nan , np.nan
    
    try:
        #these are the important ones!
        rates_Doubles = df.loc[('Rates', 'Doubles Net (D-D0)')]
        rates_Doubles_u = float(rates_Doubles.Uncertainty)
        rates_Doubles = float(rates_Doubles.Value)
    except (KeyError, ValueError):
        return None
    
    try:
        rates_Triples = df.loc[('Rates', 'Triples Net (T-T0)')]
        rates_Triples_u = float(rates_Triples.Uncertainty)
        rates_Triples = float(rates_Triples.Value)
    except (KeyError, ValueError):
        return None
        
    try:
        #TDCR
        TDCR = df.loc[('Rates', 'TDCR')]
        TDCR_u = float(TDCR.Uncertainty)
        TDCR = float(TDCR.Value)
    except (KeyError, ValueError):
        return None
    
    #corrections
    try:
        acc_correction = df.loc[('Rates', 'Accidental Corr.')]
        acc_correction = _get_time_value(acc_correction)
    except KeyError:
        acc_correction = 0.0
    
    try:
        meas_dead_time = df.loc[('Rates', 'Measurement Dead Time')]
        meas_dead_time = _get_time_value(meas_dead_time)
    except (KeyError, TypeError):
        return None
    
    try:
        #single uncertainties
        u_rel_counting_statistics = df.loc[('Uncertainty Budget', 'u_rel Counting Stat')]
        u_rel_counting_statistics = _get_time_value(u_rel_counting_statistics)
    except KeyError:
        u_rel_counting_statistics = 0.0
    
    try:
        u_rel_decay_correction = df.loc[('Uncertainty Budget', 'u_rel Decay Corr')]
        u_rel_decay_correction = _get_time_value(u_rel_decay_correction)
    except KeyError:
        u_rel_decay_correction = 0.0
    
    try:
        u_rel_dead_time = df.loc[('Uncertainty Budget', 'u_rel Dead Time')]
        u_rel_dead_time = _get_time_value(u_rel_dead_time)
    except KeyError:
        u_rel_dead_time = 0.0
    
    try:
        u_rel_coinc_w = df.loc[('Uncertainty Budget', 'u_rel Coinc Window')]
        u_rel_coinc_w = _get_time_value(u_rel_coinc_w)
    except KeyError:
        u_rel_coinc_w = 0.0
        
    try:
        birks = df.loc[('Model Stochastic', 'Birks Constant (kB)')]
        float(birks.Value)
        birks = birks.Value
    except (KeyError, TypeError):
        return None
    
    #efficiency
    #despite different efficiencies are reported in the csv file
    #only the one evaluated through the Stochastic model in considered
    MSefficiencyD, MSefficiencyD_u = _get_activity(df, 'Model Stochastic', secondary_key='Efficiency')
        
    if np.isnan(MSefficiencyD):
        return None
    
    #Calculated activity at measurement time
    Cactivity = rates_Doubles / MSefficiencyD
    #GUF for uncertainty propagation
    Cactivity_u = Cactivity * np.sqrt(np.power(rates_Doubles_u/rates_Doubles,2) + np.power(MSefficiencyD_u/MSefficiencyD,2))
    
    data = {'nuclide':nuclide, 'source':sourcelab, 'meas_date':meas_date, 'meas_time':meas_time, 'meas_dead_time':meas_dead_time, 'device':device, 'ext_dead_time':dead_time, 'coincidence_w':coincidence_w, 'ndfilter':ndfilter, 'cocktail':cocktail, 'volume':volume, 
    'rates_DoublesRaw':rates_DoublesRaw, 'rates_DoublesRaw_u':rates_DoublesRaw_u, 'rates_TriplesRaw':rates_TriplesRaw, 'rates_TriplesRaw_u':rates_TriplesRaw_u,
    'rates_Doubles':rates_Doubles, 'rates_Doubles_u':rates_Doubles_u, 'rates_Triples':rates_Triples, 'rates_Triples_u':rates_Triples_u, 'TDCR':TDCR, 'TDCR_u':TDCR_u, 'MSefficiencyD':MSefficiencyD, 'MSefficiencyD_u':MSefficiencyD_u, 'C_Activity':Cactivity, 'C_Activity_u':Cactivity_u,
    'birks':birks, 'filename':filename}
    
    return pd.Series(data)

def _get_filenames(folder):
    """Retrieve filepaths of files in the corresponding folder (down to 2 levels)
    
    Parameters
    ----------
    folder : str (Path)
        directory name where measurement files are found, subdirectories are inspected recursively for 2 levels 
    
    Return
    ------
    second_level_names : list
        list of filepaths from the folder directory
    """
    names = [os.path.join(folder, name) for name in os.listdir(folder)]
    
    first_level_names = []
    for name in names:
        if not os.path.isfile(name):
            lower_level_names = [os.path.join(name, subname) for subname in os.listdir(name)]
            first_level_names += lower_level_names
        else:
            first_level_names.append(name)
            
    second_level_names = []
    for name in first_level_names:
        if not os.path.isfile(name):
            lower_level_names = [os.path.join(name, subname) for subname in os.listdir(name)]
            second_level_names += lower_level_names
        else:
            second_level_names.append(name)
            
    #cleanup (filter folders and return only files)
    second_level_names = [filename for filename in second_level_names if os.path.isfile(filename)]

    return second_level_names
    
def _get_birks_best_value(data, autoplot=False):
    """Inspect the data and find the best value of birks constant from measurements performed with different ndfilters
    
    Parameters
    ----------
    data : pandas.DataFrame
        original dataset
    autoplot : bool
        whether automatically display a plot (default False)
    
    Return
    ------
    birks_labels : pandas.Series
        Series identifying the label for the best birks constant value
    (data_results, data_uncertainties) : tuple
        tuple containing elaborated results and uncertainties
    """
    #sort dataframe in chronological order (just to be sure)
    data.sort_values('meas_date', axis=0, ascending=True, inplace=True, ignore_index=False)
    
    _birksindexes = tuple(product(data['ext_dead_time'].unique(), data['coincidence_w'].unique(), data['cocktail'].unique(), data['source'].unique()))
    
    data_results = {}
    data_uncertainties = {}
    
    for item in _birksindexes:
        _filter_with_ndfilter = (data['ext_dead_time'] == item[0]) & (data['coincidence_w'] == item[1]) & (data['cocktail'] == item[2]) & (data['source'] == item[3]) & (data['ndfilter'] != '0')
        _filter_without_ndfilter = (data['ext_dead_time'] == item[0]) & (data['coincidence_w'] == item[1]) & (data['cocktail'] == item[2]) & (data['source'] == item[3]) & (data['ndfilter'] == '0')
        if np.sum(_filter_with_ndfilter) > 0:
            birks_subsample = data[_filter_with_ndfilter]
            normalization_subsample = data[_filter_without_ndfilter]
            
            if autoplot:
                PlotBirks = visualization.Plotbirks(title=f'EDT: {item[0]}, COI: {item[1]}, CT: {item[2]}, ID: {item[3]}')
            
            for birks_value in birks_subsample['birks'].unique():
                
                ndfilter_selected_birks = birks_subsample[birks_subsample['birks'] == birks_value]
                withoutndfilter_selected_birks = normalization_subsample[normalization_subsample['birks'] == birks_value]
                
                #arrays approach
                series_measurement_date = []
                series_calc_activity = []
                series_calc_uncertainty = []
                series_tdcr = []
                for _ndfilter in ndfilter_selected_birks['ndfilter'].unique():
                    #disposable lists for ordering relevant data
                    _md = []
                    _ca = []
                    _cu = []
                    _td = []
                    for measurement_date, calc_activity, calc_uncertainty, TDRC_value in zip(ndfilter_selected_birks[ndfilter_selected_birks['ndfilter'] == _ndfilter]['meas_date'], ndfilter_selected_birks[ndfilter_selected_birks['ndfilter'] == _ndfilter]['C_Activity'], ndfilter_selected_birks[ndfilter_selected_birks['ndfilter'] == _ndfilter]['C_Activity_u'], ndfilter_selected_birks[ndfilter_selected_birks['ndfilter'] == _ndfilter]['TDCR']):
                        
                        _md.append(measurement_date)
                        _ca.append(calc_activity)
                        _cu.append(calc_uncertainty)
                        _td.append(TDRC_value)

                    series_measurement_date.append(_md)
                    series_calc_activity.append(_ca)
                    series_calc_uncertainty.append(_cu)
                    series_tdcr.append(_td)

                #normalization row
                _md = []
                _ca = []
                _cu = []
                _td = []
                for _ in range(len(series_measurement_date[0])):
                    #average date to find the closest 0nd filter datum (previous date is checked first than closest)
                    dates_ = [row[_] for row in series_measurement_date]
                    target_measurement_date = dates_[0] + (np.max(dates_) - np.min(dates_))/2

                    #find the previous measurement date among the data taken without ndfilter
                    support_column = withoutndfilter_selected_birks['meas_date'] - target_measurement_date
                    support_column = support_column[support_column < pd.Timedelta(days=0)]
                    if len(support_column) > 0:
                        support_column.abs().idxmin()
                        previous_idx = support_column.abs().idxmin()
                        selected_row = withoutndfilter_selected_birks.loc[previous_idx]
                    
                    else:
                        #find the closest measurement date among the data taken without ndfilter
                        closest_idx = (withoutndfilter_selected_birks['meas_date'] - target_measurement_date).abs().idxmin()
                        selected_row = withoutndfilter_selected_birks.loc[closest_idx]
                    
                    _md.append(selected_row['meas_date'])
                    _ca.append(selected_row['C_Activity'])
                    _cu.append(selected_row['C_Activity_u'])
                    _td.append(selected_row['TDCR'])
                    
                series_measurement_date.append(_md)
                series_calc_activity.append(_ca)
                series_calc_uncertainty.append(_cu)
                series_tdcr.append(_td)
                
                series_measurement_date, series_calc_activity, series_calc_uncertainty, series_tdcr = np.array(series_measurement_date), np.array(series_calc_activity), np.array(series_calc_uncertainty), np.array(series_tdcr)

                for meas_date_col, C_Activity_col, C_Activity_u_col, TDCR_col in zip(series_measurement_date.T, series_calc_activity.T, series_calc_uncertainty.T, series_tdcr.T):
                    average_measurement_date = meas_date_col.min() + (meas_date_col.max() - meas_date_col.min())/2
                    
                    if (item[0], item[1], item[2], average_measurement_date) not in data_results.keys():
                        data_results[(item[0], item[1], item[2], average_measurement_date)] = {}
                        data_uncertainties[(item[0], item[1], item[2], average_measurement_date)] = {}
                
                    #calculations
                    r_val, r_unc = _linear_fitting_procedure_birks(TDCR_col, C_Activity_col / C_Activity_col[-1], C_Activity_u_col / C_Activity_col[-1])
                    data_results[(item[0], item[1], item[2], average_measurement_date)][birks_value] = r_val
                    data_uncertainties[(item[0], item[1], item[2], average_measurement_date)][birks_value] = r_unc
                    
                    slope = (1 - r_val) / (np.max(TDCR_col) - 0)
                    residuals = C_Activity_col / C_Activity_col[-1] - (slope * TDCR_col + r_val)
                    
                    if autoplot:
                        PlotBirks.add_dataset(TDCR_col, C_Activity_col / C_Activity_col[-1], C_Activity_u_col / C_Activity_col[-1], r_val, r_unc, residuals, label=birks_value)

            if autoplot:
                PlotBirks.show()

    data_results = pd.DataFrame(data_results).T.sort_index()
    data_uncertainties = pd.DataFrame(data_uncertainties).T.sort_index()

    #which birks minimizes
    birks_labels = np.abs(data_results - 1).idxmin(axis=1)

    return birks_labels, (data_results, data_uncertainties)
    
def renormalize_data(reduced_data, half_life=None, half_life_uncertainty=None):
    """Normalize the dataset based on the updated half-life value
    
    Parameters
    ----------
    reduced_data : pandas.DataFrame
        reduced dataframe to be updated and normalized
    half_life : float
        half-life value in d for decay correction; skipped if None (default None)
    half_life_uncertainty : float
        half-life standard uncertainty in d for decay correction uncertainty; skipped if None (default None)
    
    Return
    ------
    reduced_data_normalized : pandas.DataFrame
        reduced dataframe normalized
    """
    reduced_data['meas_time'] = reduced_data['meas_time'].astype(float)
    #normalized data
    if half_life is not None and half_life_uncertainty is not None:
        decay_const = np.log(2) / (half_life * 86400)
        reduced_data['decay_corr'] = (1 - np.exp(-decay_const * reduced_data['meas_time'])) / (decay_const * reduced_data['meas_time'])
        #linear numeric approximation (sensitivity coefficient)
        decay_const_p = np.log(2) / ((half_life + half_life_uncertainty) * 86400)
        decay_const_m = np.log(2) / ((half_life - half_life_uncertainty) * 86400)
        cs_array = ((1 - np.exp(-decay_const_p * reduced_data['meas_time'])) / (decay_const_p * reduced_data['meas_time']) - (1 - np.exp(-decay_const_m * reduced_data['meas_time'])) / (decay_const_m * reduced_data['meas_time'])) / (2 * half_life_uncertainty * 86400 + 1E-24)
        reduced_data['decay_corr_u'] = cs_array * half_life_uncertainty * 86400
    else:
        reduced_data['decay_corr'] = [1.0 for _ in range(len(reduced_data))]
        reduced_data['decay_corr_u'] = [0.0 for _ in range(len(reduced_data))]
    reduced_data['F_Act'] = reduced_data['Act'] / reduced_data['decay_corr']
    reduced_data['F_Act_u'] = reduced_data['F_Act'] * (np.power(reduced_data['Act_u']/reduced_data['Act'], 2) + np.power(reduced_data['decay_corr_u']/reduced_data['decay_corr'], 2))**0.5
    reduced_data['norm_Act'] = [1.0 for _ in range(len(reduced_data))]
    reduced_data['norm_date'] = reduced_data['meas_date']

    _indexes = tuple(product(reduced_data['ext_dead_time'].unique(), reduced_data['coincidence_w'].unique(), reduced_data['cocktail'].unique(), reduced_data['source'].unique()))
    
    for item in _indexes:
        _filter = (reduced_data['ext_dead_time'] == item[0]) & (reduced_data['coincidence_w'] == item[1]) & (reduced_data['cocktail'] == item[2]) & (reduced_data['source'] == item[3])
        if np.sum(_filter) == 1:
            reduced_data.loc[_filter,'norm_Act'] = reduced_data.loc[_filter,'F_Act']
            reduced_data.loc[_filter,'norm_date'] = reduced_data.loc[_filter,'meas_date']
        elif np.sum(_filter) > 0:
            earliest_datum = reduced_data[_filter]['meas_date'].idxmin(axis='index')
            latest_datum = reduced_data[_filter]['meas_date'].idxmax(axis='index')
            reduced_data.loc[_filter,'norm_Act'] = np.exp(np.log(reduced_data.loc[latest_datum,'F_Act']) + (np.log(reduced_data.loc[earliest_datum,'F_Act']) - np.log(reduced_data.loc[latest_datum,'F_Act']))/2)
            reduced_data.loc[_filter,'norm_date'] = reduced_data.loc[earliest_datum,'meas_date'] + (reduced_data.loc[latest_datum,'meas_date'] - reduced_data.loc[earliest_datum,'meas_date'])/2
            
    reduced_data['meas_date'] = pd.to_datetime(reduced_data['meas_date'])
    reduced_data['norm_date'] = pd.to_datetime(reduced_data['norm_date'])
    
    return reduced_data

def get_HL_data_from_dir(directory, nuclide=None, autoplot=False):
    """Retrieve relevant half-life information from csv files in a directory
    
    Parameters
    ----------
    directory : str (Path)
        directory name where measurement files are found
    nuclide : str
        safety parameter to avoid conflicts due to multiple nuclides being present (default None)
    
    Return
    ------
    data : pandas.DataFrame
        whole dataframe
    reduced_data : pandas.DataFrame
        data reduction of measurements taken in the same conditions
    information : dict
        dictionary containing useful information about input files and birks elaboration
    """
    
    files = _get_filenames(directory)
    
    datafiles = []
    rejected_files = []
    
    for file in files:
        data_from_file = open_result_file(file)
        if data_from_file is not None:
            datafiles.append(data_from_file)
        else:
            rejected_files.append(file)
    
    data = pd.concat(datafiles, axis=1).T
    information = {}
    information[('input', 'open errors')] = rejected_files

    if len(data) < len(files):
        print(f'{len(files) - len(data)} were rejected!')

    if len(data['nuclide'].unique()) > 1 and isinstance(nuclide, str):
        nuclide_filter = data['nuclide'] == nuclide
        data = data[nuclide_filter]
    elif len(data['nuclide'].unique()) > 1:
        raise KeyError('WARNING! multiple nuclides in the dataset, use the keyword argument nuclide to select the one to keep')

    nuclide_name = data.iloc[0,0]
    data['meas_date'] = pd.to_datetime(data['meas_date'])
    
    #find best birks constant value
    birks_labels, _ = _get_birks_best_value(data, autoplot=autoplot)
    information[('birks', 'selected values')] = birks_labels
    information[('birks', 'intercept values')] = _[0]
    information[('birks', 'intercept uncertainties')] = _[1]

    #measurement data
    #data reduction
    #filter measurements taken without ndfilter
    measurement_data = data[data['ndfilter'] == '0']
    _indexes = tuple(product(measurement_data['meas_date'].unique(), measurement_data['ext_dead_time'].unique(), measurement_data['coincidence_w'].unique(), measurement_data['cocktail'].unique(), measurement_data['source'].unique()))
    #group data from the same measurement (same start date and settings, nd-filter)
    
    distilled_data = []
    for item in _indexes:
        _filter = (measurement_data['meas_date'] == item[0]) & (measurement_data['ext_dead_time'] == item[1]) & (measurement_data['coincidence_w'] == item[2]) & (measurement_data['cocktail'] == item[3]) & (measurement_data['source'] == item[4])
        if np.sum(_filter) > 0:
            subsample = measurement_data[_filter]
            
            #selection of line based on the birks value
            try:
                key_birks = birks_labels.loc[item[1], item[2], item[3]]
            except KeyError:
                #no cocktail condition, assumes birks is the same regardless of cocktail used
                key_birks = birks_labels.loc[item[1], item[2]]

            label_birks = ''
            if len(key_birks) == 1:
                label_birks = key_birks[key_birks.index[0]]
            else:
                iloc_idx = key_birks.index.get_indexer([item[0]], method='nearest')
                label_birks = key_birks.iloc[iloc_idx].values[0]
            
            subsample = subsample[subsample['birks'] == label_birks]

            Act_value = subsample['C_Activity'].values[0]
            Act_unc = subsample['C_Activity_u'].values[0]
            fname = subsample['filename'].values[0]

            condensed_info = {'nuclide':nuclide_name, 'meas_date':item[0], 'meas_time':subsample['meas_time'].values[0], 'ext_dead_time':item[1], 'coincidence_w':item[2], 'ndfilter':'0', 'cocktail':item[3], 'volume':subsample['volume'].values[0], 'source':item[4], 'Act':Act_value, 'Act_u':Act_unc, 'birks':label_birks, 'filename':fname}
            
            distilled_data.append(pd.Series(condensed_info))
    
    reduced_data = pd.concat(distilled_data, axis=1).T
    reduced_data = renormalize_data(reduced_data)
    
    return data, reduced_data, information

def _exp(x, y, _lambda):
    """Exponential function accounting for radioactive decay
    
    Parameters
    ----------
    x : numpy.array or float
        independent variable
    y : float
        intercept
    _lambda : float
        decay constant

    Return
    ------
    Y : numpy.array or float
        result of exponential funtion 'y * np.exp(-_lambda * x)'
    """
    return y * np.exp(-_lambda * x)

def _exponential_fitting_procedure(x, y, uy=None, autoplot=False, title=''):
    """Perform exponential fit on the dataset (x, lny) to obtain estimate of the decay constant (slope)
    
    Parameters
    ----------
    x : numpy.array
        array containing a (0-indexed) normalized time series
    y : numpy.array
        array containing the (0-indexed) normalized activity
    uy : numpy.array
        array containing the (0-indexed) normalized activity uncertainties (default None)
    autoplot : bool
        whether automatically display a plot (default False)
    title : str
        title of the plot (default '')

    Return
    ------
    estimated_halflife : float
        half-life value (in days) resulting from performed fit
    estimated_uncertainty : float
        standard uncertaity associated with fitted half-life
    residuals : numpy.array
        residuals of the fit
    """
    if uy is not None:
        weig = 'weighted '
    else:
        weig = ''
    
    popt, pcov = curve_fit(_exp, x, y, p0=[y[0], 0], sigma=uy)
    np.sqrt(np.diag(pcov))[1]
    
    estimated_halflife = np.log(2)/popt[1]
    estimated_uncertainty = estimated_halflife * np.sqrt(np.diag(pcov))[1]/popt[1]
    
    residuals = y - _exp(x, *popt)
    
    if autoplot:
        fit_x = np.linspace(np.min(x), np.max(x), 1000)
        fit_y = _exp(fit_x, *popt)
        if title == '':
            stitle = f'{weig}exponential fit'
        else:
            stitle = title + f'\n{weig}exponential fit'
        visualization._fit_plot(x, y, fit_x, fit_y, residuals, uy=uy, suptitle=stitle, fitted_HL=f'{estimated_halflife:.1f}')

    return estimated_halflife, estimated_uncertainty, residuals
    
def _linear_fitting_procedure_M(X, Y, UY=None, autoplot=False, title='', k_limit=2.5):
    """Perform linear fit on the dataset (x, lny) to obtain estimate of the decay constant (slope) but allows outliers rejection
    
    Parameters
    ----------
    X : numpy.array
        array containing a normalized time series
    Y : numpy.array
        array containing the natural log of a normalized activity
    UY : numpy.array
        array containing the normalized activity uncertainties at different TDCR values (default None)
    autoplot : bool
        whether automatically display a plot (default False)
    title : str
        title of the plot (default '')
    k_limit : float
        number of standard deviations to check for statistical consistency (default 2.5)

    Return
    ------
    estimated_halflife : float
        half-life value (in days) resulting from performed fit
    estimated_uncertainty : float
        standard uncertainty associated with fitted half-life
    masque : list
        list of bool indicationg True at the index of elaborated data and False at the index of rejected ones
    """
    
    masque = [True] * len(X)
    satisfaction = False
    
    while not satisfaction:
        MX = np.array([xi for xi,mi in zip(X,masque) if mi == True])
        MY = np.array([yi for yi,mi in zip(Y,masque) if mi == True])
        
        #creation of independent variables matrix
        W = MX[:, np.newaxis]**[0,1]
    
        if UY is not None:
            MUY = np.array([uyi for uyi,mi in zip(UY,masque) if mi == True])
            V = np.linalg.inv(np.diag(np.power(MUY, 2)))#weighted on variances
        else:
            V = np.identity(W.shape[0])

        #parameters estimation
        parameters = np.linalg.inv(W.T@V@W)@(W.T@V@MY)
    
        #residuals
        residuals = MY - parameters@W.T
    
        #covariance matrix
        n, k = MY.shape[0], W.shape[1]
        mcov = np.true_divide(1, n-k) * (residuals.T@V@residuals) * np.linalg.inv(W.T@V@W)

        #outliers based on residuals and standard uncertainties
        if UY is not None:
            de = k_limit * MUY
        else:
            de = np.exp(parameters[0]) * np.sqrt(np.diag(mcov))[0] * k_limit
        _norm_res = np.abs(residuals/de)
        check = _norm_res > 1

        if np.sum(check) > 0:
            if np.sum(masque) < 4:
                return np.nan, np.nan, masque
            
            idx_f = _norm_res.argmax()

            #updated mask
            pos = None
            for _nn, item in enumerate(masque):
                if item == True:
                    try:
                        pos += 1
                    except TypeError:
                        pos = 0

                if pos == idx_f:
                    masque[_nn] = False
                    break
        else:
            satisfaction = True

    estimated_halflife = -np.log(2)/parameters[1]
    estimated_uncertainty = estimated_halflife * np.sqrt(np.diag(mcov))[1]
    
    if autoplot:
        fit_x = np.linspace(np.min(X), np.max(X), 1000)
        W2 = fit_x[:, np.newaxis]**[0,1]
        fit_y = parameters@W2.T
        if title == '':
            stitle = 'weighted linear fit'
        else:
            stitle = title + '\nweighted linear fit'
        visualization._fit_plot_M(X, Y, fit_x, fit_y, residuals, masque, uy=UY, suptitle=stitle, fitted_HL=f'{estimated_halflife:.1f}')
    
    return estimated_halflife, estimated_uncertainty, masque, residuals #return residuals too and see what happens

def _linear_fitting_procedure(X, Y, UY=None, autoplot=False, title=''):
    """Perform linear fit on the dataset (x, lny) to obtain estimate of the decay constant (slope)
    
    Parameters
    ----------
    X : numpy.array
        array containing a normalized time series
    Y : numpy.array
        array containing the natural log of a normalized activity
    UY : numpy.array
        array containing the normalized activity uncertainties at different TDCR values (default None)
    autoplot : bool
        whether automatically display a plot (default False)
    title : str
        title of the plot (default '')

    Return
    ------
    estimated_halflife : float
        half-life value (in days) resulting from performed fit
    estimated_uncertainty : float
        standard uncertaity associated with fitted half-life
    residuals : numpy.array
        residuals of the fit
    """
    
    #creation of independent variables matrix
    W = X[:, np.newaxis]**[0,1]
    
    if UY is not None:
        V = np.linalg.inv(np.diag(np.power(UY, 2)))#weighted on variances
        weig = 'weighted '
    else:
        V = np.identity(W.shape[0])
        weig = ''

    #parameters estimation
    parameters = np.linalg.inv(W.T@V@W)@(W.T@V@Y)
    
    #residuals
    residuals = Y - parameters@W.T
    
    #covariance matrix
    n, k = Y.shape[0], W.shape[1]
    mcov = np.true_divide(1, n-k) * (residuals.T@V@residuals) * np.linalg.inv(W.T@V@W)

    estimated_halflife = -np.log(2)/parameters[1]
    estimated_uncertainty = estimated_halflife * np.sqrt(np.diag(mcov))[1]
    
    if autoplot:
        fit_x = np.linspace(np.min(X), np.max(X), 1000)
        W2 = fit_x[:, np.newaxis]**[0,1]
        fit_y = parameters@W2.T
        if title == '':
            stitle = f'{weig}linear fit'
        else:
            stitle = title + f'\n{weig}linear fit'
        visualization._fit_plot(X, Y, fit_x, fit_y, residuals, uy=UY, suptitle=stitle, fitted_HL=f'{estimated_halflife:.1f}')
    
    return estimated_halflife, estimated_uncertainty, residuals
    
def _montecarlo_fitting_procedure(X, Y, UY, N=1000, masque=None, autoplot=False, title='', linear=True):
    """Perform fit on the dataset (x, y) N times from random normal draws
    
    Parameters
    ----------
    X : numpy.array
        array containing a normalized time series
    Y : numpy.array
        array containing the normalized activity
    UY : numpy.array
        array containing the normalized activity uncertainties
    N : int
        number of trials (default 1000)
    masque : list
        list of the same lenght of x indicating data to be rejected (default None)
    autoplot : bool
        whether automatically display a plot (default False)
    title : str
        title of the plot (default '')
    linear : bool
        whether perform a linear or an exponential fit (default True)

    Return
    ------
    estimated_halflife : float
        half-life value (in days) resulting from performed fit
    estimated_uncertainty : float
        standard uncertaity associated with fitted half-life
    residuals : numpy.array
        averaged residuals of the fit
    """

    results = []
    residuals = []
    
    if masque is not None:
        X = X[masque]
        Y = Y[masque]
        UY = UY[masque]
    
    for _ in range(N):
        ry = np.random.normal(loc=Y, scale=UY)

        if linear:
            i_halflife, i_uncertainty, _ = _linear_fitting_procedure(X, np.log(ry))
        else:
            i_halflife, i_uncertainty, _ = _exponential_fitting_procedure(X - X[0], ry)
        
        results.append(i_halflife)
        residuals.append(_)

    residuals = np.average(residuals, axis=0)
    residuals = residuals.flatten()

    estimated_halflife, estimated_uncertainty = np.average(results), np.std(results)

    if autoplot:
        if linear:
            _fit = 'linear '
        else:
            _fit = 'exponential '
        
        if title == '':
            stitle = f'montecarlo {_fit}fit'
        else:
            stitle = title + f'\nmontecarlo {_fit}fit'
        visualization._distribution_plot(results, suptitle=stitle, fitted_HL=f'{estimated_halflife:.1f}')
    
    return estimated_halflife, estimated_uncertainty, residuals

def _get_autocorrelation(residuals, alpha=0.05):
    """Perform autocorrelation calculation using Ljung–Box procedure on the residuals of fit
    
    Parameters
    ----------
    residuals : numpy.array
        evalauted residuals of fit
    alpha : float
        significance level (default 0.05)
    
    Return
    ------
    autocorrelation : bool
        whether autocorrelation is identified within the input residuals
    """
    res_res = residuals - np.mean(residuals)
    res_norm = np.sum(res_res**2)
    auto_corr = np.correlate(res_res, res_res, mode='full') / res_norm
    auto_corr = auto_corr[int(len(auto_corr)/2):]

    #autocorrelation calculation with Ljung–Box
    kappa = np.arange(1, len(residuals), 1)
    Q = len(residuals) * (len(residuals) + 2) * np.sum(np.power(auto_corr[1:], 2) / (len(residuals) - kappa))

    autocorrelation = chi2.cdf(Q, len(residuals)-1) > 1 - alpha
    
    return autocorrelation

def fit_data(dataset, data_threshold=3, MC_trials=10000, autoplot=False, fit='all'):#options...
    """Perform fitting on the dataset to obtain a compilation of evaluated half-lives of every sub-dataset
    
    Parameters
    ----------
    dataset : pandas.DataFrame
        dataframe containing data to fit
    data_threshold : int
        minimum number of data to perform the fit (default 3)
    MC_trials : int
        number of MC trials for the MonteCarlo fitting method (default 10000)
    autoplot : bool
        whether automatically display a plot (default False)
    fit : str
        string defining a specific fit to perform (default 'all')
        accepted input (case insensitive):
            'weighted linear', 'wl'                     to perform Weighted Linear fit
            'linear', 'l'                               to perform Linear fit
            'weighted exponential', 'wexp'              to perform Weighted Exponential fit
            'exponential', 'exp'                        to perform Exponential fit
            'montecarlo linear', 'mcl'                  to perform MonteCarlo Linear fit
            'montecarlo exponential', 'mcexp'           to perform MonteCarlo Exponential fit
            'montecarlo', 'mcm'                         to perform only MonteCarlo fits
            'weighted', 'w'                             to perform only weighted fits
            'nonweighted', 'nw'                         to perform only non-weighted fits
            'all'                                       to perform all the fits
    
    Return
    ------
    result : pandas.DataFrame
        result of the performed fit
    """
    
    #fallback to 'all' if invalid string
    if fit.lower() not in ('weighted linear', 'wl', 'montecarlo', 'mcm', 'weighted', 'w',
    'linear', 'l', 'weighted exponential', 'wexp', 'exponential', 'exp', 'nonweighted', 'nw',
    'montecarlo linear', 'mcl', 'montecarlo exponential', 'mcexp', 'all'):
        fit = 'all'
    
    #get normalized times (in days) array
    n_times = np.array([(MEAS - NORM).total_seconds()/86400 for MEAS, NORM in zip(dataset['meas_date'], dataset['norm_date'])])
    
    #divide the dataset into multiple subsets
    #depending on the experimental parameters:
    # 'ext_dead_time'
    # 'coincidence_w'
    # 'cocktail'
    # 'source'
    _indexes = tuple(product(dataset['ext_dead_time'].unique(), dataset['coincidence_w'].unique(), dataset['cocktail'].unique(), dataset['source'].unique()))
    
    series = []
    fit_results = {}
    information = {}
    
    for item in _indexes:
        _filter = (dataset['ext_dead_time'] == item[0]) & (dataset['coincidence_w'] == item[1]) & (dataset['cocktail'] == item[2]) & (dataset['source'] == item[3])
        if np.sum(_filter) > data_threshold:#minimum number of data to perform the fit
            local_dataset = dataset[_filter]
            local_time = n_times[_filter]
            setup = f'EDT: {item[0]}, COI: {item[1]}\nCT: {item[2]}, ID: {item[3]}, BIRKS: {local_dataset['birks'].values[0]}'
            y, uy, ruy = (local_dataset['Act']/local_dataset['norm_Act']).to_numpy(dtype='float64'), (local_dataset['Act_u']/local_dataset['norm_Act']).to_numpy(dtype='float64'), (local_dataset['Act_u']/local_dataset['Act']).to_numpy(dtype='float64')
            #lazy data rejection based on weighted linear fit and used for every other fitting procedure
            #weighted linear fit
            if fit.lower() not in ('weighted linear', 'wl', 'weighted', 'w', 'all'):
                local_apt = False
            else:
                local_apt = autoplot
            linmodel_halflife, linmodel_uncertainty, masque, linmodel_residuals = _linear_fitting_procedure_M(local_time, np.log(y), ruy, autoplot=local_apt, title=setup)

            #additional contribution to uncertainties other than that of the fit are grouped in high frequency '_uHF', medium frequency '_uMF' and low frequency '_uLF'
            #evaluation of high frequency uncertainty (Pomme), constant part (assuming fixed relative uncertainty)
            if masque is not None:
                total_time = local_time[masque]
            else:
                total_time = local_time
            _sum = np.sum(np.abs(total_time / np.max(total_time) - 1))
            constant_part = 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(1 / _sum) * np.mean(ruy[masque])
            
            #evaluation of medium frequency uncertainty (Pomme) through autocorrelation
            #also source positioning for measurement? (it should be statistic as well)

            #low frequency
            #stability of the solution, background correction?
            
            #choice of fit to perform
            if fit.lower() in ('weighted linear', 'wl', 'weighted', 'w', 'all'):
                linmodel_autocorrelation = _get_autocorrelation(linmodel_residuals)
                fit_results['wl_HL'] = linmodel_halflife
                fit_results['wl_uHL'] = linmodel_uncertainty
                fit_results['wl_uHF'] = np.power(linmodel_halflife, 2) * constant_part
                if linmodel_autocorrelation:
                    #manage uncertainty contribution from medium frequency effects
                    fit_results['wl_uMF'] = np.power(linmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(linmodel_residuals)+1)) * np.mean(ruy[masque]) # test, Pomme
                else:
                    fit_results['wl_uMF'] = 0.0
                fit_results['wl_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'weighted linear fit', 'autocorrelation')] = (bool(linmodel_autocorrelation), 1-0.05, linmodel_residuals)
            
            #linear fit
            if fit.lower() in ('linear', 'l', 'nonweighted', 'nw', 'all'):
                alinmodel_halflife, alinmodel_uncertainty, alinmodel_residuals = _linear_fitting_procedure(local_time[masque], np.log(y[masque]), autoplot=autoplot, title=setup)
                fit_results['l_HL'] = alinmodel_halflife
                fit_results['l_uHL'] = alinmodel_uncertainty
                fit_results['l_uHF'] = np.power(alinmodel_halflife, 2) * constant_part
                alinmodel_autocorrelation = _get_autocorrelation(alinmodel_residuals)
                if alinmodel_autocorrelation:
                    #manage uncertainty contribution from medium frequency effects
                    fit_results['l_uMF'] = np.power(alinmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(alinmodel_residuals)+1)) * np.mean(ruy[masque])
                else:
                    fit_results['l_uMF'] = 0.0
                fit_results['l_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'linear fit', 'autocorrelation')] = (bool(alinmodel_autocorrelation), 1-0.05, alinmodel_residuals)
            
            #weighted exponential fit
            if fit.lower() in ('weighted exponential', 'wexp', 'weighted', 'w', 'all'):
                wexpmodel_halflife, wexpmodel_uncertainty, wexpmodel_residuals = _exponential_fitting_procedure(local_time[masque] - local_time[masque][0], y[masque], uy[masque], autoplot=autoplot, title=setup)
                fit_results['wexp_HL'] = wexpmodel_halflife
                fit_results['wexp_uHL'] = wexpmodel_uncertainty
                fit_results['wexp_uHF'] = np.power(wexpmodel_halflife, 2) * constant_part
                wexpmodel_autocorrelation = _get_autocorrelation(wexpmodel_residuals)
                if wexpmodel_autocorrelation:
                    #manage uncertainty contribution from medium frequency effects
                    fit_results['wexp_uMF'] = np.power(wexpmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(wexpmodel_residuals)+1)) * np.mean(ruy[masque])
                else:
                    fit_results['wexp_uMF'] = 0.0
                fit_results['wexp_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'weighted exponential fit', 'autocorrelation')] = (bool(wexpmodel_autocorrelation), 1-0.05, wexpmodel_residuals)
                
            #exponential fit
            if fit.lower() in ('exponential', 'exp', 'nonweighted', 'nw', 'all'):
                expmodel_halflife, expmodel_uncertainty, expmodel_residuals = _exponential_fitting_procedure(local_time[masque] - local_time[masque][0], y[masque], autoplot=autoplot, title=setup)
                fit_results['exp_HL'] = expmodel_halflife
                fit_results['exp_uHL'] = expmodel_uncertainty
                fit_results['exp_uHF'] = np.power(expmodel_halflife, 2) * constant_part
                expmodel_autocorrelation = _get_autocorrelation(expmodel_residuals)
                if expmodel_autocorrelation:
                    fit_results['exp_uMF'] = np.power(expmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(expmodel_residuals)+1)) * np.mean(ruy[masque])
                else:
                    fit_results['exp_uMF'] = 0.0
                fit_results['exp_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'exponential fit', 'autocorrelation')] = (bool(expmodel_autocorrelation), 1-0.05, expmodel_residuals)
                
            #montecarlo linear fit
            if fit.lower() in ('montecarlo linear', 'mcl', 'montecarlo', 'mcm', 'all'):
                mcmodel_halflife, mcmodel_uncertainty, mcmodel_residuals = _montecarlo_fitting_procedure(local_time, y, uy, autoplot=autoplot, masque=masque, N=MC_trials, title=setup, linear=True)
                fit_results['mcl_HL'] = mcmodel_halflife
                fit_results['mcl_uHL'] = mcmodel_uncertainty
                fit_results['mcl_uHF'] = 0.0 #montecarlo should already consider this contribution
                mcmodel_autocorrelation = _get_autocorrelation(mcmodel_residuals)
                if mcmodel_autocorrelation:
                    fit_results['mcl_uMF'] = np.power(mcmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(mcmodel_residuals)+1)) * np.mean(ruy[masque])
                else:
                    fit_results['mcl_uMF'] = 0.0
                fit_results['mcl_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'montecarlo linear fit', 'autocorrelation')] = (bool(mcmodel_autocorrelation), 1-0.05, mcmodel_residuals)
                
            #montecarlo exponential fit
            if fit.lower() in ('montecarlo exponential', 'mcexp', 'montecarlo', 'mcm', 'all'):
                mcmodel_halflife, mcmodel_uncertainty, mcmodel_residuals = _montecarlo_fitting_procedure(local_time, y, uy, autoplot=autoplot, masque=masque, N=MC_trials, title=setup, linear=False)
                fit_results['mcexp_HL'] = mcmodel_halflife
                fit_results['mcexp_uHL'] = mcmodel_uncertainty
                fit_results['mcexp_uHF'] = 0.0 #montecarlo should already consider this contribution
                mcmodel_autocorrelation = _get_autocorrelation(mcmodel_residuals)
                if mcmodel_autocorrelation:
                    fit_results['mcexp_uMF'] = np.power(mcmodel_halflife, 2) * 2 / np.log(2) * (1 / (np.max(total_time) * 2)) * np.sqrt(2/(len(mcmodel_residuals)+1)) * np.mean(ruy[masque])
                else:
                    fit_results['mcexp_uMF'] = 0.0
                fit_results['mcexp_uLF'] = 0.0
                information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'montecarlo exponential fit', 'autocorrelation')] = (bool(mcmodel_autocorrelation), 1-0.05, mcmodel_residuals)
            
            information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'accepted datapoints')] = masque
            
            _fdata = {'ext_dead_time':item[0], 'coincidence_w':item[1], 'cocktail':item[2], 'birks':local_dataset['birks'].values[0], 'source':item[3], 
            **fit_results}
            series.append(pd.Series(_fdata))

        else:
            information[(f'{item[0]} {item[1]} {item[2]} {item[3]}', 'accepted datapoints')] = 'DataSeries rejected'

    fitted_data = pd.concat(series, axis=1).T
    return fitted_data, information
    
def PMM_method(x, v, a=None):
    """Return a half-life value (with uncertainty) performing a power moderated mean
    
    Parameters
    ----------
    x : numpy.array (pandas.Series)
        values
    v : numpy.array (pandas.Series)
        variances of x
    a : numeric
        alpha parameter of the PMM method; if None converts to 2-3/N with N=len(x) (default None)
    
    Return
    ------
    PMM_half_life : float
        power moderated mean for half_life
    PMM_uncertainty : float
        power moderated uncertainty
    a : float
        alpha parameter of the PMM method
    wi : numpy.array
        adopted weights
    """

    obs_chisq = 2
    s2 = None
    
    SAFETY_NET = 50000
    _n = 0
    while obs_chisq > 1:
        if s2 is not None:
            s2 += np.average(v) / 100
        else:
            s2 = 0

        w = 1 / (v + s2)
        xmp = np.sum(x * w) / np.sum(w)
        obs_chisq = 1 / (len(x)-1) * np.sum(np.power(x - xmp, 2) / (v + s2))
        
        _n += 1
        if _n == SAFETY_NET:
            print(f'did not converge after {SAFETY_NET} iterations')
            break
    
    #value of a depending on the reliability of the uncertianties
    #a = 0          uncertainties are uninformative or weakly proportional
    #a = 2 - 3/N    uncertainties are informative but with a tendency to be underestimated
    #a = 2          uncertainties are informative and accurate, the dataset is large and consistent
    if not isinstance(a, (int, float)):
        if len(x) > 60:
            a = 2
        else:
            a = 2 - 3/len(x)
    elif a < 0:
        a = 0
    elif a > 2:
        a = 2
    
    x_ave = np.average(x)
    u2x = np.sum(np.power(x - x_ave, 2)) / (len(x) * (len(x) - 1))
    u2xmp = 1 / np.sum(1 / (v + s2))
    S = np.sqrt(len(x) * np.max((u2x, u2xmp)))
    
    var_xref = 1 / np.sum(1 / ((v + s2)**(a/2) * S**(2-a)))
    wi = var_xref / ((v + s2)**(a/2) * S**(2-a))
    
    return np.sum(x*wi), np.sqrt(var_xref), a, wi

def BirgeAdjust(x, v, mod=True):
    """Calculate a reference value using the Birge Adjustment
    
    Parameters
    ----------
    x : numpy.array (pandas.Series)
        values
    v : numpy.array (pandas.Series)
        variances of x
    mod : bool
        whether to use the Modified Birge [Bodnar & Elster 2014 Metrologia 51 516] (default True)
    
    Return
    ------
    M : float
        Birge weighted average
    Birge_uncertainty : float
        Birge Adjustment uncertainty
    sBirge : float
        Birge ratio
    """

    n = len(x)
    sBirge = 1
    N = 1
    sigma = v ** 0.5
    k = 1 / v
    M = np.average(x, weights=k)
    uM = (np.sum(k))**-0.5

    if n > 2:
        sBirge = np.sqrt(1/(n-1) * np.sum(np.power(x-M, 2)/np.power(sigma + 1e-11, 2)))

    if n > 3 and mod:
        N = np.sqrt((n-1)/(n-3))
        #uB = N * sBirge * uM, with N = np.sqrt((n-1)/(n-3))   # Modified Birge [Bodnar & Elster 2014 Metrologia 51 516]
        #uB = N * sBirge * uM, with N = 1                      # unmodified Birge [Birge 1932 Phys. Rev. 40 207]

    return M, N * sBirge * uM, float(sBirge)

def DerSimonianLairdp(x, v):
    """Calculate a reference value using a DerSimonian and Laird procedure
   
    see refs.
    DerSimonian R and Laird N 1986 Meta-analysis in clinical
    trials Control. Clin. Trials 7 177–88
    DerSimonian R and Kacker R 2007 Random-effects model for
    meta-analysis of clinical trials: an update Contemp. Clin.
    Trials 28 105–14
    
    Parameters
    ----------
    x : numpy.array (pandas.Series)
        values
    v : numpy.array (pandas.Series)
        variances of x
    
    Return
    ------
    M : float
        DerSimonian-Laird weighted average
    uM : float
        DerSimonian-Laird uncertainty
    tauDL : float
        additional variance
    """
    
    n = len(x)
    Q = np.sum(1/v * np.power(x - np.average(x), 2))    
    tauM2 = (Q - n + 1) / (np.sum(1/v) - np.sum(v**-2)/np.sum(1/v))
    if tauM2 > 0:
        tauDL = np.sqrt(tauM2)
    else:
        tauDL = 0

    w = 1 / (tauDL**2 + v)
    M = np.average(x, weights=w)
    uM = np.sqrt(1/np.sum(w))

    return M, uM, tauDL

def CoxProcedureA(x, v, alpha=0.05, noFilter=False, k=2.5):
    """Calculate the reference value following the procedure A in Cox 2002 Metrologia 39 589
   
    Parameters
    ----------
    x : numpy.array (pandas.Series)
        values
    v : numpy.array (pandas.Series)
        variances of x
    alpha : float
        risk of first species (defaul = 0.05)
    noFilter: bool
        inactivate or activate the filtering of outliers (default False)
    k : float
        coverage factor for expanded uncertainty on degrees of equivalance (DoE) (default 2.5)

    Return
    ------
    y : float
        reference value
    u_y : float
        standard uncertainty of the reference value
    chi2TestResult : bool
        result of the consistency check
    d : list of floats
        list of DoE
    U_di : list of floats
        expanded uncertainty of the DoE
    _filter : np.array
        array of bool indicating which element was retained for calculation
    """

    N = len(x)
    w = 1 / v

    y = np.sum(x*w) / np.sum(w)
    u_y = 1 / np.sqrt(np.sum(w))

    khi2_obs = np.sum(np.power(x - y, 2) / v)
    chi2TestResult = chi2.cdf(khi2_obs, N-1) < 1-alpha

    d = x - y
    U_di = k * (v - u_y**2)**0.5
    _filter = v >= 0

    if not chi2TestResult and not noFilter:
        _filter = np.abs(d) < U_di
        filtered_x = x[_filter]
        filtered_v = v[_filter]

        w = 1 / filtered_v

        y = np.sum(filtered_x*w) / np.sum(w)
        u_y = 1 / np.sqrt(np.sum(w))

        d = filtered_x - y
        U_di = 2 * (filtered_v - u_y**2)**0.5
    
    return y, u_y, chi2TestResult, d, U_di, _filter

def CoxProcedureB(x, u, M=100000):
    """Calculate the reference value following the procedure B in Cox 2002 Metrologia 39 589

    Parameters
    ----------
    x : numpy.array (pandas.Series)
        values
    u : numpy.array (pandas.Series)
        uncertainties of x
    M : int
        number of gaussian draws (default 100000)

    Return
    ------
    y : float
        reference value
    u_y : float
        standard uncertainty of the reference value
    d : numpy.array
        degrees of equivalence (DoE)
    U_di : numpy.array
        expanded uncertainty of DoE
    """
    
    N = len(x)

    r_norm = np.random.normal(x, u, size=(M, N))
    medianVector = np.median(r_norm, axis=1)

    y = np.mean(medianVector)
    u_y = np.std(medianVector)

    d = x - y
    V = u**2 - u_y**2
    V = np.array([np.max((0, value)) for value in V])
    U_di = 2 * np.sqrt(V)

    return y, u_y, d, U_di
    
def get_result(fitted_data, method='all', iterative=False):
    """Return single half-life values (with uncertainty) considering all the independent results from fitted_data
    
    Parameters
    ----------
    fitted_data : pandas.DataFrame
        file path of the file from which the information are recalled
    method : str
        string defining the method to be adopted to average results (default 'all')
        accepted input (case insensitive):
            'arithmetic average', 'aa'          to perform Arithmetic Average
            'weighted average', 'wa'            to perform Weighted Average
            'power-moderated mean', 'pmm'       to perform Power-Moderated Mean
            'genetic algorithm', 'csg', 'ga'    to perform Genetic Algorithm
            'birge adjustment', 'birge', 'ba'   to perform Birge Adjustment
            'dersimonianlaird', 'dsl'           to perform DerSimonian-Laird procedure
            'cox procedure a', 'coxa'           to perform Cox procedure A
            'cox procedure b', 'coxb'           to perform Cox procedure B
            'all'                               to perform all of the above
    iterative : bool
        whether results are calculated iteratively (default False)

    Return
    ------
    results : dict
        dictionary containing dictionaries for each fitting method including tuples (half_life value, half_life uncertainty) for each averaging method
    information : dict
        dictionary with additional information concerning the averaging procedures
    """

    #fallback to 'all' if invalid string
    if method.lower() not in ('arithmetic average', 'aa', 
    'weighted average', 'wa', 'power-moderated mean', 'pmm', 'genetic algorithm', 
    'csg', 'ga', 'birge adjustment', 'birge', 'ba', 'dersimonianlaird', 'dsl', 
    'cox procedure a', 'coxa', 'cox procedure b', 'coxb', 'all'):
        method = 'all'

    #result columns
    result_columns = list(fitted_data.columns)
    for item in ('ext_dead_time', 'coincidence_w', 'cocktail', 'birks', 'source'):
        try:
            result_columns.remove(item)
        except ValueError:
            pass
    
    list_of_elaborated_data = []
    for item in set([item.split('_')[0] for item in result_columns]):
        if {f'{item}_HL', f'{item}_uHL', f'{item}_uHF', f'{item}_uMF', f'{item}_uLF'}.issubset(result_columns):
            list_of_elaborated_data.append(item)

    elaboration_method_dictionary = {'l':'linear fit', 'wl':'weigthed linear fit', 'wexp':'weighted exponential fit', 'exp':'exponential fit',
    'mcl':'montecarlo linear fit', 'mcexp':'montecarlo exponential fit'}
    results = {}
    information = {}

    #explore different averages
    for elaboration_method in list_of_elaborated_data:
        results[elaboration_method_dictionary[elaboration_method]] = {}
        #calculate combined variances
        fitted_data[f'{elaboration_method}_combined_variance'] = np.power(fitted_data[f'{elaboration_method}_uHL'], 2) + np.power(fitted_data[f'{elaboration_method}_uHF'], 2) + np.power(fitted_data[f'{elaboration_method}_uMF'], 2) + np.power(fitted_data[f'{elaboration_method}_uLF'], 2)
        
        #drop np.nan if any (for safety)
        data_subset = fitted_data.dropna(subset=[f'{elaboration_method}_HL', f'{elaboration_method}_uHL', f'{elaboration_method}_uHF'], inplace=False)
        
        #arithmetic average
        if method.lower() in ('arithmetic average', 'aa', 'all'):
            AA, uAA = np.average(data_subset[f'{elaboration_method}_HL']), np.std(data_subset[f'{elaboration_method}_HL'])
            results[elaboration_method_dictionary[elaboration_method]]['arithmetic average'] = (AA, uAA)
            information[(elaboration_method_dictionary[elaboration_method], 'arithmetic average')] = {'N':len(data_subset[f'{elaboration_method}_HL'])}
    
        #weighted average
        if method.lower() in ('weighted average', 'wa', 'all'):
            w = 1 / data_subset[f'{elaboration_method}_combined_variance']
            WA, uWA = np.sum(data_subset[f'{elaboration_method}_HL'] * w) / np.sum(w), np.sqrt(1 / np.sum(w))
            results[elaboration_method_dictionary[elaboration_method]]['weighted average'] = (WA, uWA)
            information[(elaboration_method_dictionary[elaboration_method], 'weighted average')] = {'weights':w}

        #PMM
        if method.lower() in ('power-moderated mean', 'pmm', 'all'):
            PMM, uPMM, PMM_a, PMM_wi = PMM_method(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance'])
            results[elaboration_method_dictionary[elaboration_method]]['power-moderated mean'] = (PMM, uPMM)
            information[(elaboration_method_dictionary[elaboration_method], 'power-moderated mean')] = {'alpha':PMM_a, 'weights':PMM_wi}

        #consensusGen
        if method.lower() in ('genetic algorithm', 'csg', 'ga', 'all'):
            _resultcsg = csg.consensusGen.consensusGen(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance']**0.5, model='normal')#ng=1, ni=100000, threshold=0.01, model="normal", df=20
            results[elaboration_method_dictionary[elaboration_method]]['genetic algorithm'] = (_resultcsg[0][-1], _resultcsg[1][-1])
            information[(elaboration_method_dictionary[elaboration_method], 'genetic algorithm')] = {'csg[2]':_resultcsg[2], 'csg[3]':_resultcsg[3], 'csg[4]':_resultcsg[4], 'csg[5]':_resultcsg[5], 'csg[6]':_resultcsg[6], 'csg[7]':_resultcsg[7]}

        #Birge Adjustment
        if method.lower() in ('birge adjustment', 'birge', 'ba', 'all'):
            BADJ, uBADJ, sBirge = BirgeAdjust(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance'], mod=True)
            results[elaboration_method_dictionary[elaboration_method]]['Birge adjustment'] = (BADJ, uBADJ)
            information[(elaboration_method_dictionary[elaboration_method], 'Birge adjustment')] = {'Birge parameter':sBirge}
            
        #DerSimonianLairdp
        if method.lower() in ('dersimonianlaird', 'dsl', 'all'):
            DSLp, uDSLp, sDSLp_DL = DerSimonianLairdp(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance'])
            results[elaboration_method_dictionary[elaboration_method]]['DerSimonian-Laird'] = (DSLp, uDSLp)
            information[(elaboration_method_dictionary[elaboration_method], 'DerSimonian-Laird')] = {'DerSimonian-Laird additional variance':sDSLp_DL}
            
        #CoxProcedureA
        if method.lower() in ('cox procedure a', 'coxa', 'all'):
            COXA, uCOXA, COXA_chitest, COXA_DOE, COXA_UDOE, COXA_filter = CoxProcedureA(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance'], noFilter=False)
            results[elaboration_method_dictionary[elaboration_method]]['Cox procedure A'] = (COXA, uCOXA)
            information[(elaboration_method_dictionary[elaboration_method], 'Cox procedure A')] = {'chi-square test passed':bool(COXA_chitest), 'DoE':COXA_DOE, 'U(DoE)':COXA_UDOE, 'filter':COXA_filter}

        #CoxProcedureB
        if method.lower() in ('cox procedure b', 'coxb', 'all'):
            COXB, uCOXB, COXB_DOE, COXB_UDOE = CoxProcedureB(data_subset[f'{elaboration_method}_HL'], data_subset[f'{elaboration_method}_combined_variance']**0.5)
            results[elaboration_method_dictionary[elaboration_method]]['Cox procedure B'] = (COXB, uCOXB)
            information[(elaboration_method_dictionary[elaboration_method], 'Cox procedure B')] = {'DoE':COXB_DOE, 'U(DoE)':COXB_UDOE}

        #and also further statistical tests?

        #always show the result plot if not in iterative mode
        if not iterative:
            visualization.plot_results(fitted_data, f'{elaboration_method}_HL', f'{elaboration_method}_combined_variance', averages=results[elaboration_method_dictionary[elaboration_method]], title=elaboration_method_dictionary[elaboration_method])

    contributions = fitted_data.copy(deep=True)
    for elaboration_method in list_of_elaborated_data:
        contributions.drop(labels=f'{elaboration_method}_HL', axis=1, index=None, columns=None, level=None, inplace=True, errors='raise')
        contributions[f'{elaboration_method}_uHL'] = contributions[f'{elaboration_method}_uHL']**2 / contributions[f'{elaboration_method}_combined_variance']
        contributions[f'{elaboration_method}_uHF'] = contributions[f'{elaboration_method}_uHF']**2 / contributions[f'{elaboration_method}_combined_variance']
        contributions[f'{elaboration_method}_uMF'] = contributions[f'{elaboration_method}_uMF']**2 / contributions[f'{elaboration_method}_combined_variance']
        contributions[f'{elaboration_method}_uLF'] = contributions[f'{elaboration_method}_uLF']**2 / contributions[f'{elaboration_method}_combined_variance']
    information[('uncertainty', 'contribution')] = contributions
    information[('uncertainty', 'contribution table labels')] = list_of_elaborated_data
    
    return results, information

def iterative_procedure(dfr, MC_trials=10000, fit='all', method='all', max_iterations=10, control_check=0.01):
    """Iterative procedure for half-life determination
    useful for shorter half-lives where counting_time << half-life doesn't hold and decay correction during counting has to be considered
    
    Parameters
    ----------
    dfr : pandas.DataFrame
        dataframe containing data to fit
    MC_trials : int
        number of montecarlo trials (default 10000)
    fit : str
        fitting procedure of choice (default 'all')
    method : str
        averaging method of choice (default 'all')
    max_iterations : int
        maximum number of iterations before returning the result even if no convergence is achieved (default 10)
    control_check : float
        ratio of measurement uncertainty used as check to stop the iteration (default 0.01)
    
    Return
    ------
    fitted_data : pandas.DataFrame
        elaborated dataset
    half_life_results : pandas.DataFrame
        half_life value
    fitting_information : dict
        dictionary with additional information concerning fitting procedures, including iterations
    averaging_information : dict
        dictionary with additional information concerning the averaging procedures
    """
    #preparation stuff
    if fit not in ('weighted linear', 'wl', 'linear', 'l', 'weighted exponential', 'wexp', 'exponential', 'exp',
    'montecarlo linear', 'mcl', 'montecarlo exponential', 'mcexp'):
        fit = 'montecarlo linear'
    if method not in ('arithmetic average', 'aa', 'weighted average', 'wa', 'power-moderated mean', 'pmm', 'genetic algorithm', 
    'csg', 'ga', 'birge adjustment', 'birge', 'ba', 'dersimonianlaird', 'dsl', 'cox procedure a', 'coxa', 'cox procedure b', 'coxb'):
        method = 'ga'
    if MC_trials < 25000:
        MC_trials = 25000
    
    n_iter = 0
    iteration_information = {}
    previous_HL, previous_uHL = None, None
    while True:
        print(f'- iteration {n_iter + 1}')
        if previous_HL is not None and previous_uHL is not None:
            dfr = renormalize_data(dfr, previous_HL, previous_uHL)
        fitted_data, fitting_information = fit_data(dfr, autoplot=False, MC_trials=MC_trials, fit=fit)
        half_life_results, averaging_information = get_result(fitted_data, method=method, iterative=True)

        k0 = list(half_life_results.keys())[0]
        k1 = list(half_life_results[k0].keys())[0]
        _HL, _uHL = half_life_results[k0][k1][0], half_life_results[k0][k1][1]
        #all checks
        if previous_HL is not None:
            check = _HL / previous_HL - 1
        else:
            check = '-'
        iteration_information[('iteration', f'{n_iter + 1}')] = (_HL, _uHL, check)
        previous_HL, previous_uHL = _HL, _uHL
        if check != '-' and np.abs(check) < np.sqrt(2) * _uHL/_HL * control_check:
            break
        if n_iter >= max_iterations:
            break
        n_iter += 1
    
    ids = {'weighted linear':'wl', 'linear':'l', 'weighted exponential':'wexp', 'exponential':'exp', 'montecarlo linear':'mcl', 'montecarlo exponential':'mcexp'}
    reverse_ids = {'l':'linear fit', 'wl':'weigthed linear fit', 'wexp':'weighted exponential fit', 'exp':'exponential fit', 'mcl':'montecarlo linear fit', 'mcexp':'montecarlo exponential fit'}
    visualization.plot_results(fitted_data, f'{ids.get(fit, fit)}_HL', f'{ids.get(fit, fit)}_combined_variance', averages=half_life_results[reverse_ids[ids.get(fit, fit)]], title=reverse_ids[ids.get(fit, fit)])
    fitting_information = {**fitting_information, **iteration_information}
        
    return fitted_data, half_life_results, fitting_information, averaging_information

def elaboration(path, apt=False, nuclide=None, write_csv=False, MC_trials=10000, fit='all', method='all', output_path='', iterative=False):
    """Comprehensive function returning a dictionary with half-life values (with uncertainties) from name of the folder where the data are found
    
    Parameters
    ----------
    path : str (Path)
        directory name where measurement files are found
    apt : bool
        show plots during elaboration (default False)
    write_csv : bool
        write csv files at various stages of the elaboration (default False)
    MC_trials : int
        number of montecarlo trials (default 10000)
    fit : str
        fitting procedure of choice (default 'all')
    method : str
        averaging method of choice (default 'all')
    output_path : str (or Path)
        destination to save results files, if invalid defaults to {path}/elaboration (default '')
    iterative : bool
        whether to perform an iterative elaboration; useful for shorter half-lives to correct for counting decay (default False)
        In case this argument is True, 'fit' and 'method' arguments need to be selected and cannot be set to all
    
    Return
    ------
    half_life_results : dict
        half_life value
    information : dict
        dictionary with useful information
    """
    #get data in pandas.DataFrames format
    print('\ngathering information from files...')
    df, dfr, data_information = get_HL_data_from_dir(path, nuclide=nuclide, autoplot=apt)
    print('...done!\n\nfitting on data...')

    if iterative:
        fitted_data, half_life_results, fitting_information, averaging_information = iterative_procedure(dfr, MC_trials=MC_trials, fit=fit, method=method)
            
    else:
        fitted_data, fitting_information = fit_data(dfr, autoplot=apt, MC_trials=MC_trials, fit=fit)
        print('...done!\n\nperforming averages of results...')
    
        half_life_results, averaging_information = get_result(fitted_data, method=method)
        print('...done!')
        
    information = {**data_information, **fitting_information, **averaging_information}

    if write_csv:
        if not os.path.exists(output_path):
            output_path = os.path.join(path, 'elaboration')
            if not os.path.exists(output_path):
                os.makedirs(output_path)
        #csv output of the datasets
        date = datetime.datetime.today().strftime("%Y-%m-%d %H_%M_%S")
        df.to_csv(os.path.join(output_path, f'original_dataset_{date}.csv'))
        dfr.to_csv(os.path.join(output_path, f'elaborated_dataset_{date}.csv'))
        fitted_data.to_csv(os.path.join(output_path, f'fitted_dataset_{date}.csv'))
        
        r_df = pd.DataFrame(half_life_results)
        for column_name in r_df.columns:
            r_df[f'{column_name}_value / d'], r_df[f'{column_name}_uncertainty / d'] = zip(*r_df[column_name])
            r_df.drop(labels=column_name, axis=1, index=None, columns=None, level=None, inplace=True, errors='raise')
        r_df.to_csv(os.path.join(output_path, f'results_{date}.csv'))
        
        with open(os.path.join(output_path, f'additional_information_{date}.dict'), 'wb') as info_file:
            pickle.dump(information, info_file)

    print('COMPLETED!')
    
    return half_life_results, information
    
def read_info(info_file):
    """Load and return the information dictionary from pickled file
    
    Parameters
    ----------
    info_file : str (Path)
        name or path pointing to the .dict file
    
    Return
    ------
    information : dict
        dictionary with detailed information about the whole procedure
    """
    with open(info_file, 'rb') as read_info_file:
        return pickle.load(read_info_file)

def load_config(config_file):
    """Load and return the configuration
    
    Parameters
    ----------
    config_file : str (Path)
        file name or path where settings are found
    
    Return
    ------
    settings : dict
        settings dictionary
        
    
    Configuration file should be structured liked reported below, all settings have to be listed after the [Elaboration] section:

    [Elaboration]
    autoplot = True
    nuclide = None
    write_csv = True
    MC_trials = 20000
    fit = all
    method = all
    output_path = 
    iterative = False
    
    
    """
    if not os.path.exists(config_file):
        #create_default_config()
        print('no configuration file found, fallback to default configuration')
        return {}
        
    settings = {}

    config = configparser.ConfigParser()
    config.read(config_file)
    
    try:
        settings['apt'] = config.getboolean('Elaboration', 'autoplot')
    except Exception:
        pass
    try:
        nuclide_safety = config.get('Elaboration', 'nuclide')
        if nuclide_safety.lower() != 'none':
            settings['nuclide'] = nuclide_safety
    except Exception:
        pass
    try:
        settings['write_csv'] = config.getboolean('Elaboration', 'write_csv')
    except Exception:
        pass
    try:
        MC_trials = config.getint('Elaboration', 'MC_trials')
        if MC_trials > 100 and MC_trials < 10000000:
            settings['MC_trials'] = MC_trials
    except Exception:
        pass
    try:
        settings['fit'] = config.get('Elaboration', 'fit')
    except Exception:
        pass
    try:
        settings['method'] = config.get('Elaboration', 'method')
    except Exception:
        pass
    try:
        settings['output_path'] = config.get('Elaboration', 'output_path')
    except Exception:
        pass
    try:
        settings['iterative'] = config.getboolean('Elaboration', 'iterative')
    except Exception:
        pass
    
    return settings