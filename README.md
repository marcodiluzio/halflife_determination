# halflife_determination

## General
halflife_determination is a package for evaluation of half-life value (and corresponding uncertainty) of radionuclides which activity is monitored by collecting information on beta decay through a series of TDCR (Triple-to-Double Coincidence Ratio) liquid scintillation measurements. It is written in python 3 and relies on scientific modules such as *numpy*, *pandas*, *scipy*, *consensusGen*.  
At present, the only input file format allowed is the csv resulting from TDCR elaboration performed with package [TDCRPy](https://github.com/RomainCoulon/TDCRPy). All input files to be considered for elaboration need to be inside a main folder; they can, however, be organized into subfolders.  
Additional formats and data from other measurement techniques will be investigated in the future.

## Structure
The package is composed of the following modules:
- *hl_elaboration* (measurement elaboration)
- *visualization*  (graphical visualization of the results)

*hl_elaboration* contains functions following a modular structure as it is composed of three main functions (*get_HL_data_from_dir*, *fit_data*, *get_result*), each one returning the input that feeds the following function until the final output is obtained. Initial input is a folder name containing the measurement results. After each function call, a wealth of output is returned, including *pandas.DataFrames* containing partial results along with dictionaries providing useful supplementary information. Finally, a single value of half-life is reported together with the corresponding estimated uncertainty, evaluated through sum of variances from independent contributions or MonteCarlo simulation.  
A convenience function (called *elaboration*) grouping the three main functions in a single call is available. It can also be accessed from the terminal command line and manages keyword arguments passed to the intermediate functions via an optional configuration file.  

## Measurement data management
Liquid scintillation measurements performed with TDCR method need to be taken over time for the investigated radionuclide. A minimum of 4 measurements for dataset are necessary, where a dataset indicates all measurements acquired for a sample in the same experimental conditions (extended dead time, coincidence window, scintillation cocktail). Additional measurements to evaluate Birks constant are also required. A so-called Birks evaluation consists in a tight series of measurements acquired for each dataset where the sample is measured with and wihtout various ND filters to artificially variate the TDCR value and record the changes to the evaluated activity. at the beginning (or at regular intervals) of the measurement campaign is also supported.  
All files have to be of .csv format returned from '' and saved into a main folder (which will be recalled by the software). No special filename structure is required and files can be organized in subfolders (down to 2 levels) within the main folder and still be recognized by the code. Showed below it's an example of a valid way to manage your measurement files.  

```
Pm-147/
├── 2025-06/
│   ├── TDCR_Results_Pm-147_LNHB_2025_PS1_0.1_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_PS1_0.1_nanoTDCR_10_100_12.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_PS1_0.1_nanoTDCR_10_100_5.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_UG4_0_nanoTDCR_50_50_8.csv
├── 2025-08/
│   ├── TDCR_Results_Pm-147_LNHB_2025_August_PS1_0_nanoTDCR_10_100_10.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_August_PS1_0_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_August_PS1_0_nanoTDCR_10_100_12.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_August_UG4_0_nanoTDCR_50_50_9.csv
├── 2025-09/
│   ├── TDCR_Results_Pm-147_LNHB_2025_September_PS1_0_nanoTDCR_10_100_10.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_September_PS1_0_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_September_PS1_0_nanoTDCR_10_100_12.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_September_UG4_0_nanoTDCR_50_50_9.csv
├── 2025-10/
│   ├── TDCR_Results_Pm-147_LNHB_2025_October_PS1_0_nanoTDCR_10_100_10.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_October_PS1_0_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_October_PS1_0_nanoTDCR_10_100_12.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_October_UG4_0_nanoTDCR_50_50_9.csv
├── 2025-11/
│   ├── TDCR_Results_Pm-147_LNHB_2025_November_PS1_0_nanoTDCR_10_100_10.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_November_PS1_0_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_November_PS1_0_nanoTDCR_10_100_12.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_November_UG4_0_nanoTDCR_50_50_9.csv
├── 2025-12/
│   ├── TDCR_Results_Pm-147_LNHB_2025_December_PS1_0_nanoTDCR_10_100_10.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_December_PS1_0_nanoTDCR_10_100_11.csv
│   ├── TDCR_Results_Pm-147_LNHB_2025_December_PS1_0_nanoTDCR_10_100_12.csv
│   ├── ...
│   └── TDCR_Results_Pm-147_LNHB_2025_December_UG4_0_nanoTDCR_50_50_9.csv
```

## Quick start
### Install
- install the package from pip
  ```shell
  pip install halflife_determination
  ``` 

### Run from script (option 1)
- import modules in python script
  ```python
  from halflife_determination import hl_elaboration as hle   #main script with all relevant functions  
  from halflife_determination import visualization           #optional, only for custom visualization
  ```

- call elaboration function
  ```python
  kwargs = {'apt':False, 'nuclide':None, 'write_csv':True, 'MC_trials':20000, 'fit':'all', 'method':'all', 'output_path':'', 'iterative':False}  
  results, information = hle.elaboration('path_to_data_folder', **kwargs)
  ```

- or, for a finer control, call main functions directly
  ```python
  original_dataset, elaborated_dataset, information = hle.get_HL_data_from_dir('path_to_data_folder', nuclide=None, autoplot=False)  
  fitted_data, information = hle.fit_data(elaborated_dataset, data_threshold=3, MC_trials=10000, autoplot=False, fit='all')  
  results, information = hle.get_result(fitted_data, method='all', iterative=False)
  ```

### Run from terminal (option 2)
- run directly from the command line (terminal or PowerShell)
  ```shell
  python -m halflife_determination 'path_to_data_folder' 'path_to_configuration_file'{optional}  
  ```
  (which is equivalent of calling elaboration)
  ```python
  hl_elaboration.elaboration('path_to_data_folder', **kwargs)
  ```
  if the 'path_to_configuration_file' points to a file that does not exist, a default configuration file will be created at that destination  


## Contacts
Marco Di Luzio  
m.diluzio@inrim.it
