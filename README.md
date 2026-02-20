# halflife_determination

## General
halflife_determination is a package for evaluation of half-life value (and corresponding uncertainty) of radionuclides which activity is monitored by collecting information on beta decay through a series of TDCR (Triple-to-Double Coincidence Ratio) liquid scintillation measurements. It is written in python 3 and relies on scientific modules such as *numpy*, *pandas*, *scipy*, *consensusGen*.  
At present, the only input file format allowed is the csv resulting from TDCR elaboration performed with package [TDCRPy](https://github.com/RomainCoulon/TDCRPy). All input files to be considered for elaboration need to be inside the same main folder; they can, however, be organized into subfolders.  
Additional formats and data from other measurement techniques will be investigated in the future.

## Structure
The package is composed of the following modules:
- *hl_elaboration* (measurement elaboration)
- *visualization*  (graphical visualization of the results)

*hl_elaboration* contains functions following a modular structure as it is composed of three main functions (*get_HL_data_from_dir*, *fit_data*, *get_result*), each one returning the input that feeds the following function until the final output is obtained. Initial input is a folder name containing the measurement results. After each function call, a wealth of output is returned, including *pandas.DataFrames* containing partial results along with dictionaries providing useful supplementary information. Finally, a single value of half-life is reported together with the corresponding estimated uncertainty, evaluated through sum of variances from independent contributions or MonteCarlo simulation.  
A convenience function (called *elaboration*) grouping the three main functions in a single call is available. It can also be accessed from the terminal command line and manages keyword arguments passed to the intermediate functions via an optional configuration file.  

## Input data management
Liquid scintillation measurements performed with TDCR method need to be taken over time for the investigated radionuclide. A minimum of 4 measurements for dataset are necessary, where a dataset indicates all measurements acquired for a sample in the same experimental conditions (extended dead time, coincidence window, scintillation cocktail). Additional measurements to evaluate Birks constant are also required. A so-called Birks evaluation consists in a tight series of measurements acquired for each dataset where the sample is measured with and wihtout various ND filters to artificially variate the TDCR value and record the changes to the evaluated activity. At least one Birks evaluation is required for each dataset, however multiple evaluations held at regular intervals during the measurement campaign are also supported. Showed below it's an example of an allowed way to structure the measurements where t0, ... tk are dates and Sample measurment(t) is the acquisition performed at the corresponding date  

```
t0  ── Sample measurement(t0) ────────────────
|   ── Sample measurement(t0) + ND filter 1 ──
|   ── Sample measurement(t0) + ND filter 2 ──
|   ── Sample measurement(t0) + ND filter n ──
|
t1  ── Sample measurement(t1) ────────────────
|
t2  ── Sample measurement(t2) ────────────────
|
...
|
tn  ── Sample measurement(tn) ────────────────
|
tj  ── Sample measurement(tj) ────────────────
|   ── Sample measurement(tj) + ND filter 1 ──
|   ── Sample measurement(tj) + ND filter 2 ──
|   ── Sample measurement(tj) + ND filter n ──
|
...
|
tk  ── Sample measurement(tk) ────────────────
```

All files (measurement for halflife determination and Birks evalaution) have to be of .csv format returned from '[REFERENCE!]' module and saved into a main folder (which will be recalled by the software).  
Measurement filenames have to be structured in a precise way reported below  

```
TDCR_Results_{nuclide}_{lab}_{year}_{month}_{cocktail and sample number}_{ND filter}_{instrument code}_{extended dead time}_{coincidence window}_{progressive integer}.csv
or
TDCR_Results_{nuclide}_{lab}_{year}_{cocktail and sample number}_{ND filter}_{instrument code}_{extended dead time}_{coincidence window}_{progressive integer}.csv

filename example: TDCR_Results_Pm-147_LNHB_2025_PS1_0_nanoTDCR_10_100_12.csv
```

with words in {} are relevant measurement information and should not contain underscores.
- The {cocktail and sample number} has to be an alphanumeric string as the numeric part suffix identifies the sample code.
- If {ND filter} is 0 means no ND filter is applied to the sample and the measurement is adopted for half-life determination, it defines a Birks evaluation measurement otherwise.  
Files can be organized in subfolders (down to 2 levels) within the main folder and still be recognized by the code. Showed below it's an example of a valid way to manage your measurement files:  

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

And the typical format of a csv input file:  

```
Section;Parameter;Value;Uncertainty;Unit;Description
General;Nuclide;Pm-147;;-;Radionuclide name
General;Reference Date;2024-12-16 12:00:00;;-;Reference date
General;Measurement Date;2025-11-18 15:05:00;;-;Date of measurement
General;Half-life (T1/2);82786000.0;13000.0;days;Half-life
Source;Laboratory;LNHB;;-;Provider of the standard solution
Source;LS Cocktail;ProSafe+;;-;Liquid scintillation cocktail
Source;Aqueous Fraction;0.0625;;rel;Aqueous fraction in scintillator
Source;Grey Filter ND;0;;-;Neutral density grey filter
...
Model Analytic;Combined Rel. Unc. Activity;0.0012135880474480828;;rel;Combined relative uncertainty for activity
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
