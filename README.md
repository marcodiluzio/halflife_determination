# halflife_determination

## General
halflife_determination is a package for evaluation of half-life value (and uncertainty) of radionuclides which activity is monitored collecting information on beta decay through a series of TDCR (Triple-to-Double Coincidence Ratio) liquid scintillation measurements. It is written in python version 3.13 and relies on scientific modules such as *numpy*, *pandas*, *scipy*, *consensusGen*.

## Structure
The package is composed of the following modules:
- *hl_elaboration* (measurement elaboration)
- *visualization*  (graphical visualization of the results)

*hl_elaboration* contains functions following a modular structure as it is composed of three main functions (*get_HL_data_from_dir*, *fit_data*, *get_result*), each one returning the input that feeds the following function until the final output is obtained. Initial input is a folder name containing the measurement results. After each function call, a wealth of output is returned, including *pandas.DataFrames* containing partial results along with dictionaries providing useful supplementary information. Finally, a single value of half-life is reported together with the corresponding estimated uncertainty, evaluated through sum of variances from independent contributions or MonteCarlo simulation.  
A convenience function (called *elaboration*) grouping the three main functions in a single call is available. It can also be accessed from the terminal command line and manages keyword arguments passed to the intermediate functions via an optional configuration file.

## Quick start
### Install
- install the package from pip  
  `pip install halflife_determination`  

### Run from script (option 1)
- import modules in python script  
  `from halflife_determination import hl_elaboration as hle   #main script with all relevant functions`  
  `from halflife_determination import visualization           #optional, only for custom visualization`  

- call elaboration function  
  `kwargs = {'apt':False, 'nuclide':None, 'write_csv':True, 'MC_trials'=20000, 'fit':'all', 'method':'all', 'output_path':'', 'iterative':False}`  
  `results, information = hle.elaboration('path_to_data_folder', **kwargs)`

- for a finer control, call main function directly
  A  
  B  
  C  

### Run from terminal (option 2)
- run directly from the command line (terminal or PowerShell)  
  `python -m halflife_determination 'path_to_data_folder' 'path_to_configuration_file'{optional}` (which is equivalent of calling `hl_elaboration.elaboration('path_to_data_folder', **kwargs)`)


## Contacts
Marco Di Luzio  
m.diluzio@inrim.it

\[GitHub-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
