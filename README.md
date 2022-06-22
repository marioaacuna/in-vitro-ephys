# In-vitro Ephys
This project includes scripts to analyze in vitro ephys data
# INTRO
This pipeline runs on MATLAB. it will help you analyse in vitro electrophysiological data obtainedin the Nevian Lab. This pipeline was created and is mantained by Mario A. Acuña. If you have questions, please ask him (mario.acuna@unibe.ch)
So far we can analyse:
# 
- SAG properties
- AP features
- Pairing analysis
- Evoked EPSPs
- Evoked EPSCs
- Paired-pulse ratios
- EPSP summation.
- In future releases:- Spontaneous EPSCs 
(the detailed information about each of them will come in a future release)

# How to set up your ephys analysis pipeline:
	- login to: https://gitlab.com/nevian_group_unibe/in-vitro-ephys
	- click on Clone -> copy https link
	- go to the folder where you have all re repos:
		- Mouse right click and select 'clone'
		- It should appear what you copied already as the default https link to clone
		- And as simple as that you will have the ephys analysis on your repo folder
		
# How to run Ephys analysis
	- open Matlab
	- Create a favorite script:
		- Home > favorite > new favorite
		- give it a name Label:'run_epyhs', or so
		Code:
		
```matlab 
% This is what you will have to write
% copy and paste this code
close all force hidden; 
clc;
clear;
cd('PATH/TO/YOUR/REPO') % To be changed (add in repository)
run_in_vitro_ephys_analysis()
```
#
	- This will run an app where you can select what to analyse
		- For AP analysis, a new app will open and you have to select the experimenter name
		- Then choose all the other parameters you see on the app
		- For all the other analysis, follow the same logic
		
# That's it! enjoy!
