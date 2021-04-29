#####################################################################################################################################################################################################
################################################################################################# READ ME ###########################################################################################
#####################################################################################################################################################################################################

Thank you for reading our paper 'A Distribution of Human Attention to Moments in Time'. The following elements will allow you to replicate the findings of our analysis.

The main "code" folder contains the following elements:


(A) A folder containing a snippet of the R code used to download and clean the dataset:
==========================================================================================
1) The first sub-folder "extract_code" contains snippets of the R-script used to download the search index data from Google Trends ("extraction_sample_code.R"). It illustrates how the raw input processed later in the analysis was obtained. In order to run the code, the reader should change the directory path specified in Line 20. It is currently set to "/home/ltmb/Dropbox/Sebastien-Leo-Sol/time writing/code/r_code/extract_code/" and should be changed to the equivalent path on the reader's machine. All subsequent paths are relative.  
The particular script "extraction_sample_code.R" downloads worldwide search data (monthly level from 2008-01-01 to 2018-12-31) for the topics corresponding to the following keywords:
"2015", "2016", "april", "september", "valentine's day", "eid al-fitr", "sunscreen", "donald trump", "cat". 
Running the code creates a file named "placebo_comparison.txt", which gets stored inside the "/extract_code/" folder. This file is identical to the one pre-contained inside the folder "/figures_code/". The latter one is used to generate Extended Figure 1 (see below). 

"extraction_sample_code.R" also creates a file named "placebo_comparison_new.txt". It gets stored inside the "code/r_code/extract_code" folder. This file is identical to the one pre-contained inside the folder "code/r_code/figures_code/". The latter is used to generate Figure 2 in the main text (see below).


(B) A folder dedicated to the R code and STATA code written for the analysis and generate the figures in the paper:
===================================================================================================================

1) Analysis
Run in this specific order the following scripts:
	1.1  deconvo_day.do (analysis for daily queries)
	1.2  deconvo_month.do (analysis for monthly queries)
	1.3  deconvo_year (analysis for yearly queries)
	1.4  fitting_smoothing_XXX.R (betas fitting and smoothing) NOTE: replace XXX by either "day", "month", "year" or "year_world"
	1.4b bootstrap_fitting.R (same as 1.4 but for bootstrap)	
	1.5  format_predictor.do (clean and reshape smoothed betas to be used in predictions)
	1.6  predictions.do (predict and generate goodness of fit measures) NOTE: calls script_r2.do 


=========================
How to reproduce figures:
=========================

## Main text

# Figure1
	N/A
# Figure2
	This figure relies on the "placebo_comparison_new.txt" file contained inside the folder "code/r_code/figures_code/". To generate this file from scratch, refer to (A) from the README.txt
	fig2_master.do NOTE: calls the following scripts "fig2_fit_day.do", "fig2_fit_month.do", "fig2_fit_year.do"	
# Figure3
	(3a) run steps 1.3 and 1.4
	(3b) run "generating_figures_code.R" contained inside "code/r_code/figures_code/". 
	NOTE 1: The reader needs to change the directory path in Line 11. It is currently set to "/home/ltmb/Dropbox/Sebastien-Leo-Sol/time writing/code/r_code/figures_code/". 
	All subsequent paths are relative.
	NOTE 2: It relies on "KAT_share_CI.csv", contained in "code/r_code/figures_code/". This file contains all the country-level shares of attention towards the past, future, and present. 
	        This file can be generated using the steps described in 1). 
	(3c) run violin_plot.do (requires steps 1.1, 1.2, and 1.3)
	(3d) run step 1.4 for "day", "month" and "year" (requires steps 1.1, 1.2, and 1.3)
# Figure4
	(4a) run "generating_figures_code.R" contained inside "code/r_code/figures_code/". Related code is contained within the lines 56-122 of "generating_figures_code.R". 
	NOTE: It relies on "KAT_share_CI.csv", contained in "code/r_code/figures_code/". It also relies on the folder "ne_10m_admin_0_countries" 
	(source: Natural Earth, url: https://www.naturalearthdata.com/downloads/10m-physical-vectors/) containing the country boundaries for the maps.These shapefile data are in the public domain. 
	(4b) See how 4a) is generated. Related code is contained within the lines 124-172 of "generating_figures_code.R". 
	(4c) run fig4.do (requires steps 1.1, 1.2, and 1.3)

## Extended figures
# Figure1
	run "generating_figures_code.R" contained inside "code/r_code/figures_code/". Related code is contained within the lines 174-194 of the script. 
	NOTE: It relies on the file "placebo_comparison.txt", contained inside "code/r_code/figures_code/". To generate this file from scratch, refer to (A) from the README.txt.
# Figure2
	N/A
# Figure3
	run ED_fig3.do (requires 1.3 and 1.4b)
# Figure4
	run violin_plot.do (requires steps 1.1, 1.2, and 1.3)
# Figure 5
	run step 1.6
# Figure 6 
	run "generating_figures_code.R" contained inside "code/r_code/figures_code/". Related code is contained within the lines 197-260 of the script. 
# Figure 7
	run plot_KAT_robust.R then panel 2x2_robust.do

