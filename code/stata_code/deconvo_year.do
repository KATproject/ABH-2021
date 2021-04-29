/*
Author: S. Annan-Phan
Date: 8/4/2020
This do file run the deconvolution model for "year" related querries
*/


//Set path to working directory
cap cd "~\Dropbox\temporal_awareness"
cap cd "/home/sphan/Research_projects/temporal_awareness"

set seed 1234 // set seed for bootstrap procedure (set to 1234 in paper)
set maxvar 10000
global BS = 0 // set # of bootstrap iteration (set to 3000 in paper)
local L = floor(250 / 7) // set number of lags/leads (250 days in paper)
set matsize 11000

//File to store results by region
tempname memhold
postfile `memhold' str18 region tune past_b past_se ///
pres_b pres_se fut_b fut_se ///
using "saves/summary_attention_results.dta", replace 
***********************************
*** data loading and formatting ***
***********************************

import delim "weekly_interest_years.txt", clear varn(1)
replace date = substr(date,1,10)
g date2 = date(date,"YMD",2000)
drop date
rename date2 date
format date %td
g week = week(date)
g year = year(date)
g month = month(date)
cap replace hits = "0" if hits == "<1"
cap destring hits, replace force
keep if year >= 2008 & key > 2008
sort geo keyw date
drop if geo == "KP"
egen t = seq(), by(keyw geo) 

g D = (year == keyw)  // gen dirac function (indicator)

egen seq_targ = seq(), by(geo keyw D)

replace D = 1  if year[_n+1] == keyw[_n+1] & seq_targ[_n+1] == 1 

drop seq_targ 
sort geo keyw date
*id
egen id = group(keyw geo)
*panel
xtset id t
drop if date == .
sort geo keyw date
egen country = group(geo)

// re-scale to adjust for outlier
tempvar med maxmed
egen `med' = median(hits), by(geo key)
egen `maxmed' = max(`med'), by(geo)	
g raw = hits
replace hits = hits * `maxmed' / `med'

// gen lags/leads of dirac function
foreach k of num 1(1)`L'{
	g L`k' = D[_n-`k'] 
	g F`k' = D[_n+`k'] 
}
******************************************************************************//



// list of 20 countries with daily/monthly querries 
// For fig 3 (c)
local subsamp = "AR AU BR CL CO DE ES FR GB HK ID IN MX NL NO NZ PE SE SG US" 

g subsamp = 0
foreach c in `subsamp' {
 replace subsamp = 1 if geo == "`c'"
}


******************
*    Analysis    *
******************
save "saves/data_year_reg", replace	

* Pooled reg
		
									
					
// code for matrix							
local LC_past = "_b[L1]"
local LC_fut = "_b[F1]"
foreach q of num 2/34 {
	local LC_past = "`LC_past' + _b[L`q']"
	local LC_fut = "`LC_fut' + _b[F`q']"
}


if $BS != 0 {
	local matL = `L' * 2 +1
	forvalues i = 1/$BS {
		di "itearation `i'"
		mat KAT`i' = J(`matL',3,0)
		preserve
			bsample, cluster(id)
			qui rreg hits L* F* D `FE', tune(6)
		restore
		// save matrix
		foreach k of num `L'(1)1{
			local posi = `L' - `k' + 1
			mat KAT`i'[`posi',1] = _b[L`k']
			mat KAT`i'[`posi',2] = _se[L`k']			
			mat KAT`i'[`posi',3] = `i'		
		}	
		local mid = `L' + 1
		mat KAT`i'[`mid',1] =  _b[D]
		mat KAT`i'[`mid',2] = _se[D]					
		mat KAT`i'[`mid',3] =  `i'
		
		
		foreach k of num 1(1)`L'{
			local pos = `k' + `mid'
			mat KAT`i'[`pos',1] = _b[F`k']
			mat KAT`i'[`pos',2] = _se[F`k']
			mat KAT`i'[`pos',3] = `i'	
		}		
	
	}	

	preserve
		forvalues k = 1/$BS {
		clear 
		svmat KAT`k'
		rename (KAT`k'1 KAT`k'2 KAT`k'3) (beta se ite)
		egen L = seq(), from(-35)
		tempfile f`k'
		save `f`k''
		}
	restore	
	
	preserve
	use `f1', clear
	forvalues k = 2/$BS {
		qui append using `f`k''
	}
	save "saves/BS_kat_1234.dta", replace
	restore
	
}

rreg hits L* F* D, tune(6) 

// save post estimates statistics with regular SE (no bootstrap)

nlcom (`LC_past') / (`LC_past' + `LC_fut' + _b[D]),post
local b_past = _b[_nl_1]
local se_past = _se[_nl_1]

qui rreg hits L* F* D , tune(6)
nlcom (`LC_fut') / (`LC_past' + `LC_fut' + _b[D]),post
local b_fut = _b[_nl_1]
local se_fut = _se[_nl_1]

qui rreg hits L* F* D , tune(6)	
nlcom (_b[D]) / (`LC_past' + `LC_fut' + _b[D]),post
local b_pres = _b[_nl_1]
local se_pres = _se[_nl_1]	

post `memhold' ("full'") (6) (`b_past') (`se_past') (`b_pres') (`se_pres') ///
 (`b_fut') (`se_fut')

rreg hits L* F* D, tune(6)
// matrix to store full set of betas (raw KAT)
local L = 35
foreach SUBSAMPLE of num 0 1 { // repeat for (0) full sample and (1) sub sample of 20 countries
	if `SUBSAMPLE' == 1 {
		rreg hits L* F* D if subsamp == 1, tune(6)
	}
	else{
		rreg hits L* F* D, tune(6)
	}
	mat lags = J(`L',4,0)

	foreach k of num 1(1)`L'{
		mat lags[`k',1] = _b[L`k']
		mat lags[`k',2] = _se[L`k']
		mat lags[`k',3] = _b[F`k']
		mat lags[`k',4] = _se[F`k']
	}
	// graph betas
	preserve
		clear 
		svmat lags
		rename (lags1 lags2 lags3 lags4) (beta0 se0 beta1 se1)
		egen L = seq()
		reshape long beta se, i(L) j(type)
		replace L = L * -1 if type == 0
		drop type
		sort L
		local n_obs = _N + 1
		set obs `n_obs' 	
		replace L = 0 if L == .
		replace beta = _b[D] if L == 0
		replace se = _se[D] if L == 0		
		replace L = L * 7 
		g min = beta - 1.96 * se
		g max = beta + 1.96 * se	
		g state = "0"
		tempfile pooled 
		save `pooled'
		
		sort L
		g w = 1 / se

		tw rcap max min L, col(dkgreen*.50) ///
		|| scatter beta L , mc(dkgreen) xline(0, lc(black) lp(dash) lw(thin)) ///
		legend(off) yline(0, lc(black) lw(thin)) xtitle(days) ytitle(betas) ///
		name(main`SUBSAMPLE', replace) title(Weekly)
		
		if `SUBSAMPLE' == 1 {
			outsheet L beta se min max w using pooled_inter_year.csv , comma replace
		}
		else{
			outsheet L beta se min max w using pooled_inter_year_world.csv , comma replace
		}
	restore
}

// Fixed effect Robustness-test

local L = 35
local FE_count = 0
foreach FE in "" "i.country" "i.country i.year" "i.country i.year i.keyw" "i.country#i.year" { 
	rreg hits L* F* D `FE', tune(6)
	mat lags = J(`L',4,0)

	foreach k of num 1(1)`L'{
		mat lags[`k',1] = _b[L`k']
		mat lags[`k',2] = _se[L`k']
		mat lags[`k',3] = _b[F`k']
		mat lags[`k',4] = _se[F`k']
	}
	// graph betas
	preserve
		clear 
		svmat lags
		rename (lags1 lags2 lags3 lags4) (beta0 se0 beta1 se1)
		egen L = seq()
		reshape long beta se, i(L) j(type)
		replace L = L * -1 if type == 0
		drop type
		sort L
		local n_obs = _N + 1
		set obs `n_obs' 	
		replace L = 0 if L == .
		replace beta = _b[D] if L == 0
		replace se = _se[D] if L == 0		
		sort L
		g w = 1 / se
		g FE = `FE_count'
		tempfile f`FE_count'
		save `f`FE_count''
		local FE_count = `FE_count' + 1
	restore
}

preserve
	use `f0'
	append using `f1'
	append using `f2'
	append using `f3'
	append using `f4'
	outsheet L FE beta se w using pooled_inter_year_world_robustness.csv , comma replace
restore

// Spatial heterogeneity (fig 4 a-b) ************************************************

levelsof(geo), local(state_list)

foreach s in `state_list' {
	di "`s'"		
	rreg hits L* F* D  if geo == "`s'", tune(6)
	nlcom (`LC_past') / (`LC_past' + `LC_fut' + _b[D]),post
	local b_past = _b[_nl_1]
	local se_past = _se[_nl_1]

	rreg hits L* F* D  if geo == "`s'", tune(6)
	nlcom (`LC_fut') / (`LC_past' + `LC_fut' + _b[D]),post
	local b_fut = _b[_nl_1]
	local se_fut = _se[_nl_1]

	rreg hits L* F* D  if geo == "`s'", tune(6)
	nlcom (_b[D]) / (`LC_past' + `LC_fut' + _b[D]),post
	local b_pres = _b[_nl_1]
	local se_pres = _se[_nl_1]	
	post `memhold' ("`s'") (6) (`b_past') (`se_past') (`b_pres') (`se_pres') ///
	 (`b_fut') (`se_fut')
}
postclose `memhold'

preserve
	keep if subsamp == 1
	levelsof(geo), local(state_sub_list)
restore
foreach s in `state_sub_list' {
	rreg hits L* F* D  if geo == "`s'", tune(6)
	mat lags = J(`L',4,0)

	foreach k of num 1(1)`L'{
		mat lags[`k',1] = _b[L`k']
		mat lags[`k',2] = _se[L`k']
		mat lags[`k',3] = _b[F`k']
		mat lags[`k',4] = _se[F`k']
	}
	preserve
		clear 
		svmat lags
		rename (lags1 lags2 lags3 lags4) (beta0 se0 beta1 se1)
		egen L = seq()
		reshape long beta se, i(L) j(type)
		replace L = L * -1 if type == 0
		drop type
		sort L
		local n_obs = _N + 1
		set obs `n_obs' 	
		replace L = 0 if L == .
		replace beta = _b[D] if L == 0
		replace se = _se[D] if L == 0		
		replace L = L * 7 
		g min = beta - 1.96 * se
		g max = beta + 1.96 * se	
		g state = "`s'"
		tempfile beta`s'
		save `beta`s''
	restore	
}
* Save full matrix per regions (fig3)
preserve
	use `pooled', clear
	foreach s in `state_sub_list' {
		append using `beta`s''
	}
	save "saves/state_beta_year", replace	
restore

// Time heterogeneity (fig 4 c) ************************************************
local L = 35
foreach Y of num 2009/2018 {
	preserve	
	reg hits L* F* D if subsamp == 1 & keyw == `Y'

	mat lags = J(`L',4,0)

	foreach k of num 1(1)`L'{
		mat lags[`k',1] = _b[L`k']
		mat lags[`k',2] = _se[L`k']
		mat lags[`k',3] = _b[F`k']
		mat lags[`k',4] = _se[F`k']
	}
	*

	clear 
	svmat lags
	rename (lags1 lags2 lags3 lags4) (beta0 se0 beta1 se1)
	egen L = seq()
	reshape long beta se, i(L) j(type)
	replace L = L * -1 if type == 0
	drop type
	sort L
	local n_obs = _N + 1
	set obs `n_obs' 	
	replace L = 0 if L == .

	replace beta = _b[D] if L == 0
	replace se = _se[D] if L == 0		

	replace L = L * 7 
	
	*1 past , 2 pres, 3 fut
	g grp = L < 0
	replace grp = 2 if L == 0
	replace grp = 3 if L >0
	drop if abs(L)>=240
	
	collapse (sum) beta, by(grp)
	egen totint = total(beta)
	replace beta = beta / totint
	drop totint
	g year = `Y'
	****
	tempfile share`Y'
	save `share`Y''
	restore
}
preserve
	use `share2009', clear
	foreach j of num 2010/2018 {
		append using `share`j''
	}
	g tag = "year"
	tw connect beta year if grp == 3 ///
	|| connect beta year if grp == 2 ///
	|| connect beta year if grp == 1, legend(lab(1 "fut") lab(2 "pres") lab(3 "past"))
	save "saves/share_over_time_year.dta", replace
restore

* Create csv for ternary plot (fig3.b) and maps (fig4)
preserve
	use "saves/summary_attention_results.dta", clear
	keep if tune == 6 
	drop tune
	foreach period in "past" "pres" "fut" {
		g min_`period' = `period'_b - 1.96 * `period'_se
		g max_`period' = `period'_b + 1.96 * `period'_se
		drop `period'_se
	}
	outsheet * using "saves/KAT_share_CI.csv", replace
restore


