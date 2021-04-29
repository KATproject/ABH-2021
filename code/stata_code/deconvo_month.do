/*
Author: S. Annan-Phan
Date: 8/4/2020
This do file run the deconvolution model for "month" related querries
*/


//Set path to working directory
cap cd "~\Dropbox\temporal_awareness"
cap cd "/home/sphan/Research_projects/temporal_awareness"
set seed 1234 

*Number of lag/leads (weeks)
local L =  90
********************************************************************************
*data cleaning
********************************************************************************

import delim "international_daily_interest_months_noyear.txt", clear varn(1)
egen t = seq(), from(1) to(270)
egen id = seq(), block(270)
drop keyword
rename keyword_e keyword

*clean date
replace date = substr(date,1,10)
g date2 = date(date,"YMD",2000)
drop date
rename date2 date
format date %td
g week = week(date)
g year = year(date)
g month = month(date)
egen state = group(geo)

*Decompose the "search" into month and year
g keyword1 = keyword
g year_s = 0

*clean month translation
g month_search = 1 if keyword1 == "january" | keyword1 == "januar" | keyword1 == "janvier" | keyword1 == "enero"	
replace month_search = 2 if keyword1 == "february" | keyword1 == "febrero" | keyword1 == "februar" | keyword1 == "fÃ©vrier"
replace month_search = 3 if keyword1 == "march" | keyword1 == "mars" | keyword1 == "marzo" | keyword1 == "mÃ¤rz"
replace month_search = 4 if keyword1 == "april" | keyword1 == "abril" | keyword1 == "avril"
replace month_search = 5 if keyword1 == "may" | keyword1 == "mai" | keyword1 == "mayo"
replace month_search = 6 if keyword1 == "june" | keyword1 == "juni" | keyword1 == "junio" | keyword1 == "juin"
replace month_search = 7 if keyword1 == "july" | keyword1 == "julio" | keyword1 == "juli" | keyword1 == "juillet"
replace month_search = 8 if keyword1 == "august" | keyword1 == "agosto" | keyword1 == "aout"
replace month_search = 9 if keyword1 == "september" | keyword1 == "septembre" | keyword1 == "septiembre"
replace month_search = 10 if keyword1 == "october" | keyword1 == "octobre" | keyword1 == "octubre" | keyword1 == "oktober"
replace month_search = 11 if keyword1 == "november" | keyword1 == "novembre" | keyword1 == "noviembre"
replace month_search = 12 if keyword1 == "december" | keyword1 == "dezember" | keyword1 == "dÃ©cembre" | keyword1 == "diciembre"
drop keyword1

replace year_s = year
replace year_s = year + 1 if (month_s == 1 | month_s ==2 | month_s ==3 | month == 4) & month == 12 
replace year_s = year + 1 if (month_s == 1 | month_s ==2 | month_s ==3) & month == 11 
replace year_s = year + 1 if (month_s == 1 | month_s ==2) & month == 10 
replace year_s = year + 1 if (month_s == 1) & month == 9

replace year_s = year - 1 if (month_s == 12 | month_s ==11 | month_s ==10 | month == 9) & month == 1
replace year_s = year - 1 if (month_s == 12 | month_s ==11 | month_s ==10) & month == 2
replace year_s = year - 1 if (month_s == 12 | month_s ==11) & month == 3
replace year_s = year - 1 if (month_s == 12) & month == 4

*clean the score variable
replace hits = "0" if hits == "<1"
destring hits, replace force
save "international_daily_interest_months_noyear.dta", replace
keep if year >= 2008
sort geo keyw date
*time

********************************************************************************
* generate panel of events
********************************************************************************

sort geo keyw date
*id

*panel
xtset id t
*Dummy if t is in S
drop if date == .
sort geo keyw date
g dow = dow(date)
g dom = day(date)
g we = (dow == 0 | dow == 6)
egen lastday = max(dom), by(month)
g D = (month_s == month)

*Correct for higher amplitude first day of the month
merge m:1 geo using "amplitude_by_country_noyear_ite0", nogen
drop if geo == "PT" | geo == "PH" | geo == "IT"
g DSW = D
replace DSW = dsw if (D == 1 & dom == 1) 
drop D
rename DSW D

foreach k of num 1(1)`L'{
	g L`k' = D[_n-`k'] 
	g F`k' = D[_n+`k'] 
}

preserve
collapse hits, by(t)
save "avg_hits_month_inter.dta", replace
restore



******************
*    Analysis    *
******************

save "saves/data_month_reg", replace	

levelsof(geo), local(state_list)

*Matrix of estimates 

* Pooled reg
egen cluster_id = group(state year)
rreg hits L* F* D, tune(6)

*statistics for fig3

*compile the long command
local LC_past = "_b[L1]"
local LC_fut = "_b[F1]"
foreach q of num 2/34 {
	local LC_past = "`LC_past' + _b[L`q']"
	local LC_fut = "`LC_fut' + _b[F`q']"
}

*share of attention to past
nlcom (`LC_past') / (`LC_past' + `LC_fut' + _b[D]),post
local b_past = _b[_nl_1]
local se_past = _se[_nl_1]
*share of attention to futur
qui rreg hits L* F* D , tune(6)
nlcom (`LC_fut') / (`LC_past' + `LC_fut' + _b[D]),post
local b_fut = _b[_nl_1]
local se_fut = _se[_nl_1]
*share of attention to present
qui rreg hits L* F* D , tune(6)	
nlcom (_b[D]) / (`LC_past' + `LC_fut' + _b[D]),post
local b_pres = _b[_nl_1]
local se_pres = _se[_nl_1]	


rreg hits L* F* D , tune(6)	//re-run the pooled regression to save betas and graph
mat lags = J(`L',4,0)

foreach k of num 1(1)`L'{
	mat lags[`k',1] = _b[L`k']
	mat lags[`k',2] = _se[L`k']
	mat lags[`k',3] = _b[F`k']
	mat lags[`k',4] = _se[F`k']
}
*graph
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
	g min = beta - 1.96 * se
	g max = beta + 1.96 * se	
	g state = "0"

	tempfile pooled
	save `pooled'
	sort L
	tw rarea max min L , col(green*.50) || line beta L , lc(dkgreen) xline(0, lc(black) lp(dash) lw(thin)) ///
	legend(off) yline(0, lc(black) lw(thin)) xtitle(days) ytitle(betas) name(dum2, replace)	
	g w = 1/se
	
	outsheet L beta se min max w using pooled_inter_month.csv , comma replace
	save "saves/pooled_beta_inter", replace	
restore

*Heterogeneity over space
foreach s in `state_list' {
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
* Save full matrix per regions (ED fig 4)
preserve
	use `pooled', clear
	foreach s in `state_list' {
		append using `beta`s''
	}
	save "saves/state_beta_month", replace	
restore

*Heterogeneity over time

foreach Y of num 2008/2018 {
quietly {
	cap drop phase_cons
	reg hits L* F* D  if year_s == `Y' , nocons	
	
	mat lags = J(`L',5,0)
	foreach k of num `step'(`step')`L'{
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
		
		g min = beta - 1.96 * se
		g max = beta + 1.96 * se	
		g preiod = "`Y'"
		tempfile beta`Y'
		save `beta`Y''		
	restore
}
}

use `pooled', clear
foreach Y of num 2008/2018 {
	append using `beta`Y''
}
save "saves/time_beta_step1", replace	

rename period year
destring year, replace
drop if year == 0
g grp = L < 0
replace grp = 2 if L == 0
replace grp = 3 if L >0	
drop if abs(L) == 90
collapse (sum) beta, by(year grp)
egen totint = total(beta), by(year)
replace beta = beta / totint
g tag = "month" 

save "saves/share_over_time_month.dta", replace

