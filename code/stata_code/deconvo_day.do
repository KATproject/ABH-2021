/*
Author: S. Annan-Phan
Date: 8/4/2020
This do file run the deconvolution model for "day" related querries
*/

//Set path to working directory
cd "~\Dropbox\temporal_awareness"
*Number of lag/leads (weeks), lot of colinearity above 25
local L = 90
local level = "international_holidays"

*data preparation
import delim "daily_interest_international_holidays.txt", clear varn(1)

egen id = group(geo keyw search_t)
*clean date
replace date = substr(date,1,10)
g date2 = date(date,"YMD",2000)
drop date
rename date2 date
format date %td
g week = week(date)
g year = year(date)
g month = month(date)
egen country = group(geo)

g date_search = date(search_t,"YMD",2000)
format date_s %td
cap replace hits = "0" if hits == "<1"
cap destring hits, replace force
drop search_t
save "international_holidays.dta", replace

drop if geo == "PT" | geo == "PH" | geo == "IT" 
keep if year >= 2008

sort geo keyw date

*time
egen t = seq(), by(keyw geo)
sort geo keyw date
*panel
xtset id t
*Dummy if t is in S (dirac function)
g D = date == date_s
sort geo keyw date
g dow = dow(date)
* lag and leads of the dirac function
quietly{
foreach k of num 1(1)`L'{
	g L`k' = D[_n-`k']
	g F`k' = D[_n+`k']
}
}

******************
*    Analysis    *
******************
save "saves/data_day_reg", replace	

*Matrix of estimates 

* Pooled reg
reg hits L* F* D , nocons

*create command of cumulative attention
local LC_past = "_b[L1]"
local LC_fut = "_b[F1]"
foreach q of num 2/89 {
	local LC_past = "`LC_past' + _b[L`q']"
	local LC_fut = "`LC_fut' + _b[F`q']"
}


*statistics for fig3
nlcom (`LC_past') / (`LC_past' + `LC_fut' + _b[D]),post
local b_past = _b[_nl_1]
local se_past = _se[_nl_1]

reg hits L* F* D , nocons
nlcom (`LC_fut') / (`LC_past' + `LC_fut' + _b[D]),post
local b_fut = _b[_nl_1]
local se_fut = _se[_nl_1]

reg hits L* F* D , nocons
nlcom (_b[D]) / (`LC_past' + `LC_fut' + _b[D]),post
local b_pres = _b[_nl_1]
local se_pres = _se[_nl_1]

reg hits L* F* D , nocons
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
	drop if se == 0
	
	sum beta
	g min = beta - 1.96 * se
	g max = beta + 1.96 * se	
	g state = "0"
	tempfile pooled
	save `pooled'
	sort L
	g w = 1 /se

	
	outsheet L beta se min max w using pooled_inter_day.csv , comma replace
	tw rarea max min L , col(ebblue*.30) || line beta L , lc(navy) ///
	legend(off) yline(0, lc(black) lp(dash)) xtitle(days) ytitle(betas) name(dum2, replace)	
	graph export "saves/figure_analysis/pooled_inter_holidays.png", replace
restore

tempfile basedata
save `basedata'

* Space heterogeneity
levelsof(geo), local(state_list)
foreach s in `state_list'{
	reg hits L* F* D if geo == "`s'", nocons
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
		drop if se == 0
		
		sum beta
		g min = beta - 1.96 * se
		g max = beta + 1.96 * se	
		g state = "`s'"
		tempfile beta`s'
		save `beta`s''
	restore
}
preserve
	use `pooled', clear
	foreach s in `state_list' {
		append using `beta`s''
	}
	save "saves/state_beta_day", replace	
restore

* Time heterogeneity
foreach Y of num 2008/2018 {
quietly {
	reg hits L* F* D i.country if year == `Y', nocons

	mat lags = J(`L',5,0)
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
		drop if se == 0
		sum beta
		g min = beta - 1.96 * se
		g max = beta + 1.96 * se	
		g state = "`Y'"
		tempfile beta`Y'
		save `beta`Y''		
	restore
}
}


use `pooled', clear
local start_year = 2008 + 1
foreach Y of num 2008/2018{
	append using `beta`Y''
}
rename state period

g grp = L < 0
replace grp = 2 if L == 0
replace grp = 3 if L >0	
collapse (sum) beta, by(period grp)
egen totint = total(beta), by(period)
replace beta = beta / totint

destring period, replace

rename period year
drop if year == 0
tw scatter beta year if grp == 1, mc(red) connect(line) lc(red) lp(dash) ///
|| lfit beta year if grp == 1, lc(red) ///
|| scatter beta year if grp == 3, mc(ebblue) connect(line) lc(ebblue) lp(dash) ///
|| lfit beta year if grp == 3, lc(ebblue) ///
|| scatter beta year if grp == 2, mc(black) connect(line) lc(black) lp(dash) ///
|| lfit beta year if grp == 2, lc(black) xtitle("") ytitle("share of interest") ///
legend(lab(1 "past") lab(3 "futur") lab(5 "pres") order(1 3 5) ///
 region(lstyle(none)) rows(1)) title("daily")

graph export "saves/figure_analysis/share_over_time_daily.png", replace


g tag = "day" 

save "saves/share_over_time_day.dta", replace

