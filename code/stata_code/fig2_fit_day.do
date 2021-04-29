use "saves/fitted_predictors", clear
keep if I == "daily"

keep pred_rational*
rename pred_rational_* fit*
rename fitD0 fitD
g j = 1
tempfile predictor
save `predictor'

import delim "placebo_comparison_new.txt",clear
g date2 = date(date,"YMD",2000)
drop date 
rename date2 date
format date %td
g year = year(date)
g month = month(date)
g day = day(date)
keep if keyw == "valentine's day" | keyw == "diwali"
egen max = max(hits), by(keyw)
replace hits = (hits / max) * 100
g D = 0
replace D = 1 if day == 14 & month == 2 & keyw == "valentine's day"
replace D = 1 if date == mdy(10,28,2008) & keyw == "diwali" 
replace D = 1 if date == mdy(10,17,2009) & keyw == "diwali" 
replace D = 1 if date == mdy(11,05,2010) & keyw == "diwali" 
replace D = 1 if date == mdy(10,26,2011) & keyw == "diwali" 
replace D = 1 if date == mdy(11,13,2012) & keyw == "diwali" 
replace D = 1 if date == mdy(11,03,2013) & keyw == "diwali" 
replace D = 1 if date == mdy(10,23,2014) & keyw == "diwali" 
replace D = 1 if date == mdy(11,11,2015) & keyw == "diwali" 
replace D = 1 if date == mdy(10,30,2016) & keyw == "diwali" 
replace D = 1 if date == mdy(10,19,2017) & keyw == "diwali" 
replace D = 1 if date == mdy(11,07,2018) & keyw == "diwali" 

tempfile baseD
save `baseD'

local W = 90

egen i = group(keyw)
xtset i date

foreach k of num 1(1)`W'{
	g L`k' = D[_n-`k'] 
	g F`k' = D[_n+`k'] 
	replace F`k' = 0 if  F`k' == .
}

reg hits L* F* D,nocons
predict pred
replace pred = . if pred == 0

g j = 1
merge m:1 j using `predictor', nogen

local cmd = "D * fitD"
foreach i of num 1/`W' {
	local cmd = "`cmd' + F`i' * fitF`i'"
	local cmd = "`cmd' + L`i' * fitL`i'"	
}
g pred_smooth = `cmd'
replace pred_smooth = . if pred == .
		
egen ID = group(keyw year) if pred_smooth != .
qui sum ID
local M = r(max)

mat scals = J(`M',2,0)
local k = 1
foreach ID of num 1/`M' {
	preserve
	keep if ID == `ID'
	foreach scalar of num 0.05(0.05)1.5{
		local name = floor(`scalar' * 100)
		g predi_`name' = pred_smooth * `scalar'
		g rs`name' = (hits - predi_`name')^2
		egen rss`name' = sum(rs`name')
		drop rs`name'
	}

		
		keep rss*
		keep if _n == 1
		g i = 1
		reshape long rss, i(i) j(scal)
		g scal_graph = scal / 100
		replace rss = sqrt(rss/_N)
		egen min = min(rss)
		g min_id = rss == min
		sum scal if min_id == 1
		mat scals[`k',1] = `ID'
		mat scals[`k',2] = r(mean)
		local k = `k' + 1
	restore
}

preserve
	clear
	svmat scals
	rename (scals1 scals2) (ID sc_factor)
	tempfile SCL
	save `SCL'
restore
merge m:1 ID using `SCL'
replace sc_f = sc_f / 100
replace sc_f = 1 if sc_f == .

// GET R2 with rescaled fact
preserve
keep if pred_smooth != .

replace pred_smooth = sc_f * pred_smooth

egen y_b = mean(hits)
g ssr = (hits - pred_smooth)^2
g sst = (hits - y_b)^2
egen SSR = sum(ssr)
egen SST = sum(sst)
g r2 = 1 - (SSR/SST)
di r2[1]
restore

g pred_rescale = pred_s * sc_f
sort keyw y date

replace pred_r = 0 if pred_r == .
replace pred_r = 100 if pred_r > 100

format date %tdCY
tw line hits date if i == 1 ///
|| line pred_rescale date if i == 1, cmissing(n) lc(red) lw(vthin) ///
 legend(off) xtitle("")  name(g1, replace) 

 tw line hits date if i == 2 ///
|| line pred_rescale date if i == 2, cmissing(n) lc(red) lw(vthin) ///
 legend(off) xtitle("")  name(g2, replace) 
 
graph combine g1 g2, rows(2) xcomm 

format date %tdN-D-CY
keep date keyw hits pred_rescale
 
outsheet * using "saves/fig2_prediction_day.csv", comma replace
 