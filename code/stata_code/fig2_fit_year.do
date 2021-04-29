use "saves/fitted_predictors", clear
keep if I == "yearly"

keep pred_rational*
rename pred_rational_* fit*
rename fitD0 fitD
g j = 1
tempfile predictor
save `predictor'



foreach yy of num 2015 2016 {
import delim "placebo_comparison_new.txt",clear
g date2 = date(date,"YMD",2000)
drop date 
rename date2 date
format date %td
g year = year(date)

keep if keyw == "`yy'"  
egen max = max(hits)
replace hits = (hits / max) * 100
g D = year == `yy'

tempfile baseD
save `baseD'

local K = 7
local W = floor(245/ `K')
di `W'

egen w = seq(), block(`K')
collapse D hits date, by(w)
replace date = floor(date)
tset date
cap drop L* F*

foreach k of num 1(1)`W'{
	g L`k' = D[_n-`k'] 
	g F`k' = D[_n+`k'] 
	replace F`k' = 0 if  F`k' == .
}

reg hits L* F* D,nocons
cap drop pred
predict pred

replace pred = . if pred == 0
format date %tdCY
tempfile Pred
save `Pred'


use `baseD', clear
merge 1:1 date using `Pred'

preserve
use `Pred', clear
g j = 1
merge m:1 j using `predictor', nogen

local cmd = "D * fitD"
foreach i of num 1/34 {
	local cmd = "`cmd' + F`i' * fitF`i'"
	local cmd = "`cmd' + L`i' * fitL`i'"	
}

g pred_smooth = `cmd'

save `Pred', replace
restore

merge 1:1 date using `Pred', nogen
replace pred_smooth = . if pred == .

foreach scalar of num 0.05(0.05)1.5{
	local name = floor(`scalar' * 100)
	g pred_s`name' = pred_smooth * `scalar'
	g rs`name' = (hits - pred_s`name')^2
	egen rss`name' = sum(rs`name')
	drop rs`name'
}

preserve
	keep rss*
	keep if _n == 1
	g i = 1
	reshape long rss, i(i) j(scal)
	g scal_graph = scal / 100
	replace rss = sqrt(rss/_N)
	egen min = min(rss)
	g min_id = rss == min
	sum scal if min_id == 1
	local best_scalar = r(mean)
restore


replace pred_smooth = pred_s`best_scalar'	
g RSS = rss`best_scalar'

preserve
	keep if pred_smooth != .
	keep RSS date hits pred_smooth
	egen y_bar = mean(hits)
	g ss = (hits - y_bar)^2
	egen SST = sum(ss) 
	g R2 = 1 - (RSS/SST)
	qui sum R2
	local r2 = round(r(mean),0.01)
restore

di `r2'

replace pred_smooth = 0 if pred_smooth == . & (L35 == 0 | F35 == 0)
tw line hits date || line pred_smooth date,  lc(red) lw(vthin)  ///
 legend(off) xtitle("")  name(proj`yy', replace) 

format date %tdN-D-CY
rename pred_smooth pred_rescale
keep date keyw hits pred_rescale
tempfile f`yy'
save `f`yy'', replace 
 
}

use `f2015', clear
append using `f2016'

outsheet * using "saves/fig2_prediction_year.csv", comma replace
