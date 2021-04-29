/*
Author: S. Annan-Phan
Date: 8/4/2020
do-file generating model prediction post regression and smoothed kernels
ED figure5
*/
cd "~\Dropbox\temporal_awareness"

****************************
*** First panel - Daily ****
****************************

use "saves/data_day_reg", clear	
reg hits L* F* D 
predict pred
g I = "daily"
merge m:1 I using "saves/fitted_predictors", keep(3) nogen

foreach ff in "poly" "hype" "expo" "rational" { 
	g pred_`ff' = D * pred_`ff'_D0
	foreach theta of num 1/90 {
		replace pred_`ff' = pred_`ff' + L`theta' * pred_`ff'_L`theta' + F`theta' * pred_`ff'_F`theta'
	}
}

*RMSE
preserve
	foreach ff in "poly" "hype" "expo" "rational" { 
		g s2_`ff' = (hits - pred_`ff')^2
	}
	keep s2*
	save "saves/SE_day", replace
restore

*Prediction fig SI
preserve
	g T = t * D 
	egen ttemp = max(T), by(id)
	replace T = t - ttemp
	keep if abs(T) <= 90
	collapse hits pred*, by(T)
	egen t = seq()
	tw scatter hits t, mc(black) ///
	|| line pred_hype t, lc(midblue) /// 	
	|| line pred_rational t, lc(midgreen) ///
	|| line pred_poly t, lc(gold) /// 
	|| line pred_expo t, lc(red) xtitle(Days) ///
	legend(rows(3) ring(0) pos(2) lab(1 "Hits") lab(2 "Rational") ///
	lab(3 "Hyperb") lab(4 "Poly") lab(5 "Expo") region(fc(gs15))) ///
	title(Holidays) name(fit, replace)		
	graph export "pred_daily.pdf", replace
restore
* R2
do code_stata/script_r2.do

******************************
*** Second panel - Monthly ***
******************************

use "saves/data_month_reg", clear	
rreg hits L* F* D , tune(6)	
predict pred
g I = "monthly"
merge m:1 I using "saves/fitted_predictors", keep(3) nogen

foreach ff in "poly" "hype" "expo" "rational" { 
	g pred_`ff' = D * pred_`ff'_D0 + _b[_cons] 
	foreach theta of num 1/90 {
		replace pred_`ff' = pred_`ff' + L`theta' * pred_`ff'_L`theta' + F`theta' * pred_`ff'_F`theta'
	}
}
*RMSE
preserve
	foreach ff in "poly" "hype" "expo" "rational" { 
		g s2_`ff' = (hits - pred_`ff')^2
	}
	keep s2*
	save "saves/SE_month", replace
restore

preserve
egen time_low = min(t) if D != 0, by(id)
egen t_low = mean(time_low), by(id)
egen time_up = max(t) if D != 0, by(id)
egen t_up = mean(time_up), by(id)

drop if t < t_low - 90
drop if t > t_up + 90
egen T = seq(), by(id)	
collapse hits pred pred_rational pred_hype pred_poly pred_expo, by(T)											
tw scatter hits T, mc(black) ///
|| line pred_rational T, lc(midgreen) ///
|| line pred_hype T, lc(midblue) /// 
|| line pred_poly T, lc(gold) /// 
|| line pred_expo T, lc(red) xtitle(Days) ///
legend(cols(1) ring(0) pos(2) lab(1 "Hits") lab(2 "Rational") ///
lab(3 "Hyperb") lab(4 "Poly") lab(5 "Expo") region(fc(gs15))) ///
title(Month) 
graph export "pred_monthly.pdf", replace
restore

* R2
do code_stata/script_r2.do

****************************
*** Third panel - Yearly ***
****************************
use "saves/data_year_reg", replace	
drop __*
keep if subsamp == 1
rreg hits L* F* D, tune(6)
predict pred
g I = "yearly"
merge m:1 I using "saves/fitted_predictors", keep(3) nogen
foreach ff in "poly" "hype" "expo" "rational" { 
	g pred_`ff' = D * pred_`ff'_D0  + _b[_cons]
	foreach theta of num 1/30 {
		replace pred_`ff' = pred_`ff' + L`theta' * pred_`ff'_L`theta' + F`theta' * pred_`ff'_F`theta'
	}
}

*RMSE
preserve
	foreach ff in "poly" "hype" "expo" "rational" { 
		g s2_`ff' = (hits - pred_`ff')^2
	}
	keep s2*
	save "saves/SE_year", replace
restore


preserve
	collapse hits pred pred_rational pred_hype pred_poly pred_expo, by(t)
	replace t = t * 7
	tw scatter hits t, mc(black)  ///
	|| line pred_rational t, lc(midgreen) lw(medthick) ///
	|| line pred_hype t, lc(midblue) lw(medthick) /// 
	|| line pred_poly t, lc(gold) lw(medthick) /// 
	|| line pred_expo t, lc(red) lw(medthick) xtitle(Days) /// 
	legend(rows(3) ring(0) pos(2) lab(1 "Hits") lab(2 "Rational") ///
	lab(3 "Hyperb") lab(4 "Poly") lab(5 "Expo") region(fc(gs15))) ///
	xtitle(days) ///
	legend(rows(2) ring(0) pos(2) lab(1 "Hits") lab(2 "Prediction") region(fc(gs15))) ///
	title(Year) 
	graph export "pred_yearly.pdf", replace
restore

*Oveall RMSE
use "saves/SE_day.dta", clear
append using "saves/SE_month"
append using "saves/SE_year"
collapse s2*
foreach fit in "poly" "hype" "expo" "rational" {
	g RMSE_`fit' = sqrt(s2_`fit')
}

* R2
do code_stata/script_r2.do



