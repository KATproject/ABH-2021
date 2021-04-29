
cap cd "~\Dropbox\temporal_awareness"
use "saves/share_over_time_day.dta", clear

global fig_folder = "~/Dropbox/Sebastien-Leo-Sol/time writing/figures"

append using "saves/share_over_time_month.dta"
append using "saves/share_over_time_year.dta"

tempname memhold
postfile `memhold' str18 panel slope SE using "saves/fig4_estimates.dta", replace


local lab1 = "Past"
local lab2 = "Present"
local lab3 = "Future"
preserve
keep if year >= 2009



reg beta year if grp == 1
post `memhold' ("past") (_b[year]) (_se[year]) 
reg beta year if grp == 2
post `memhold' ("present") (_b[year]) (_se[year]) 
reg beta year if grp == 3
post `memhold' ("futur") (_b[year]) (_se[year]) 
postclose `memhold'

foreach grp of num 1 2 3 {
	tw scatter beta year if grp == `grp' & tag == "day", mc(ebblue*.75)  lc(black) lp(dash) ///
	|| lfit beta year if grp == `grp' & tag == "day", lc(ebblue) lw(medthick) ///
	|| scatter beta year if grp == `grp' & tag == "month", mc(gold*.75)  lc(black) lp(dash) ///
	|| lfit beta year if grp == `grp' & tag == "month", lc(gold) lw(medthick) ///
	|| scatter beta year if grp == `grp' & tag == "year", mc(midgreen*.75) lc(black) lp(dash) ///
	|| lfit beta year if grp == `grp' & tag == "year", lc(midgreen) lw(medthick) ///
	|| scatter beta year if grp == `grp' & tag == "world", mc(black*.75) lc(black) lp(dash) ///
	|| lfit beta year if grp == `grp' & tag == "world", lc(black) lw(medthick) title(`lab`grp'') xtitle("") ytitle("share of attention") ///
	legend(lab(1 "Day") lab(3 "Month") lab(5 "year") order(1 3 5) ///
	 region(lstyle(none)) rows(1)) name(`lab`grp'',replace)
 }
graph combine Past Present Future, rows(1) ycomm xsize(14)
graph export "$fig_folder/fig4/fig4_time_serie_stata_raw.pdf", replace


use "saves/fig4_estimates.dta", clear
g CI_lb = slope - 1.965*SE
g CI_ub = slope + 1.965*SE
outsheet * using "$fig_folder/../results_details/fig4_trend_slope.csv", comma replace
