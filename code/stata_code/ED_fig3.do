use "saves/BS_kat_1234.dta", clear

g cat = L < 0
replace cat = 2 if L == 0
replace cat = 3 if L > 0

drop if abs(L) == 35

collapse (sum) beta, by(ite cat)
reshape wide beta, i(ite) j(cat)
rename (beta1 beta2 beta3) (past pres fut)
g total = past + pres + fut
replace past = past / tot
replace pres = pres / tot
replace fut = fut / tot

tw hist past, bc(ebblue) || hist fut, bc(red) || hist pres, bc(dkgreen) xtitle(attention %) ///
legend(lab(1 "past") lab(2 "future") lab(3 "present") order(3 2 1) row(1) region(ls(none))) 
graph export "BS_share.pdf", replace


foreach period in "past" "fut" "pres" {
qui sum `period'
local m = round(r(mean),0.001)
local `period'CI5 = round(r(mean) - 1.96*r(sd),0.001)
local `period'CI95 = round(r(mean) + 1.96*r(sd),0.001)	
di "`period' avg `m' CI [``period'CI5' , ``period'CI95']"
}



import delim "saves/fitted_kernel_daily_BS.csv", clear
drop v1

collapse pred*, by(l)

reshape long pred, i(l) j(ite)

egen max = max(pred), by(ite)

replace pred = max / pred


collapse ratio = pred (sd) sd = pred, by(l)
g min = ratio - 1.96*sd
g max = ratio + 1.96*sd

g L = abs(l)

tw rarea max min l if l >=0 & abs(l) <200, col(maroon*.5) fi(100)  ///
|| rarea max min L if l <= 0 & abs(l) <200, col(ebblue*.5)  fi(100)   ///
|| line ratio l if l >=0 & abs(l) <200, lc(maroon) xtitle(Days, height(5)) ///
|| line ratio L if l <= 0 & abs(l) <200, lc(ebblue) ytitle(Inverse ratio of attention to the present) ///
legend(pos(11) ring(0) order(3 4) lab(3 "Future") lab(4 "Past") region(ls(none))) 

graph export "inverse_ratio.pdf", replace

