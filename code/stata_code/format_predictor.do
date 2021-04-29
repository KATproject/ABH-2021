cap cd "~\Dropbox\temporal_awareness"

local i = 1
foreach target in "yearly" "daily" "monthly" {
	import delim "saves/fitted_kernel_`target'.csv", clear varn(1)
	rename pred pred_rational
	collapse pred* beta, by(l)
	if "`target'" == "yearly" {
		replace l = l / 7
	}
	g name = "_F"
	replace name = "_L" if l < 0
	replace name = "_D" if l == 0
	replace l = abs(l)
	tostring l, replace
	replace name = name + l
	drop l
	g I = "`target'"
	reshape wide pred* beta, i(I) j(name) string
	tempfile f`i'
	save `f`i''
	local i = `i' + 1
}

use `f1', clear
append using `f2'
append using `f3'


save "saves/fitted_predictors.dta", replace
