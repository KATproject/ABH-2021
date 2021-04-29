
/*
Author: S. Annan-Phan
Date: 8/4/2020
Script producing R2 with rescaling factor
*/

* Identify region-by-search insample
egen ID = group(keyw geo) if pred_rational != . 
sum ID
local M = r(max)

mat scals = J(`M',2,0)
local k = 1
* Loop over each ID to find optimal rescaling factors (minimizing RSS)
forvalues ID = 1/`M' {
	di round(`ID'/`M',0.01)
	quietly {
	preserve
		keep if ID == `ID'
		foreach scalar of num 0.05(0.05)1.5{
			local name = floor(`scalar' * 100)
			g pred_`name' = pred_rational * `scalar'
			g rs`name' = (hits - pred_`name')^2
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
}
* Save matrix of scalar as dataset
preserve
	clear
	svmat scals
	rename (scals1 scals2) (ID sc_factor)
	tempfile SCL
	save `SCL'
restore
* Merge and re-scale
merge m:1 ID using `SCL'
replace sc_f = sc_f / 100
replace sc_f = 1 if sc_f == .

* Compute R2 with rescaled fact
preserve
	keep if pred_rational != .
	g pred_scl = sc_f * pred_rational 

	egen y_b = mean(hits)
	g ssr = (hits - pred_scl)^2
	g sst = (hits - y_b)^2
	egen SSR = sum(ssr)
	egen SST = sum(sst)
	g r2 = 1 - (SSR/SST)
	local rsquare = r2[1]
	di "The R2 for rescaled prediction is `rsquare'"
restore
