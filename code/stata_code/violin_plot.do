/*
Author: S. Annan-Phan
Date: 8/4/2020
Script to generate series of country level KAT 
*/

cap cd "~\Dropbox\temporal_awareness"
cap cd "/shares/lab/sphan/Research_projects/temporal_awareness"

*Set width of KAT for each type of querry
local width_month = 90 
local width_day = 90
local width_year = 240

foreach lab in "day" "month" "year" {
	
	*Crate empty matrix for predictions
	clear all
	local obs = `width_`lab'' * 2 + 1
	set obs `obs'
	egen L = seq()
	replace L = L - `width_`lab'' - 1
	tempfile basis
	save `basis'
	*Load estimated betas to be smoothed
	use "saves/state_beta_`lab'", clear
	g w = 1 / (se)
	drop if state == "0"
	levelsof state, local(list_state)
	local i = 1
	*Optimal smoothing for each state (restricted rational to ensure smoothness)
	foreach s in `list_state' {
		preserve
			keep if state == "`s'"
			append using `basis'
			*Past		
			nl (beta = ({a})/(1 + {k} * (L))) if L <= 0 [aweight=w]
			predictnl past = (_coef[/a]) / (1 + _coef[/k]*L), ci(low_p top_p)
			replace past = . if L > 0
			*Future
			nl (beta = ({a})/(1 + {k} * (L))) if L >= 0 [aweight=w]
			predictnl future = (_coef[/a]) / (1 + _coef[/k]*L), ci(low_f top_f)
			replace future = . if L < 0	
			*Clean and save
			sort L
			g fit = past
			replace fit = fut if L>0
			replace fit = 0 if fit < 0	
			keep L beta fit* state se min max
			replace state = "`s'" if state == ""
			tempfile f`i'
			save `f`i''
			local i = `i' + 1
		restore
	}

	use `f1', clear
	local end = `i' - 1
	foreach j of num 2/`end'{
		append using `f`j''
	}


	keep if beta == .
	keep L fit state

	*Rescale the raw value to plot each graph on top of each others
	sort state L	
	egen seq = seq(), by(state)
	egen top = max(fit) 
	g bottom = 0

	replace top = . if seq != 1
	g scal = sum(top)
	replace fit = fit + scal 
	replace bottom = bottom + scal
	egen labval = min(fit), by(state)

	levelsof state, local(list_s)
	levelsof labval, local(list_labval)
	local list_lab = ""
	foreach s in `list_s' {
		sum labval if state == "`s'"
		local lval = r(mean)
		local list_lab = `"`list_lab' `lval' "`s'" "'
	}

	foreach s in `list_s' {
		local gph `gph' || rarea fit bottom L if state == "`s'" & abs(L) < `width_`lab'', lc(black) fc(black) fi(100)
	}

	sum scal
	local top = r(max)
	tw rarea fit bottom L if state == "DE" & abs(L) < `width_`lab'' `gph' legend(off) ///
	ysize(20) ylabel(`list_lab', angle(0) noticks) ytitle("") xtitle("days") ///
	name(all, replace)
	graph export "violin_plot_`lab'.pdf", replace

}
