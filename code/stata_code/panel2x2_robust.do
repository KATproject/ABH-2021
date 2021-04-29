cap cd "~\Dropbox\temporal_awareness"
import delim "saves/fitted_kernel_daily_robustness.csv", clear
drop v1
set scheme covid19

sum pred0
g rescale0 = 100 * pred0 / r(max)

foreach fe of num 1/4 {
sum pred`fe'
g rescale`fe' = 100 * pred`fe' / r(max)
tw (line rescale`fe' l, lc(red) lw(vthin)) ///
(line rescale0 l, lc(black) lw(vthin) legend(off) name(g`fe', replace) xtitle(""))
}

graph combine g1 g2 g3 g4, rows(2) ycomm

graph export "paper/figures/robust_kernel.pdf", replace
