cd "~/Dropbox/temporal_awareness"
do "code_stata/fig2_fit_day"
do "code_stata/fig2_fit_month"
do "code_stata/fig2_fit_year"

graph combine proj2015 proj2016 m1 m2 g1 g2, rows(6) xcomm  ysize(10)
graph export "~/Dropbox/Sebastien-Leo-Sol/time writing/figures/fig2/fig2_raw.pdf",replace
graph drop _all
