********************************************************************************
**																			  **
**				Assigntment # 4 - Causal Inference & Research Design		  **
**																			  **
********************************************************************************

cd "C:\Users\Nicolás Urdaneta\Documents\GitHub\RDD"
graph set window fontface "Times new roman"

use "https://github.com/scunning1975/causal-inference-class/raw/master/hansen_dwi", clear


*	____________________________________________________________________________
*
*									Replication
*	____________________________________________________________________________


label var male "Male"
label var aged "Age"
label var white "White"
label var acc "Accidents"

********************************************************************************
**								1. bac1 dummy								  **
********************************************************************************

gen d_bac1=(bac1>=0.08) & bac1!=.

********************************************************************************
**							2. Histogram, fig1								  **
********************************************************************************



hist bac1, addplot(pci 0 0.08 2000 0.08) bin(800) color(gs10) sch(tufte) 		///
 xtitle(BAC) frequency legend(off)
gr export "Figures/Hist.png", replace


* McCrary density test

*net install rddensity, from(https://sites.google.com/site/rdpackages/rddensity/stata) replace
*net install lpdensity, from(https://sites.google.com/site/nppackages/lpdensity/stata) replace

matrix McCrary = J(3,2,.)
local z=1
foreach x in 0.04 0.06 0.08 {
	rddensity bac1, c(`x') all
	matrix McCrary[`z',1] =`x'
	matrix McCrary[`z',2] = round(e(pv_q),0.0001)
	local z = `z'+1
}
mat l McCrary
estout matrix(McCrary) using "Tables/McCrary.tex", replace


********************************************************************************
**						  3. Balance test, table 2							  **
********************************************************************************

gen bac1_c = bac1-0.08

global covariates male white age acc

foreach var of global covariates{

	rdrobust `var' bac1 if abs(bac1_c)<=0.05, c(0.08) h(0.05 0.05) vce(hc0) p(1) kernel(uniform)
	est store reg_`var'
	sum `var' if round(bac1,0.001)==0.079
	estadd r(mean)
	*ereturn list

} 

#delimit ;
	cap n estout reg_male reg_white reg_age reg_acc using "Tables/table1.tex",
		style(tex) label notype
		cells((b(star fmt(%9.3f))) (se(fmt(%9.3f)par))) 		
		stats(N mean,
			labels("Observations" "Mean (at 0.079)")
			fmt(%9.0fc %9.2fc 2)) rename(RD_Estimate DUI)
			replace noabbrev starlevels(* 0.10 ** 0.05 *** 0.01) 
			title("Regression Disctontinuity Estimates for the Effect \\ of Exceeding BAC Thresholds on Predetermined Characteristics")   
			collabels(none) eqlabels(none) mlabels(none) mgroups(none) substitute(_ \_)
			prehead("\begin{center}" "\begin{table}[H] \centering" "\captionsetup{justification=centering}" "\caption{@title}" "\begin{threeparttable}" "\begingroup"
					"\setlength{\tabcolsep}{10pt}" 
					"\renewcommand{\arraystretch}{1.5}"" % Default value: 1" "\begin{tabular}{l*{@E}{c}}"
					"\toprule"
					"& \multicolumn{3}{>{\centering\arraybackslash}m{150pt}}{Driver demographic characteristics} & \\ \cline{2-4}"
					"& Male & White & Age & Accident \\"
					"Characteristics & (1) & (2) & (3) & (4) \\")
			prefoot("\\")
		posthead("\hline")
		postfoot("\hline \hline" "\end{tabular}" "\endgroup" "\begin{tablenotes} \centering" "\small" 
			"\item This table contains egression discontinuity based estimates of 
		 the effect of having BAC above the DUI threshold on recidivism. 
		 Panel A contains estimates with a bandwidth of 0.05 while Panel B 
		 has a bandwidth of 0.025, with all regressions utilizing a rectangular
		  kernel for weighting. 
		  Robust standard errors. * p$<$0.10, ** p$<$0.05, *** p$<$0.01" 
		  "\end{tablenotes}" "\end{threeparttable}" "\end{table}" "\end{center}");
#delimit cr

********************************************************************************
**						  4. Balance test, Fig 2							  **
********************************************************************************
	save temp, replace

	
foreach var of global covariates{
	forvalues z=1/2{

	cap drop ciupper`z'_`var'
	cap drop cilower`z'_`var'

	cap drop x*
	cap drop s*
	cap drop se*

	use temp, clear
	qui lpoly `var' bac1 if bac1<0.08 & bac1>=0.03, kernel(rectangle) bwidth(.05) degree(`z')  generate(x s) se(se) 
	keep x s se 
	drop if x==.
	gen ciupper`z'_`var'=s+1.96*se
	gen cilower`z'_`var'=s-1.96*se
	save RD`var'`z', replace


	use temp, clear
	qui lpoly `var' bac1 if bac1>=0.08 & bac1<0.15, kernel(rectangle) bwidth(.05) degree(`z')  generate(x s) se(se) 
	keep x s se 
	drop if x==.
	gen ciupper`z'_`var'=s+1.96*se
	gen cilower`z'_`var'=s-1.96*se
	append using RD`var'`z'
	save RD`var'`z', replace

	use temp, clear
	qui lpoly `var' bac1 if bac1>=0.15 & bac1<=0.2, kernel(rectangle) bwidth(.05) degree(`z')  generate(x s) se(se) 
	keep x s se 
	drop if x==.
	gen ciupper`z'_`var'=s+1.96*se
	gen cilower`z'_`var'=s-1.96*se
	append using RD`var'`z'
	save RD`var'`z', replace

	}
}


use temp, clear

keep if bac1>=0.03 & bac1<=0.2


	gen bin=.
	foreach X of num 0.03(.001).2 {
		di "`X'"
		qui replace bin=`X' if (bac1>=`X' & bac1<(`X'+.001) )
	}

	gen count=1
	collapse $covariates bac1 (sum) count, by(bin)

foreach var of global covariates{
	forvalues z=1/2{
		append using RD`var'`z'
	}
}

foreach var of global covariates{

	if "`var'"=="acc" local title = "Accident at scene"
	if "`var'"=="male" local title = "Male"
	if "`var'"=="age" local title = "Age"
	if "`var'"=="white" local title = "White"

	tw (scatter `var' bac1, mc(gs8) msize(small)) 											///
	(lfit `var' bac1 if bac1<0.08 [aw=count], lp(solid)) 						///
	(lfit `var' bac1 if bac1>=0.08 & bac1<0.15 [aw=count], lp(solid))			///
	(lfit `var' bac1 if bac1>=0.15 [aw=count], lp(solid))						///
	(connected ciupper1_`var' cilower1_`var' x if x<0.08,  sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper1_`var' cilower1_`var' x if x>=0.08 & x<0.15,  sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper1_`var' cilower1_`var' x if x>=0.15, sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	, xline(0.08) xline(0.15) nodraw xtitle("BAC")		///
	name(`var'_p1, replace) sch(tufte) title(`title') ///
	legend(order(2 "Linear fit" 1 "`title'") r(1)) ///
	xscale(range(0.028 0.202)) xlabel(0.05 0.1 0.15 0.2)

	tw (scatter `var' bac1, mc(gs8) msize(small)) 											///
	(qfit `var' bac1 if bac1<0.08 [aw=count], lp(solid) lc("20 20 140")) 						///
	(qfit `var' bac1 if bac1>=0.08 & bac1<0.15 [aw=count], lp(solid) lc("20 20 140")) 			///
	(qfit `var' bac1 if bac1>=0.15 [aw=count], lp(solid) lc("20 20 140"))     	///
	(connected ciupper2_`var' cilower2_`var' x if x<0.08,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper2_`var' cilower2_`var' x if x>=0.08 & x<0.15,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper2_`var' cilower2_`var' x if x>=0.15,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	, xline(0.08) xline(0.15) nodraw xtitle("BAC")		///
	name(`var'_p2, replace) sch(tufte) title(`title') ///
	legend(order(2 "Quadratic fit" 1 "`title'") r(1)) ///
	xscale(range(0.028 0.202)) xlabel(0.05 0.1 0.15 0.2)

}


gr combine acc_p1 male_p1 white_p1 age_p1, sch(tufte) name(Linear, replace)
gr export "Figures/Fig1a.png", replace
gr combine acc_p2 male_p2 white_p2 age_p2, sch(tufte) name(quadratic, replace)
gr export "Figures/Fig1b.png", replace



********************************************************************************
**						  5. Estimate with recid							  **
********************************************************************************
use temp, clear
global controls i.male i.white aged i.year

gen bac1c_2=(bac1-0.08)^2

reg recidivism d_bac1 bac1 $controls if abs(bac1_c)<=0.05, r
est store reg_A1
sum recidivism
estadd r(mean)

reg recidivism d_bac1 bac1_c i.d_bac1#c.bac1_c $controls if abs(bac1_c)<=0.05, r
est store reg_A2
sum recidivism
estadd r(mean)

reg recidivism d_bac1 bac1_c bac1c_2 i.d_bac1#c.bac1_c i.d_bac1#c.bac1c_2 $controls if abs(bac1_c)<=0.05, r
est store reg_A3
sum recidivism
estadd r(mean)


reg recidivism d_bac1 bac1_c $controls if abs(bac1_c)<=0.025, r
est store reg_B1
sum recidivism
estadd r(mean)

reg recidivism d_bac1 bac1_c i.d_bac1#c.bac1_c $controls if abs(bac1_c)<=0.025, r
est store reg_B2
sum recidivism
estadd r(mean)

reg recidivism d_bac1 bac1_c bac1c_2 i.d_bac1#c.bac1_c i.d_bac1#c.bac1c_2 $controls if abs(bac1_c)<=0.025, r
est store reg_B3
sum recidivism
estadd r(mean)


#delimit ;
	cap n estout reg_A1 reg_A2 reg_A3 using "Tables/table2a.tex",
		style(tex) label notype
		cells((b(star fmt(%9.3f))) (se(fmt(%9.3f)par))) 		
		stats(N mean,
			labels("Observations" "Mean")
			fmt(%9.0fc %9.2fc 2)) rename(d_bac1 DUI) keep(DUI)
			replace noabbrev starlevels(* 0.10 ** 0.05 *** 0.01) 
			title("Regression Disctontinuity Estimates for the Effect \\ of Exceeding BAC Thresholds on Recidivism")   
			collabels(none) eqlabels(none) mlabels(none) mgroups(none) substitute(_ \_)
			prehead("\begin{center}" "\begin{table}[H] \centering" "\captionsetup{justification=centering}" "\caption{@title}" "\begin{threeparttable}" "\begingroup"
					" % Default value: 1" "\begin{tabular}{l*{@E}{>{\centering\arraybackslash}m{3cm}}}"
					"\toprule"
					"& Linear BAC & Linear BAC with interaction & Quadratic BAC with interaction \\"
					" & (1) & (2) & (3) \\ \hline"
					"Panel A. BAC $\in$ [0.03,0.013] & & & \\")
			prefoot("\\");
#delimit cr

#delimit ;
	cap n estout reg_B1 reg_B2 reg_B3 using "Tables/table2b.tex",
		style(tex) label notype
		cells((b(star fmt(%9.3f))) (se(fmt(%9.3f)par))) 		
		stats(N mean,
			labels("Observations" "Mean")
			fmt(%9.0fc %9.2fc 2)) rename(d_bac1 DUI) keep(DUI)
			replace noabbrev starlevels(* 0.10 ** 0.05 *** 0.01) 
			title("Regression Disctontinuity Estimates for the Effect \\ of Exceeding BAC Thresholds on Recidivism")   
			collabels(none) eqlabels(none) mlabels(none) mgroups(none) substitute(_ \_)
			prehead("\\ "
					"Panel B. BAC $\in$ [0.055,0.105] & & & \\ ")
			prefoot("\\")
		postfoot("\hline \hline" "\end{tabular}" "\endgroup" "\begin{tablenotes} \centering" "\small"
		 "\item This table contains egression discontinuity based estimates of 
		 the effect of having BAC above the DUI threshold on recidivism. 
		 Panel A contains estimates with a bandwidth of 0.05 while Panel B 
		 has a bandwidth of 0.025, with all regressions utilizing a rectangular
		  kernel for weighting. Controls include indicators for year, race, 
		  gender, and age of the offender (not county). 
		  Robust standard errors in parenthesis. 
		  * p$<$0.10, ** p$<$0.05, *** p$<$0.01" "\end{tablenotes}" "\end{threeparttable}" "\end{table}" "\end{center}");
#delimit cr


********************************************************************************
**						  6. Outcomes, Fig 3							  **
********************************************************************************
use temp, clear
forvalues z=1/2{

	cap drop ciupper
	cap drop cilower

	cap drop x*
	cap drop s*
	cap drop se*

	use temp, clear
	qui lpoly recidivism bac1 if bac1<0.08 & bac1>=0.03, kernel(rectangle) bwidth(.05) degree(`z')  generate(x s) se(se) 
	keep x s se 
	drop if x==.
	gen ciupper`z'=s+1.96*se
	gen cilower`z'=s-1.96*se
	save RD`z', replace


	use temp, clear
	qui lpoly recidivism bac1 if bac1>=0.08 & bac1<0.15, kernel(rectangle) bwidth(.05) degree(`z')  generate(x s) se(se) 
	keep x s se 
	drop if x==.
	gen ciupper`z'=s+1.96*se
	gen cilower`z'=s-1.96*se
	append using RD`z'
	save RD`z', replace

}



use temp, clear

keep if bac1>=0.03 & bac1<=0.15


	gen bin=.
	foreach X of num 0.03(.001).15 {
		di "`X'"
		qui replace bin=`X' if (bac1>=`X' & bac1<(`X'+.001) )
	}

	gen count=1
	collapse recidivism bac1 (sum) count, by(bin)

forvalues z=1/2{
	append using RD`z'
}

	tw (scatter recidivism bac1, mc(gs8) msize(small)) 											///
	(lfit recidivism bac1 if bac1<0.08 [aw=count], lp(solid)) 						///
	(lfit recidivism bac1 if bac1>=0.08 & bac1<0.15 [aw=count], lp(solid))			///
	(lfit recidivism bac1 if bac1>=0.15 [aw=count], lp(solid))						///
	(connected ciupper1 cilower1 x if x<0.08,  sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper1 cilower1 x if x>=0.08 & x<0.15,  sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper1 cilower1 x if x>=0.15, sort msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	, xline(0.08)  xtitle("BAC") title("Linear fit")		///
	name(recidivism_p1, replace) sch(tufte) ///
	legend(order(2 "Linear fit" 1 "Recidivism") r(1)) ///
	xscale(range(0.028 0.152)) xlabel(0.05 0.1 0.15) ///
	ylabel(0.05 0.075 0.1 0.125 0.15)

	tw (scatter recidivism bac1, mc(gs8) msize(small)) 											///
	(qfit recidivism bac1 if bac1<0.08 [aw=count], lp(solid) lc("20 20 140")) 						///
	(qfit recidivism bac1 if bac1>=0.08 & bac1<0.15 [aw=count], lp(solid) lc("20 20 140")) 			///
	(qfit recidivism bac1 if bac1>=0.15 [aw=count], lp(solid) lc("20 20 140"))     	///
	(connected ciupper2 cilower2 x if x<0.08,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper2 cilower2 x if x>=0.08 & x<0.15,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	(connected ciupper2 cilower2 x if x>=0.15,  msymbol(none none) clcolor(gs2 gs2) clpat(shortdash shortdash)) ///
	, xline(0.08)  xtitle("BAC") title("Quadratic fit")		///
	name(recidivism_p2, replace) sch(tufte)  ///
	legend(order(2 "Quadratic fit" 1 "Recidivism") r(1)) ///
	xscale(range(0.028 0.152)) xlabel(0.05 0.1 0.15) ///
	ylabel(0.05 0.075 0.1 0.125 0.15)


gr combine recidivism_p1 recidivism_p2, sch(tufte) 
gr display, xsize(18) ysize(8)
gr export "Figures/Fig2.png", replace
