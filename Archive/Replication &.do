/*******************************************************************************
Author: Ken Norris

Date: Nov 2018

Last Update: Nov 2018

Description:

a.	Replications: I’d also like to see if we can replicate the Li et al 
	2014 and the Coglianese et al 2017 results with our data.

b.  Comparisons:  comparing elasticities from consumption data and traffic count
	data

c.  Extensions for various kinds of shocks: hurricane shocks and the like
	
	
*******************************************************************************/

/* A) Li et al replication

Baseline:

reg consumption on tax exclusive gasoline price &(1 + total state and federal
gasoline tax) including a vector of state level observables and state, year FE

Under the assumption that consumers bear the entire tax, 
the tax-exclusive price is not affected by a change in the tax rate, 
and dp /dr is equal to zero. (page 12 Li et al)

IV:



B) Coglianese et al. replication

Event Study around a tax increase/decrease comparing delta q(i,t) [consumption 
of gasoline] vs tax increase


*/

clear all
set more off

global dropbox "C:/Users/`c(username)'/Dropbox/Fuel Price Elasticity"

cd "$dropbox"

set scheme s1mono
set more off

/*

*Creating dta file for seasonal unemployment data from BLS here*
* See https://download.bls.gov/pub/time.series/la/la.txt for details on how to
* work with directory

import delimited "$dropbox/Data/la.data.3.txt", clear

replace series_id = subinstr(series_id, "LASST", "", 1)
replace period = subinstr(period, "M", "", 1)

keep if strpos(series_id, "00003")!=0
gen fips = substr(series_id, 1, 2)

drop series_id footnote_codes 
rename value unemployment_rate
rename period month
destring month, replace
destring fips, replace

save "$dropbox/Data/BLS Unemployment 1976 - 2018", replace

*/

* Census Population data here 

import delimited "$dropbox/Data/co-est00int-tot.csv", clear

tempfile pop20
save "`pop20'"

import delimited "$dropbox/Data/co-est2017-alldata", clear

merge 1:1 state county using `pop20', keep(matched) nogen

*2010 data has birth, death, residual, etc. data, just want to keep popn

keep state county stname ctyname popestim*

drop if county == 0

gen str5 fipscnty = string(state, "%02.0f") + string(county, "%03.0f")

tempfile population
save "`population'"

use "$dropbox/Data/cbsa2fipsxw.dta", clear

gen fipscnty = fipsstatecode + fipscountycode
keep cbsacode cbsatitle countycountyequivalent fipscnty

rename countycountyequivalent ctyname
rename cbsacode cbsa

*Keeping matched as rural counties with no cbsa 
merge 1:1 fipscnty using `population', keep(matched) nogen

rename state fips 

collapse (sum) popestim*, by(cbsa cbsatitle)
reshape long popestimate, i(cbsa) j(year)
rename popestimate cbsa_population

tempfile pop
save "`pop'"

***************************************************************************
*Cleaning traffic data here to then use in part B

use "$dropbox/Data/traffic_data.dta", clear

*First step w traffic data is to get top three sites by volume per cbsa
bys sitenum: egen total_vol = sum(traffic_volume)

*Negative here to get descending order to then keep top three
gen negtotal_vol = -total_vol

bys cbsa (negtotal_vol): gen ranked_volume = sum(negtotal_vol != negtotal_vol[_n-1])

keep if ranked_volume < 4
drop negtotal_vol ranked_volume

*state fips
replace fips = floor(fips/1000)

* State for Boston is NH? Means traffic sensor in NH?? Similarly, Chicago = Indiana, 
* but that one makes a bit more sense 

*population weight by CBSA? 

tempfile traffic
save "`traffic'"

*Cleaning hurricane waiver data and tempfile it to merge later

use "$dropbox/Data/gas_prices_waiver_panel_w_weather.dta", clear

drop _merge

*There are lots of lag variables in above dataset. Could be used later

keep date waiver ypc waiver_rfg waiver_rvp rvp_active treat_rfg mon_after_rvp ///
	 two_mon_after_rvp mon_after_rfg two_mon_after_rfg state city cbsa gas_p
rename gas_p gas_p_waiver
tostring cbsa, replace

tempfile h_waivers
save "`h_waivers'"


********************************************************************************

******************************************************************************

use "$dropbox/Data/Gasoline Prices 1998-2012.dta", clear

tempfile gas_price

*Data at CBSA level. Going to take state average to then merge to state level data
*Just want this for Li et al replication. Use CBSA level for our later analysis
bys year month fips: egen mean_gas_p = mean(gas_p)
bys year month fips: egen mean_lngas_p = mean(gas_p)

keep year month mean_gas_p mean_lngas_p fips 
rename mean_gas_p gas_p
rename mean_lngas_p lngas_p

quietly bys year month fips: gen dup = cond(_n==1, 0, _n)

drop if dup>1
drop dup

save "`gas_price'"

use "$dropbox/Data/Gasoline Tax & Consumption 1998-2012.dta", clear

merge 1:1 fips year month using `gas_price', nogen

merge m:1 fips year month using "$dropbox/Data/BLS Unemployment 1976 - 2018.dta", ///
		  keep(master match) nogen


/*	Li et al reg specification :
	Regress gasoline consumption per person by state and year on p, the tax 
	exclusive gasoline price, t the total state and federal gasoline tax, 
	all log specifications! */
	
*Need population data (do we?) to then get per capita measures for the following reg

gen lngas = ln(gasoline)
gen ln_nomtax = ln(gas_tax_nom)
gen interact_pt = ln(1 + gas_tax/gas_p)

label var lngas "log gasoline price"
label var interact_pt "log (1 + tax ratio)"

tempfile tmpfile

reg lngas lngas_p interact_pt i.year i.fips

regsave using "`tmpfile'"

preserve

use "`tmpfile'", clear
keep if _n ==1 | _n ==2
replace var = subinstr(var, "lngas_p", "log gasoline price", 1)
replace var = subinstr(var, "interact_pt", "log (1 + tax ratio)", 1)

local fn "A */** next to coefficient indicates significance at the 5/1% level."

texsave using "$dropbox/Replication Tables/Li_1.tex", title(Li et al Elasticity OLS) ///
autonumber footnote("`fn'") varlabels replace

restore

/*	Coglianese et al reg specification :
	Event study graphs 
	reg month-to-month change in log gasoline puirhcases on event study indicator
	variables for the month before, during and after gasoline tax*/

* Lag variables for months around tax
gen month_before_tax = 0
gen month_of_tax = 0
gen month_after_tax = 0

sort fips year month

* For tax indicator variables (assuming increase)*
bys fips: replace month_of_tax = 1 if gas_tax_nom[_n] > gas_tax_nom[_n-1] & year !=1998
*For negative tax (ie nominal tax decreases)
*bys fips: replace month_of_tax = 1 if gas_tax_nom[_n] < gas_tax_nom[_n-1] & year !=1998

bys fips: replace month_before_tax = 1 if month_of_tax[_n+1] == 1
bys fips: replace month_after_tax = 1 if month_of_tax[_n-1] == 1

* To produce the event study graphs need to get the change in each of these vars
* by month
bys fips: gen delta_lntax = ln_nomtax[_n] - ln_nomtax[_n-1]
bys fips: gen delta_lnprice = lngas_p[_n] - lngas_p[_n-1]
bys fips: gen delta_lngas = lngas[_n-1] - lngas[_n]
by fips: gen delta_unemployment = unemployment[_n-1] - unemployment[_n]
by fips: gen delta_lninteract = interact_pt[_n] - interact_pt[_n-1]

* Getting tax years to then produce residuals from regression 
* From Coglianese et al, time normalized relative to the month of change (t=0 here)
* Regress month to month chnages in the dependent variable on changes in state level
* unemployment rate and set of month of sample FE

gen tax_in_year = 0
bys fips year (month_of_tax): replace tax_in_year = 1 if month_of_tax[_N] == 1 

bys fips year (month_of_tax): gen norm_month = month - month[_N] ///
							  if tax_in_year
							  
* Tax residuals
reg delta_lntax delta_unemployment i.month
predict resid_lntax, residuals


* Price	resids						  
reg delta_lnprice delta_unemployment i.month
predict resid_lnprice, residuals 						  

* Consumption resids						  
reg delta_lngas delta_unemployment i.month
predict resid_lngas, residuals


bys norm_month: egen mean_resid_tax = mean(resid_lntax)	
bys norm_month: egen mean_resid_price = mean(resid_lnprice)						  							  
bys norm_month: egen mean_resid_gas = mean(resid_lngas)


graph twoway (line mean_resid_tax norm_month if norm_month > - 6, lpattern(dash)) ///
(line mean_resid_price norm_month if norm_month > - 6) ///
(line mean_resid_gas norm_month if norm_month > - 6, lpattern(shortdash_dot)), ///
title("Event Study (tax @ t=0)") ///
legend(order(1 2 3) label(1 "Tax") label(2 "Gas Price") label(3 "Gas Consumption")) ///
xline(0) xtitle("") ytitle("Change in Log")

graph export "$dropbox/Replication Tables/Coglianese_EventStudyRep.png", replace

sort fips year month							 

*OLS
reg delta_lngas month_of_tax month_before_tax month_after_tax i.month


*IV Regs (Davis & Killian style)

*ivreg2 delta_lngas delta_unemployment l.delta_lnprice ///
*			   f.delta_lnprice (delta_lnprice = month_of_tax i.month)
			   
********************************************************************************			   
			   
/* Comparisons: I’d like to compare, over a common time period and a 
	quasi-common set of locations, the elasticities produced with the 
	monthly consumption data from DOT (fuel tax records) and our measure 
	of gasoline prices and the elasticities produced with our traffic 
	count data. We may want to include some measure of CBSA population 
	to scale the traffic count data, just as the monthly consumption 
	data from DOT fuel tax records have been scaled by state population. 
	In this case, the quasi-common set of locations would mean using a set 
	of states in the former analysis that maps to the states implicitly 
	included in our city-level analyses.
*/

********************************************************************************
*Save prices and futures to tempfiles to then merge 

use "$dropbox/Data/Gasoline Prices 1998-2012.dta", clear
tempfile prices
save "`prices'"

use "$dropbox/Data/Gasoline Futures 1998-2017.dta", clear
tempfile futures
save "`futures'"

use "$dropbox/Data/Gasoline Tax & Consumption 1998-2012.dta", clear

merge 1:m fips year month using `prices', nogen keep(match)
merge m:1 year month using `futures', nogen keep(match)

merge m:1 fips year month using "$dropbox/Data/BLS Unemployment 1976 - 2018.dta", ///
		  keep(master match) nogen

drop date
		  
tempfile prices_upd
save "`prices_upd'"
		  
use "`traffic'", clear 

merge m:1 cbsa year month using "`prices_upd'", keep (master matched) nogen


tostring cbsa, replace
merge m:1 cbsa year using `pop', keep(matched) nogen

*Per Capita traffic measure here (VMT/popn[CBSA])

gen traffic_pc = traffic/cbsa_population*1000
gen log_traffic_pc = log(traffic_pc)
gen ln_nomtax = log(gas_tax_nom)
gen ln_gascons = log(gasoline)

*We want to include day of week effects (fixed?) see how others do this
gen dayofweek = dow(date)
forv i=0/6 {
	local w: word `=`i'+1' of `c(Weekdays)'
	local wd "`wd' `i' `w'"
}
label def dow `wd'
label val dayofweek dow


*Run through regs here

forval i=0/1{
	if `i'==0{
		bys sitenum (date): gen delta_lntax = ln_nomtax[_n] - ln_nomtax[_n-1]
		bys sitenum (date): gen delta_lnprice = lngas_p[_n] - lngas_p[_n-1]
		bys sitenum (date): gen delta_lngascons = ln_gascons[_n] - ln_gascons[_n-1]
		bys sitenum (date): gen delta_gascons = gasoline[_n] - gasoline[_n-1]
		bys sitenum (date): gen delta_gas_p = gas_p[_n] - gas_p[_n-1]
		by sitenum (date): gen delta_unemployment = unemployment_rate[_n] - unemployment_rate[_n-1]

		encode cbsa, gen(newcbsa)

		*First compute price/quantity elasticities from tax data 
		
		*This specification I think biases towards those days we have traffic data
		*if not systematic in any way (ie. traffic counter collection days are 
		*random) shouldn't be a problem
		*reg delta_lngascons delta_lnprice delta_unemployment i.month i.year i.newcbsa
		reg ln_gascons lngas_p unemployment_rate i.dayofweek i.month i.year i.newcbsa
		
		tempfile reggas
		save "`reggas'"
		
		margins, eyex(lngas_p) atmeans nose
		
		
		****Traffic count elasticity******

		*First going to calculate elasticity at the daily level
		bys sitenum (date): gen delta_logtraffic_daily = log_traffic_volume[_n] - log_traffic_volume[_n-1]
		replace delta_logtraffic_daily = 0 if delta_logtraffic_daily ==.

		*Per capita traffic count 
		bys sitenum (date): gen delta_logtraffic_pc= log_traffic_pc[_n] - log_traffic_pc[_n-1]
		replace delta_logtraffic_pc = 0 if delta_logtraffic_pc ==.

		*reg delta_logtraffic_daily delta_lnprice delta_unemployment i.month i.year i.newcbsa 
		reg log_traffic_volume lngas_p i.dayofweek i.month i.year i.newcbsa
		/*
		summarize delta_lnprice
		local meanp = r(mean)
		local f = `meanp'*_b[delta_lnprice] + _b[_cons]

		display "price :eyex = " (`meanp'*_b[delta_lnprice])/`f'
		*/
		*margins, eyex(delta_lnprice) atmeans nose
		margins, eyex(lngas_p) atmeans nose
		
		*Checking out per capita elasticity to check if different from unadjusted 
		*reg delta_logtraffic_pc delta_lnprice delta_unemployment i.month i.year i.newcbsa

		*margins, eyex(delta_lnprice) atmeans nose
		
	}
	else if `i'==1{
		

		**Monthly level calculation (change forval i=0/0 to 0/1)

		preserve

		collapse log_traffic_volume ln_gascons lngas_p unemployment_rate, by(month year sitenum)

		bys sitenum (month year): gen delta_logtraffic =  log_traffic_volume[_n] - log_traffic_volume[_n-1]
		bys sitenum (month year): gen delta_lngas_p = lngas_p[_n] - lngas_p[_n-1]
		bys sitenum (month year): gen delta_ln_gascons = ln_gascons[_n] - ln_gascons[_n-1]
		bys sitenum (month year): gen delta_unemployment = unemployment_rate[_n] - unemployment_rate[_n-1]

		replace delta_logtraffic = 0 if delta_logtraffic ==.
		replace delta_lngas_p = 0 if delta_lngas_p ==.
		replace delta_ln_gascons = 0 if delta_ln_gascons ==.
		
		reg delta_ln_gascons delta_lngas_p i.month i.year 
		*reg ln_gascons lngas_p i.month i.year 
		
		margins, eyex(delta_lngas_p) atmeans nose

		reg delta_logtraffic delta_lngas_p i.month i.year 
		*reg delta_logtraffic delta_lngas_p delta_lninteract i.month
		*reg log_traffic_volume lngas_p

		margins, eyex (delta_lngas_p) atmeans nose
		
		
		*Now going to test if same coeff from cons regression above with 
		*log traffic monthly 
		gen byte sample = 1
		
		append using "`reggas'"
		replace sample = 2 if missing(sample)
		
		forval i==1/12{
			foreach ind in "month" "year"{
				gen `ind'_`i' = 0
				replace `ind'_`i' = 1 if `ind' == `i'
				replace `ind'_`i' = 1 if `ind' == `i' + 2000
			}
		}
			
		reg delta_lngascons delta_lnprice delta_unemployment month* year* if sample == 2
		est store reggas
		
		reg delta_logtraffic delta_lngas_p unemployment_rate month* year* if sample == 1
		est store trafmonth
		
		suest reggas trafmonth
		
		test [reggas_mean]delta_lnprice = [trafmonth_mean]delta_lngas_p

		restore
	}
}

********************************************************************************
/* Extensions for various shocks: 

	I’d like to allow the price elasticity to vary with whether a hurricane 
	shock has occurred that triggers a waiver of fuel content regulations, 
	the start of the RFG season (May 15), as well as include the price without
	interaction in the models. We may also break some of these out into event 
	studies.
	
	For these hurricane shocks and for RFG phase-ins, let’s include an 
	indicator for 1-month after the event. And let’s focus our work here 
	on the daily panel. We should also include some distributed lags of 
	prices on consumption in the daily model.
*/

********************************************************************************

rename OPIS_city city

*Going to use twp merges first city then cbsa waiver dates to traffic data

merge m:1 date city using `h_waivers', keep(master matched) nogen

merge m:1 date cbsa using `h_waivers', keep(master matched) nogen

*ivreg2 delta_logtraffic (delta_lngas_p = delta_lninteract) i.month

gen ln_gas_p_waiver = log(gas_p_waiver)
gen ln_gastax = log(gas_tax)


encode cbsa, gen(cbsacode)
drop cbsa
rename cbsacode cbsa

encode sitenum, gen(sitecode)
encode city, gen(citycode)

*Need to figure out the panel var here - Austin/Tampa, Denver/SLC, Columbus/Houston
*share a site number?
*Going to replace Tampa with a couple 1's infront of it to give it a unique ID

replace sitecode = 1129 if sitecode == 29 & citycode == 44
replace sitecode = 1147 if sitecode == 47 & citycode == 38
replace sitecode = 1121 if sitecode == 21 & citycode == 11

*Dropping these for now
drop if gas_p_waiver ==.

gen int week = wofd(date)

bys week sitenum: egen average_weekly_price_waiver = mean(ln_gas_p_waiver)
bys week sitenum: egen average_weekly_price = mean(lngas_p)

tsset sitecode date

*Because of all the gaps in the data need to generate the lags manually
gen p_waiver_lag1 = ln_gas_p_waiver[_n-1]
gen p_waiver_lag2 = ln_gas_p_waiver[_n-2]
gen p_waiver_lag3 = ln_gas_p_waiver[_n-3]
					   
*Checking the sum of the coefficient estimates and testing statistical significance
*for different price/lag variables 
/*
lincom ln_gas_p_waiver+L.ln_gas_p_waiver

lincom L.ln_gas_p_waiver+L2.ln_gas_p_waiver

lincom L2.ln_gas_p_waiver+L3.ln_gas_p_waiver

lincom L.ln_gas_p_waiver+L3.ln_gas_p_waiver

lincom ln_gas_p_waiver+L3.ln_gas_p_waiver
*/

*Saving regs to then write to a latex file here 
*Uses regsave (documentation online)
tempfile tmpfile1
local replace "replace"

rename log_traffic_volume ln_traffic
						  
foreach reg in "reg" "ivreg2"{
	foreach depvar of varlist ln_traf ln_gascons{
		foreach regressor of varlist waiver ln_gastax{
			if "`reg'" == "ivreg2"{
				`reg' `depvar' (ln_gas_p_waiver = `regressor') ///
						mon_after_rfg unemployment_rate i.year i.month i.cbsa ///
						c.month#i.cbsa, cl(sitenum)
					
				regsave using "`tmpfile1'", addlabel(IV,"`regressor'",depvar, "`depvar'") ///
				table(`reg'_`depvar'_`regressor', format(%12.2f) asterisk(5 1) parentheses(stderr)) ///
				`replace'
		
				local replace "append"
			}			
			else {			
				`reg' `depvar' ln_gas_p_waiver `regressor' ///
				mon_after_rfg unemployment_rate i.year i.month i.cbsa c.month#i.cbsa, cl(sitenum)
				
				regsave using "`tmpfile1'", addlabel(depvar, "`depvar'", regressor, "`regressor'") ///
				table(`reg'_`depvar'_`regressor', format(%12.2f) asterisk(5 1) parentheses(stderr)) ///
				`replace'
		
				local replace "append"
			}
	}
 }
}

preserve

use "`tmpfile1'", clear
replace var = subinstr(var,"_coef","",1)
replace var = "" if strpos(var,"stderr")!=0

*Tagging SE to then drop from data 
replace var = "monthse" if strpos(var[_n+1], "month")!=0
replace var = "cbsase" if strpos(var[_n+1], "cbsa")!=0
replace var = "yearse" if strpos(var[_n+1], "year")!=0
replace var = "yearse" if strpos(var[_n+1], "cons")!=0

*Dropping fixed effects from Latex output here
drop if strpos(var, "month")!=0
drop if strpos(var, "cbsa")!=0
drop if strpos(var, "year")!=0

*Rename vars for display
replace var = subinstr(var, "ln_gas_p_waiver", "log price(waiver dataset)", 1)
replace var = subinstr(var, "unemployment_rate", "unemployment", 1)
replace var = subinstr(var, "ln_gas_tax", "log gas tax", 1)
label variable var "Variable name"

local fn "A */** next to coefficient indicates significance at the 5/1% level. SEs clustered at the traffic sensor level"

texsave using "$dropbox/Replication Tables/Traffic_Regs.tex", title(OLS(1-4) IV (5-8) Estim Log Traffic Volume) ///
autonumber footnote("`fn'") varlabels replace						  

restore						  

tempfile tmpfile2
local replace "replace"
						  
foreach reg in "reg" "ivreg2"{
	foreach depvar of varlist ln_traf ln_gascons{
		if "`reg'" == "ivreg2"{
			`reg' `depvar' (ln_gas_p_waiver = p_waiver_lag1 p_waiver_lag2 p_waiver_lag3) ///
					waiver ln_gastax mon_after_rfg unemployment_rate i.year i.month i.cbsa ///
					c.month#i.cbsa, cl(sitenum)
					
			regsave using "`tmpfile2'", addlabel(IV,"p+1",depvar, "`depvar'") ///
			table(`reg'_`depvar', format(%12.2f) asterisk(5 1) parentheses(stderr)) ///
			`replace'
		
			local replace "append"
		}			
		else {			
			`reg' `depvar' ln_gas_p_waiver waiver ///
			mon_after_rfg unemployment_rate i.year i.month i.cbsa c.month#i.cbsa, cl(sitenum)
				
			regsave using "`tmpfile2'", addlabel(depvar, "`depvar'") ///
			table(`reg'_`depvar', format(%8.3f) asterisk(5 1) parentheses(stderr)) ///
			`replace'
		
			local replace "append"
		}
	
 }
}

preserve

use "`tmpfile2'", clear
replace var = subinstr(var,"_coef","",1)
replace var = "" if strpos(var,"stderr")!=0

*Tagging SE to then drop from data 
replace var = "monthse" if strpos(var[_n+1], "month")!=0
replace var = "cbsase" if strpos(var[_n+1], "cbsa")!=0
replace var = "yearse" if strpos(var[_n+1], "year")!=0
replace var = "yearse" if strpos(var[_n+1], "cons")!=0

*Dropping fixed effects from Latex output here
drop if strpos(var, "month")!=0
drop if strpos(var, "cbsa")!=0
drop if strpos(var, "year")!=0

*Rename vars for display
replace var = subinstr(var, "ln_gas_p_waiver", "log price(waiver dataset)", 1)
replace var = subinstr(var, "unemployment_rate", "unemployment", 1)
replace var = subinstr(var, "ln_gas_tax", "log gas tax", 1)
label variable var "Variable name"

local fn "A */** next to coefficient indicates significance at the 5/1% level. SEs clustered at the traffic sensor level"

texsave using "$dropbox/Replication Tables/Traffic_Regs_plag.tex", title(OLS(1-2) IV (3-4) Estim Log Traffic Volume) ///
autonumber footnote("`fn'") varlabels replace						  

restore		


rename  ln_traffic log_traffic_volume
* All of above significant w/o cluster, add cluster at cbsa and none significant
* what does this mean? 


* Try an event study graph for waivers here, details on set up earlier (243-)
preserve

*Want to keep it at daily level, get daily change and then generate event windows 
collapse (sum) ln_gas_p_waiver log_traffic_volume unemployment waiver mon_after_rfg, by(sitecode date)

gen month = month(date)
gen year = year(date)

bys sitecode: gen delta_waiver = waiver[_n] - waiver[_n-1]
bys sitecode: gen delta_lntraffic = log_traffic_volume[_n] - log_traffic_volume[_n-1]
bys sitecode: gen delta_lnprice = ln_gas_p_waiver[_n] - ln_gas_p_waiver[_n-1]
bys sitecode: gen delta_unemployment = unemployment_rate[_n] - unemployment_rate[_n-1]

*bys sitecode: gen delta_rfg_lag = mon_after_rfg[_n] - mon_after_rfg[_n-1]
*bys sitecode: gen delta_rfg_lag2 = two_mon_after_rfg[_n] - two_mon_after_rfg[_n-1]

bys sitecode month: gen datenum = _n
bys sitecode month: gen target=datenum if waiver == 1
egen td = min(target), by(sitecode month)
drop target
gen diff = datenum-td

* Waiver residuals
reg delta_waiver delta_unemployment i.month
predict resid_waiver, residuals


* Price	resids						  
reg delta_lnprice delta_unemployment i.month
predict resid_lnprice, residuals 						  


* Traffic resids						  
reg delta_lntraffic delta_unemployment i.month
predict resid_lntraffic, residuals


bys diff: egen mean_resid_waiver = mean(resid_waiver)	
bys diff: egen mean_resid_price = mean(resid_lnprice)						  							  
bys diff: egen mean_resid_traffic = mean(resid_lntraffic)

drop if diff > 10
drop if diff < -10

graph twoway (line mean_resid_waiver diff, lpattern(dash)) ///
(line mean_resid_price diff) ///
(line mean_resid_traffic diff, lpattern(shortdash_dot)), ///
title("Daily Changes(Waiver @ t=0)") ///
legend(order(1 2 3) label(1 "Waiver") label(2 "Gas Price") label(3 "Monthly Traffic (CBSA)")) ///
xline(0) xtitle("") ytitle("")

graph export "$dropbox/Replication Tables/Waiver_DailyEventStudy.png", replace
						  
restore
						  
save "$dropbox/Data/Fuel Elasticity Panel.dta", replace

