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
rename fips state_fips

save "$dropbox/Data/BLS Unemployment 1976 - 2018", replace

*/

* Census Population data here 
/*
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
*/
***************************************************************************

*Save prices and futures to tempfiles to then merge 

use "$dropbox/Data/Gasoline Futures 1998-2017.dta", clear
tempfile futures
save "`futures'"

*Pull monthly consumption from this dataset and then merge to traffic data for 
*replication purposes
use "$dropbox/Data/state_fuel_tax_consumption_1980_2016.dta", clear

keep if year > 1997
keep gasoline gas_tax year month state

statastates, name(state)
drop _merge state
rename state_abbrev state

tempfile monthly_consump 
save "`monthly_consump'"

********************************************************************************
*Construct panel here, merging in some external data sources not in traffic/waiver

use "$dropbox/Data/Gasoline Prices 1998-2012.dta", clear

rename fips state_fips

merge m:1 state_fips year month using "$dropbox/Data/BLS Unemployment 1976 - 2018.dta", ///
		  keep(master match) nogen

merge m:1 state_fips year month using `monthly_consump', keep(master matched) nogen



/*	Li et al reg specification :
	Regress gasoline consumption per person by state and year on p, the tax 
	exclusive gasoline price, t the total state and federal gasoline tax, 
	all log specifications! */
	

gen lngas = ln(gasoline)
gen log_gas_p = ln(gas_p)
gen ln_nomtax = ln(gas_tax)
gen interact_pt = ln(1 + gas_tax/gas_p)

label var lngas "log gasoline consumption"
label var interact_pt "log (1 + tax ratio)"

tempfile tmpfile

reg lngas log_gas_p interact_pt i.year i.month i.state_fips i.cbsa

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
	
*preserve

**** Removing Traffic for now (too many gaps in the data when merged to traffic)

*Collapsing to region (can change geographic unit here)
foreach region in state_fips{
	
	collapse gas_tax ln_nomtax log_gas_p lngas unemployment ///
			 interact_pt, by(`region' month year)
		
	* Lag variables for months around tax
	gen month_before_tax = 0
	gen month_of_tax = 0
	gen month_after_tax = 0
	
	gen date = ym(year, month)

	* For tax indicator variables (assuming increase)*
	bys `region' (date): replace month_of_tax = 1 if gas_tax > gas_tax[_n-1] & year !=1998
	*For negative tax (ie nominal tax decreases)
	*bys `region' (date): replace month_of_tax = 1 if gas_tax_nom[_n] < gas_tax_nom[_n-1] & year !=1998

	bys `region' (date): replace month_before_tax = 1 if month_of_tax[_n+1] == 1
	bys `region' (date): replace month_after_tax = 1 if month_of_tax[_n-1] == 1

	* To produce the event study graphs need to get the change in each of these vars
	* by month
	bys `region' (date): gen delta_lntax = ln_nomtax - ln_nomtax[_n-1]
	bys `region' (date): gen delta_lnprice = log_gas_p - log_gas_p[_n-1]
	bys `region' (date): gen delta_lngas = lngas[_n-1] - lngas
	*bys `region' (date): gen delta_traffic = log_traffic_volume[_n-1] - log_traffic_volume
	bys `region' (date): gen delta_unemployment = unemployment - unemployment[_n-1]
	bys `region' (date): gen delta_lninteract = interact_pt - interact_pt[_n-1]

	/*
	
	 Getting tax years to then produce residuals from regression 
	
	 From Coglianese et al, time normalized relative to the month of change (t=0 here)
	
	* Overlapping events could cause problems for the graph, 
	  going to create bins of months for each region 
	  
	*/
	
	egen number_taxchanges = total(month_of_tax), by(`region')
	
	egen bin = cut(date), at(456(10)636)
	
	
	gen tax_in_year = 0
	bys `region': replace tax_in_year = 1 if month_of_tax[_N] == 1

	bys `region' bin (month_of_tax): gen norm_month = date - date[_N] if month_of_tax[_N] == 1
		  
	* Regress month to month changes in the dependent variable on changes in state level
	* unemployment rate and set of month of sample FE
								  
	* Tax residuals
	reg delta_lntax delta_unemployment i.month
	predict resid_lntax, residuals


	* Price	resids						  
	reg delta_lnprice delta_unemployment i.month
	predict resid_lnprice, residuals 						  

	* Consumption resids						  
	reg delta_lngas delta_unemployment i.month
	predict resid_lngas, residuals
	
	/*
	*Traffic resids
	reg delta_traffic delta_unemployment i.month
	predict resid_traffic, residuals
	*/

	bys norm_month: egen mean_resid_tax = mean(resid_lntax)	
	bys norm_month: egen mean_resid_price = mean(resid_lnprice)						  							  
	bys norm_month: egen mean_resid_gas = mean(resid_lngas)
	*bys norm_month: egen mean_resid_traffic = mean(resid_traffic)


	graph twoway (line mean_resid_tax norm_month if norm_month > - 6 & norm_month < 6, lpattern(dash)) ///
	(line mean_resid_price norm_month if norm_month > -6 & norm_month < 6) ///
	(line mean_resid_gas norm_month if norm_month > - 6 & norm_month < 6, lpattern(longdash_dot)), ///
	title("Event Study (tax @ t=0)") ///
	legend(order(1 2 3) label(1 "Tax") label(2 "Gas Price") label(3 "Gas Consumption")) ///
	xline(0) xtitle("") ytitle("Change in Log")

	graph export "$dropbox/Replication Tables/Coglianese_EventStudyRep.png", replace

	sort `region' year month							 

	*OLS
	reg delta_lngas month_of_tax month_before_tax month_after_tax i.month
	
	tempfile tmpfile1
	regsave using "`tmpfile1'"
	

	preserve

	use "`tmpfile1'", clear
	keep if _n == 1 | _n == 2 | _n== 3

	texsave using "$dropbox/Replication Tables/Month_of_Tax_Coglianese", title(Event Study Regressions(Coglianese et al Rep)) ///
	autonumber varlabels replace
	
	restore
	
	reg delta_lngas month_of_tax month_before_tax month_after_tax i.month
	
	/*
	gen date = ym(year, month)

	xtset `region' date
	gen lag_delta_lnprice = L.delta_lnprice
	gen for_delta_lnprice = F.delta_lnprice

	*IV Regs (Davis & Killian style)

	ivreg2 delta_lngas delta_unemployment lag_delta_lnprice ///
				   for_delta_lnprice (delta_lnprice = month_of_tax i.month)
	*/		   

	*restore		
}
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

*Construct panel here, merging in some external data sources not in traffic/waiver

use "$dropbox/Data/traffic_waiver_fuel_panel.dta", clear

merge m:1 state_fips year month using "$dropbox/Data/BLS Unemployment 1976 - 2018.dta", ///
		  keep(master match) nogen

merge m:1 state_fips year month using `monthly_consump', keep(master matched) nogen

merge m:1 year month using `futures', nogen keep(match)

gen lngas = ln(gasoline)
gen ln_nomtax = ln(gas_tax)

*This var appears in Li et al. 
gen interact_pt = ln(1 + gas_tax/gas_p)

********************************************************************************
*Checking out the two different consumption measures time series here
preserve 
	
collapse gas_tax log_gas_p lngas log_traffic_volume unemployment ///
			 interact_pt, by(month year)
gen date = ym(year, month)
format date %tm

gen gas2000 = lngas if year == 2000 & month == 1
gen traffic2000 = log_traffic_volume if year == 2000 & month == 1

foreach v of varlist gas2000 traffic2000{
	replace `v' = `v'[_n-1] if mi(`v')
}

* Indexing to 01/2000 here for basis of comparison
replace lngas = lngas/gas2000
replace log_traffic_volume = log_traffic_volume/traffic2000

sort date
twoway (line lngas date) (line log_traffic_volume date), legend(order(1 2) ///
label(1 "Mean Log Gas Cons") label(2 "Mean Log Traffic Volume")) note("nb: indexed to 01/2000") ///
saving("$dropbox/Replication Tables/Consumption Time Series", replace)
	
restore	
********************************************************************************
tostring cbsa, replace

*Per Capita traffic measure here (VMT/popn[CBSA])
gen traffic_pc = traffic/pop*1000
gen log_traffic_pc = log(traffic_pc)

*We want to include day of week effects (fixed?) see how others do this
gen dayofweek = dow(date)
forv i=0/6 {
	local w: word `=`i'+1' of `c(Weekdays)'
	local wd "`wd' `i' `w'"
}
label def dow `wd'
label val dayofweek dow

*May use padd to check for geographic differentiation based on waiver granted post
*hurricane [map here: https://www.eia.gov/petroleum/marketing/monthly/pdf/paddmap.pdf]

gen padd = 1
replace padd = 2 if inlist(state_fips, 17, 18, 19, 20, 21, 26, 27, 29, 31, 38, ///
									   39, 40, 46, 47, 55)
replace padd = 3 if inlist(state_fips, 1, 5, 22, 28, 35, 48)
replace padd = 4 if inlist(state_fips, 8, 16, 30, 49, 56)
replace padd = 5 if inlist(state_fips, 4, 6, 32, 41, 53)

label var padd "Petroleum Admin for Defense Districts"

rename OPIS_city city

encode cbsa, gen(cbsacode)
drop cbsa
rename cbsacode cbsa

encode sitenum, gen(sitecode)
encode city, gen(citycode)

*Need to figure out the panel var here - Austin/Tampa, Denver/SLC, Columbus/Houston
*share a site number?
*Going to replace duplicates with other site codes, don't know what to do about 
*Denver/SLC

bys sitecode date: gen dup = cond(_N==1, 0, _n)
replace sitecode = sitecode + 100000 if dup == 1
replace sitecode = sitecode + 200000 if dup == 2
drop dup

gen int week = wofd(date)

*Weighted average for the traffic volume measure to control for differing number
*of days reported

bys week sitecode: egen average_weekly_price = mean(log_gas_p)

bys week sitecode: gen days_in_week = _N
*bys month year sitecode: gen days_in_mo = _N

bys week sitecode: egen average_traffic = mean(log_traffic_volume)

*Create some weekly indicators 
gen weekday = 0
replace weekday = 1 if inlist(dayofweek, 1, 2, 3, 4, 5)


*Create some indicators for hurricane landfall and interactions of price and 
*waiver status here
gen hurricane_lf = 0
gen hurricane_lf_L = 0
gen hurricane_lf_L2 = 0
gen hurricane_month = 0
gen hurricane_month_L = 0

*From Wiki: dates of landfall of Katrina, Rita, Gustav, Ike, Isaac, Sandy, resp.
local hurricanes "29aug2005 21sep2005 01sep2008 13sep2008 29aug2012 27oct2012"
foreach h of local hurricanes{

	replace hurricane_lf = 1 if date == td(`h') 
	replace hurricane_lf_L = 1 if date == td(`h') + 1 
	replace hurricane_lf_L2 = 1 if date == td(`h') + 2
	
	*Indicators for months of hurricane and one after starting with landfall date
	replace hurricane_month = 1 if inrange(date, td(`h'), td(`h') + 30)

	replace hurricane_month_L= 1 if inrange(date, td(`h') + 30, td(`h') + 60)		
}


*Creating hurricane and price interaction vars

gen p_hurricane = log_gas_p*hurricane_month
gen p_hurricane_L = log_gas_p*hurricane_month_L

label var hurricane_lf "Hurricane Landfall date"
label var hurricane_lf_L "Hurricane Landfall date t+1"
label var hurricane_lf_L2 "Hurricane Landfall date t+2"
label var hurricane_month "Month of Hurricane"
label var hurricane_month_L "Month Lag Hurricane"
label var p_hurricane "Price Hurricane Interaction"
label var p_hurricane_L "Price Hurricane Interaction t+1"


*We want to compare expected vs unexpected price increases so first gen some vars
*for expected price increase

foreach waiv in rfg rvpI rvpII{
	gen p_waiver_`waiv'= 0
	replace p_waiver_`waiv' = log_gas_p*treat_`waiv'
	label var p_waiver_`waiv' "Price waiver*`waiv'"
}
*******************************************************************

* Regression controls here *

********************************************************************

local controls_wthr prcp snow tmax tmin 

*Replacing missing values for snow and precip
replace snow = prcp if prcp >=0 & tmax <=32
replace snow = 0 if tmin > 40

egen miss_prcp = mean(prcp), by(year cbsa)
replace prcp = miss_prcp if mi(prcp)

*Run through regs here

*If returning to changes in var for later variations see Shocks_extension.do

forval i=0/0{
	if `i'==0{

		tempfile tmpfile1
		local replace "replace"
								  
		foreach reg in "reg"{
			foreach depvar of varlist log_traffic_volume lngas{
				if "`reg'" == "ivreg2"{
					`reg' `depvar' (log_gas_p = `insert_regressor') ///
							mon_after_rfg i.year i.month i.dayofweek i.sitecode, ///
						    cl(sitenum)
							
					regsave using "`tmpfile1'", addlabel("`IV,"`insert_regressor'",depvar, "`depvar'",FE,"month, year, traffic sensor") ///
					table(`reg'_`depvar', format(%12.3f) asterisk(5 1) parentheses(stderr)) ///
					`replace'
				
					local replace "append"
				}			
				else {			
					`reg' `depvar' log_gas_p  p_waiver_rfg p_hurricane ///
					`controls_wthr' i.dayofweek i.year i.month i.cbsa, cl(sitenum)
					
					*Test equality of variables 
					
					test p_waiver_rfg = p_hurricane 
					
					local Fstat = r(F)
					
					if "`depvar'" == "log_traffic_volume" {
						regsave using "`tmpfile1'", addlabel(depvar, "`depvar'",FE,"month, year, day of week, traffic sensor", Test, "p_hurricane = p_waiver_rfg", Fstat,"`Fstat'") ///
						table(`reg'_`depvar', format(%12.3f) asterisk(5 1) parentheses(stderr)) ///
						`replace'
					}
					else {
						regsave using "`tmpfile1'", addlabel(depvar, "`depvar'",FE,"month, year, day of week", Test, "p_hurricane = p_waiver_rfg", Fstat,"`Fstat'") ///
						table(`reg'_`depvar', format(%12.3f) asterisk(5 1) parentheses(stderr)) ///
						`replace'
					}
				}
					local replace "append"
				
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
		replace var = "dayofweekse" if strpos(var[_n+1], "dayofweek")!=0

		*Dropping fixed effects from Latex output here
		drop if strpos(var, "month")!=0
		drop if strpos(var, "cbsa")!=0
		drop if strpos(var, "year")!=0
		drop if strpos(var, "week")!=0

		*Rename vars for display
		replace var = subinstr(var, "log_gas_p", "log price", 1)
		replace var = subinstr(var, "p_hurricane", "price * hurricane landfall", 1)
		replace var = subinstr(var, "p_waiver_rfg", "price * RFG waiver active", 1)

		local fn "Standard errors are clustered at traffic site number. A */** next to coefficient indicates significance at the 5/1% level. SEs clustered at the traffic sensor level"

		texsave using "$dropbox/Replication Tables/Traffic_Regs.tex", title(OLS Estimates of Traffic Volume & Gas Consumption) ///
		autonumber footnote("`fn'") varlabels replace						  

		restore	
		
		
		*reg lngas log_gas_p p_waiver_rfg p_hurricane `controls_wthr' i.dayofweek i.month i.year i.sitecode
		
		
		*reg delta_logtraffic_daily delta_lnprice delta_unemployment i.month i.year i.newcbsa 
		reg log_traffic_volume log_gas_p `controls_wthr' i.dayofweek i.month i.year i.sitecode
		
		reg log_traffic_volume log_gas_p p_waiver_rfg p_hurricane `controls_wthr' i.dayofweek i.month i.year i.sitecode
		/*
		summarize delta_lnprice
		local meanp = r(mean)
		local f = `meanp'*_b[delta_lnprice] + _b[_cons]

		display "price :eyex = " (`meanp'*_b[delta_lnprice])/`f'
		*/
		*margins, eyex(delta_lnprice) atmeans nose
		
		*Checking out per capita elasticity to check if different from unadjusted 
		*reg delta_logtraffic_pc delta_lnprice delta_unemployment i.month i.year i.newcbsa

		*margins, eyex(delta_lnprice) atmeans nose
	}
	
	else if `i'==1{
		

		**Monthly level calculation (change forval i=0/0 to 0/1)

		preserve

		collapse log_traffic_volume lngas log_gas_p unemployment_rate average_weekly_price ///
				snow prcp tmax tmin, by(week month year sitecode)
		
		reg lngas log_gas_p `controls_wthr' i.month i.year i.sitecode
		
		*Want to know what the above looks like centered at 0, so center p by 
		*subtracting the mean of p from each value 
		
		egen mean_log_gas_p = mean(log_gas_p)
		gen c_log_gas_p = log_gas_p - mean_log_gas_p
		gen c_log_gas_p2 = (c_log_gas_p)^2
		
		reg lngas c_log_gas_p `controls_wthr' i.month i.year i.sitecode
		
		*Correlation goes from 0.958 to -0.17 with centered 
		corr c_log_gas_p c_log_gas_p2

		reg log_traffic_volume log_gas_p `controls_wthr' i.month i.year i.sitecode
		*reg delta_logtraffic delta_lngas_p i.month i.year 
		*reg delta_logtraffic delta_lngas_p delta_lninteract i.month
		*reg log_traffic_volume lngas_p

		*margins, eyex (delta_lngas_p) atmeans nose
		
		/*
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
		*/

		restore
	}
}
