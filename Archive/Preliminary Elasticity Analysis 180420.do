/*

Author: Joe Aldy

Date: Apr 2018

Last Update: Apr 2018

Descripton: Merge city fuel price data and consumption data to and analyze */
/* ALL PRICE DATA CONVERTED TO 2015$ IN INITIAL CONSTRUCTION OF DATASETS */

clear
capture log close

global dropbox1  "C:/Users/`c(username)'/Dropbox/Fuel Taxes"
global dropbox2  "C:/Users/`c(username)'/Dropbox/Boutique Fuels/Stata/output/intermediate"
global dropbox3  "C:/Users/`c(username)'/Dropbox/Boutique Fuels/Stata/data/BEA Income"
global dropbox4  "C:/Users/`c(username)'/Dropbox/Boutique Fuels/Stata/output/final_panels"

set matsize 10000
log using "$dropbox1\Output\Preliminary Elasticity Analysis 180420", replace
/* CPI-U Deflator  */
cd "$dropbox3/"
import excel using "CPIAUCSL.xls", cellrange(A12:B846) clear
rename (A B) (date cpiu)
g year = year(date)
collapse cpiu, by(year)
g cpiu2015a = cpiu if year==2015
egen cpiu2015 = max(cpiu2015a)
drop cpiu2015a
drop if year == 2016 // INCOMPLETE
sort year
save "$dropbox1/Output/CPI_U", replace
clear

cd "$dropbox2/fuel"

/* GASOLINE FUTURES DATA */
/* CONSTRUCT MONTHLY MEAN AND STANDARD DEVIATION OF FRONT MONTH CONTRACT PRICES */
/* USE REGULAR GAS FUTURE CONTRACTS THROUGH END OF 2005; RBOB STARTING 2006 */
use "regular_gas_future_contracts.dta"
g year = year(date)
g month = month(date)
drop if year>2005
sort year 
merge year using "$dropbox1/Output/CPI_U"
drop _merge
drop if year < 1998
ren ny_gas_future_1 gas_fut
replace gas_fut = gas_fut*cpiu2015/cpiu
g gas_fut_sd = gas_fut
g lngas_fut = ln(gas_fut)
g lngas_fut_sd = ln(gas_fut)
collapse (mean) gas_fut (mean) lngas_fut (sd) gas_fut_sd (sd) lngas_fut_sd, by(year month)
sort year month
save "$dropbox1/Output/Gasoline Futures 1984-2017", replace
clear

use "rbob_future_contracts.dta"
g year = year(date)
g month = month(date)
drop if year<2006
sort year 
merge year using "$dropbox1/Output/CPI_U"
drop _merge
drop if year < 1998
ren ny_rbob_future_1 gas_fut
replace gas_fut = gas_fut*cpiu2015/cpiu
g gas_fut_sd = gas_fut
g lngas_fut = ln(gas_fut)
g lngas_fut_sd = ln(gas_fut)
collapse (mean) gas_fut (mean) lngas_fut (sd) gas_fut_sd (sd) lngas_fut_sd, by(year month)
sort year month
append using "$dropbox1/Output/Gasoline Futures 1984-2017"
drop if month==.
sort year month
save "$dropbox1/Output/Gasoline Futures 1998-2017", replace
clear

/* GASOLINE PRICE DATA */
/* CONSTRUCT MONTHLY MEAN AND STANDARD DEVIATION OF CITY GASOLINE PRICES */
use "OPIS_data_1998_2012.dta"
g state1 = substr(state,2,2)
drop state
ren state1 state
do "$dropbox1/code/State FIPS codes.do"
g year = year(date)
sort year
merge year using "$dropbox1/Output/CPI_U"
drop _merge
drop if year < 1998
replace gas_p = gas_p*cpiu2015/cpiu
g month = month(date)
g gas_p_sd = gas_p
g lngas_p = ln(gas_p)
g lngas_p_sd = ln(gas_p)
collapse (mean) gas_p (mean) lngas_p (sd) gas_p_sd (sd) lngas_p_sd (mean) fips (mean) date, by(year month cbsa)
sort year month cbsa
save "$dropbox1/Output/Gasoline Prices 1998-2012", replace
clear

/* GASOLINE TAX & CONSUMPTION DATA */
use "$dropbox1/Output/state_fuel_tax_consumption_1980_2016.dta"
do "$dropbox1/code/State FIPS codes.do"
sort year
merge year using "$dropbox1/Output/CPI_U"
drop _merge
drop if year < 1998 
drop diesel diesel_tax federal_diesel_tax
g gas_tax_nom = gas_tax
g federal_gas_tax_nom = federal_gas_tax 
replace gas_tax = gas_tax*cpiu2015/cpiu
replace federal_gas_tax = federal_gas_tax*cpiu2015/cpiu
drop if year>2012
sort year month fips
save "$dropbox1/Output/Gasoline Tax & Consumption 1998-2012", replace
clear

/* STATE INCOME & POP DATA */
use "$dropbox2/BEA state personal income 1969-2014.dta"
drop if year<1998 | year>2012
ren stfips fips
keep fips year y2015 ypc2015 pop
sort year fips
save "$dropbox1/Output/BEA state personal income 1998-2012.dta", replace
clear

/* REGULATIONS & WAIVERS DATA */
use "$dropbox4/gas_prices_waiver_panel.dta"
drop if gas_p~=.
ren treat_rfg rfg
ren treat_carb carb
g rvp = treat_rvpI + treat_rvpII
replace rvp = 1 if rvp>1 & rvp~=.
g onemonth = mon_after_rfg+mon_after_rvp
g twomonth = two_mon_after_rfg+two_mon_after_rvp
replace onemonth = 1 if onemonth>1 & onemonth~=.
replace twomonth = 1 if twomonth>1 & twomonth~=.
collapse (max) waiver (max) rfg (max) carb (max) rvp (max) onemonth (max) twomonth, by (year month cbsa)
save "$dropbox1/Output/Regs & Waivers 1998-2012.dta", replace
clear


/* MERGE DATASETS */
use "$dropbox1/Output/Gasoline Prices 1998-2012"
merge 1:1 year month cbsa using "$dropbox1/Output/Regs & Waivers 1998-2012.dta"
drop _merge
sort year month fips
merge m:1 year month fips using "$dropbox1/Output/Gasoline Tax & Consumption 1998-2012"
drop _merge
sort year fips
merge m:1 year fips using "$dropbox1/Output/BEA state personal income 1998-2012.dta"
drop _merge 
sort year month
merge year month using "$dropbox1/Output/Gasoline Prices 1998-2012"
drop _merge
sort year month
merge m:1 year month using "$dropbox1/Output/Gasoline Futures 1998-2017"
drop _merge

/* CLEAN UP OMITTEDS */

foreach v in waiver rfg rvp onemonth twomonth carb {
	replace `v' = 0 if `v'==.
	}

/* CONSTRUCT DEPENDENT VARIABLE */
g lncons = ln(gasoline/pop)

/* REGRESSIONS */
replace date = floor(date)
tsset cbsa date

ivreg2 lncons (lngas_p = waiver rfg carb rvp gas_tax) ypc2015 i.year i.month i.cbsa, cl(cbsa)
*reg lncons lngas_p ypc2015 i.year i.month#i.cbsa, cl(cbsa)
*reg lncons lngas_p L1.lngas_p ypc2015 i.year i.month#i.cbsa, cl(cbsa)
*reg lncons lngas_p lngas_p_sd ypc2015 i.year i.month#i.cbsa, cl(cbsa)

drop if cbsa==.
gen date_to_use = ym(year, month)

gen int week = wofd(date)

bys week cbsa: egen average_weekly_price = mean(lngas_p)

tsset cbsa date_to_use

ivreg2 lncons (lngas_p lngas_p_sd = L.lngas_p waiver rfg carb rvp onemonth twomonth lngas_fut lngas_fut_sd gas_tax) ypc2015 i.year i.month i.cbsa, cl(cbsa)

ivreg2 lncons (lngas_p lngas_p_sd = average_weekly_price waiver rfg carb rvp onemonth twomonth lngas_fut lngas_fut_sd gas_tax) ypc2015 c.month#i.cbsa i.month i.cbsa, cl(cbsa)











