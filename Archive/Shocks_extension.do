bys sitecode (month year): gen delta_logtraffic =  log_traffic_volume[_n] - log_traffic_volume[_n-1]
bys sitecode (month year): gen delta_ln_gascons = lngas[_n] - lngas[_n-1]
bys sitecode (month year): gen delta_unemployment = unemployment_rate[_n] - unemployment_rate[_n-1]
		
bys sitecode (month year): gen delta_lngas_p = log_gas_p[_n] - log_gas_p[_n-1]
bys sitecode (month year): gen delta_lngas_p2 = log_gas_p2[_n] - log_gas_p2[_n-1]
/*replace delta_logtraffic = 0 if delta_logtraffic ==.
replace delta_lngas_p = 0 if delta_lngas_p ==.
replace delta_ln_gascons = 0 if delta_ln_gascons ==.
*/
		
bys sitenum (date): gen delta_lntax = ln_nomtax[_n] - ln_nomtax[_n-1]
		bys sitenum (date): gen delta_lnprice = log_gas_p[_n] - log_gas_p[_n-1]
		bys sitenum (date): gen delta_lngascons = lngas[_n] - lngas[_n-1]
		bys sitenum (date): gen delta_gascons = gasoline[_n] - gasoline[_n-1]
		bys sitenum (date): gen delta_gas_p = gas_p[_n] - gas_p[_n-1]
		by sitenum (date): gen delta_unemployment = unemployment_rate[_n] - unemployment_rate[_n-1]
		
		****Traffic count elasticity******

		*First going to calculate elasticity at the daily level
		bys sitenum (date): gen delta_logtraffic_daily = log_traffic_volume[_n] - log_traffic_volume[_n-1]
		replace delta_logtraffic_daily = 0 if delta_logtraffic_daily ==.

		*Per capita traffic count 
		bys sitenum (date): gen delta_logtraffic_pc= log_traffic_pc[_n] - log_traffic_pc[_n-1]
		replace delta_logtraffic_pc = 0 if delta_logtraffic_pc ==.

		*reg delta_lngascons delta_lnprice delta_unemployment i.month i.year i.newcbsa		
		

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

*We want to start interacting waiver with price data, include the futures rate 
*somehow here as well
gen price_waiver = log_gas_p * waiver 

*ivreg2 delta_logtraffic (delta_lngas_p = delta_lninteract) i.month

tsset sitecode date

*Because of all the gaps in the data need to generate the lags manually
gen p_waiver_lag1 = log_gas_p[_n-1]
gen p_waiver_lag2 = log_gas_p[_n-2]
gen p_waiver_lag3 = log_gas_p[_n-3]
					   
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
rename log_traffic_volume ln_traffic


tempfile tmpfile1
local replace "replace"
						  
foreach reg in "reg" "ivreg2"{
	foreach depvar of varlist ln_traffic lngas{
		foreach regressor of varlist ln_nomtax{
			if "`reg'" == "ivreg2"{
				`reg' `depvar' (log_gas_p = `regressor') ///
						mon_after_rfg unemployment_rate i.year i.month i.cbsa ///
						c.month#i.cbsa, cl(sitenum)
					
				regsave using "`tmpfile1'", addlabel(IV,"`regressor'",depvar, "`depvar'") ///
				table(`reg'_`depvar'_`regressor', format(%12.2f) asterisk(5 1) parentheses(stderr)) ///
				`replace'
		
				local replace "append"
			}			
			else {			
				`reg' `depvar' log_gas_p `regressor' ///
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

