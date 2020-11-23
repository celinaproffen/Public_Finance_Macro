clear all

*-------------------------------------------------------------------------------
* exercise 1
*-------------------------------------------------------------------------------

* create quartely aggregated data from CEX
*-------------------------------------------------------------------------------

* change directory
cd C:\Users\david\Documents\GitHub\Public_Finance_Macro\PS2_davide

* import data
use data/raw/cexdata.dta, replace

* add label (courtesy of Celina)
*-------------------------------------------------------------------------------
label var id "HH identification number"
label var inumber "Interview number, values from 2 to 5"
label var month "Interview month"
label var quarter "Interview quarter"
label var year "Interview year"
label var hhsize "HH size, including children"
label var age "Age of the household head"
label var grinc "HH income in past 12 months before taxes deated by the CPI"
label var netinc "HH income in past 12 months, after taxes deated by the CPI"
label var food "HH Food cons. in past 3 months, deated by the price deator for food"
label var ndcons1 "Nondurable + apparel/health/educ cons. in past 3 months, deated"
label var ndcons2 "Nondurable cons. in past 3 months, deated"
label var ndconsserv "Nondurables + imputed service ows from consumer durables"
label var totcons "Total cons. in past three 3 months, deated by the CPI"

* summary statistic
summarize

* change year structure
replace year = 1900 + year
replace year = 2000 if year == 1900

* create date variable
gen date = quarterly(string(year) + " " + "Q" + string(quarter), "YQ")
format date %tq

label var date "Date"

* create average and median consumption variables
collapse (mean) totcons_a=totcons hhsize_a=hhsize age_a=age grinc_a=grinc netinc_a=netinc food_a=food ndcons1_a=ndcons1 ndcons2_a=ndcons2 ndconsserv_a=ndconsserv (median) totcons_m=totcons hhsize_m=hhsize age_m=age grinc_m=grinc netinc_m=netinc food_m=food ndcons1_m=ndcons1 ndcons2_m=ndcons2 ndconsserv_m=ndconsserv, by(date)

* save 
save data/cexdata_agg_q.dta, replace

* import quarterly GDP data
*-------------------------------------------------------------------------------

clear all

* import data
import delimited using data/raw/gdp_q.csv

* create NBER recession datahttps://www.nber.org/research/data/us-business-cycle-expansions-and-contractions
* source available here: 
gen nber_rec = cond(date == "1990-07-01" | date == "1990-10-01" | date == "1991-01-01" | date == "2000-01-01" | date == "2000-04-01" | date == "2000-07-01" | date == "2000-10-01", 1, 0)


* change date format
split date, p(-)
replace date2 = "Q1" if date2 == "01"
replace date2 = "Q2" if date2 == "04"
replace date2 = "Q3" if date2 == "07"
replace date2 = "Q4" if date2 == "10"
drop date3 date
gen date = quarterly(date1 + " " + date2, "YQ")
format date %tq

* add label
label var date "Date"
label var gdp "nominal GDP"
label var nber_rec "NBER recession"

* drop unnecessary variables
drop date1 date2

* save 
save data/gdp_q.dta, replace

* create plots
*-------------------------------------------------------------------------------

clear all

* change directory
cd C:\Users\david\Documents\GitHub\Public_Finance_Macro\PS2_davide

* import data
use data/cexdata_agg_q.dta

* merge dataframe
merge 1:1 date using data/gdp_q.dta 

* create % change data 
tsset date
gen l_netinc = ln(netinc_a)
gen g_netinc = D.l_netinc

gen l_totcons = ln(totcons_a)
gen g_totcons = D.l_totcons

gen l_gdp = ln(gdp)
gen g_gdp = D.l_gdp

replace nber_rec = 0.1*nber_rec
gen nber_rec_neg = -nber_rec

* plot consumption and income variation
twoway (line g_totcons date) (line g_netinc date) (bar nber_rec_neg date) (bar nber_rec date)
graph export graphs/inc_cons.png, replace

replace nber_rec = 0.6*nber_rec
replace nber_rec_neg = -nber_rec

* plot consumption and income variation
twoway (line g_netinc date) (line g_gdp date) (bar nber_rec_neg date) (bar nber_rec date)
graph export graphs/inc_gdp.png, replace


