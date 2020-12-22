********************************************************************************
* MACRO FINANANCE PROJECT 4
********************************************************************************

clear all   
*set maxvar 32767 
*set mem 1000m   
set more off   

* Celina/ Work
gl data "C:\Users\proffen\OneDrive\Dokumente\Econ Classes\Frankfurt\PublicFinanceMacro\projects\project4"
gl output "C:\Users\proffen\OneDrive\Dokumente\Econ Classes\Frankfurt\PublicFinanceMacro\projects\project4\stata_output"

* Laptop
*gl data "C:\Users\celin\OneDrive\Dokumente\Econ Classes\Frankfurt\PublicFinanceMacro\projects\project4"
*gl output "C:\Users\celin\OneDrive\Dokumente\Econ Classes\Frankfurt\PublicFinanceMacro\projects\project4\stata_output"

capture log close 
log using "$output\project4_CP.log",replace
use "$data\PANEL6812_HH.dta",replace
save "$output\data_project4_clean.dta",replace


********************************************************************************
* Basica data cleaning
********************************************************************************

*1) Focus on household with household heads in age range 25 to 60 according to variable age_head
 keep if 25<=age_head<=60
 
*2) Deflate income with cpi(household pre-government earnings (income) and household post-government earnings (income))
*3) Define log income (call it ln_inc) 

forvalues i=1(1)7{
 gen deflated_yhh`i' = yhh`i'/cpi
 gen ln_yhh`i' = log(deflated_yhh`i')
 }
 
 save "$output\data_project4_clean.dta",replace
 
********************************************************************************
* FIRST STAGE ESTIMATION
********************************************************************************
 
*1) Choose a degree 3 (order 4, i.e., a cubic) polynomial in age, variable age_head (I choose all coefficicents to be 1)
 gen pol3_age = (age_head + age_head^2 + age_head^3)
 gen age_head_squared = age_head^2
 gen age_head_cubed = age_head^3
 
 
*2) Choose a family size (I choose family size of 4 and create a dummy variable for it)
 
*3) Dummy for married (married)

*4) Dummy for college and high-school but no college
 gen d_coll=(Educ>=16)
 gen d_hs=(12<=Educ<16)
 
*5) interactions of the dummies d_coll and d_hs with the age polynomial
 gen d_coll_int_age  = d_coll*pol3_age
 gen d_hs_int_age    = d_hs*pol3_age
  
 tsset persnr year
 
 forvalues i=1(1)7{
	 quietly xtreg ln_yhh`i' d_coll d_hs age* famsize married i.year , i(persnr) fe robust cluster(persnr)
	 
	 *predict values of yhh`i' and residuals
	 predict ln_yhh`i'_hat
	 gen ln_yhh`i'_tilde = ln_yhh`i' - ln_yhh`i'_hat
	 
	 *residual variance for yhh`i'
	 egen var_ln_yhh`i'_tilde = sd(ln_yhh`i'_tilde)
	 replace var_ln_yhh`i'_tilde = var_ln_yhh`i'_tilde^2
 }
 
 /* Note.. i don't understand why this does not work
	quietly xtreg ln_yhh`i' famsize married i.year d_coll##age* d_hs#age* ///
	 d_coll#age_head_squared* d_hs##age_head_squared* d_coll#age_head_cubed* ///
	 d_hs#age_head_cubed* age_head_cubed* , i(persnr) fe robust cluster(persnr)
 */
 
 ********************************************************************************
 * Interpret these variances of residual log incomes!
 
 forvalues i=1(1)7{
	 display "residual variance ln_yy`i'_tilde:" as result var_ln_yhh`i'_tilde[1]
}
 
********************************************************************************
*Regress predicted values on degree 3 polynomial and a dummy for married
* there are two lines: one for married, one not married: how can I change labels/ colors?
 
  forvalues i=1(1)7{
	 reg ln_yhh`i'_hat age_head age_head_squared age_head_cubed married
	 predict smoothed_yhh`i'_hat
	 *twoway scatter smoothed_yhh`i'_hat age_head, title(Scatter plot for `: variable label yhh`i'') legend(size(medsmall))
	 graph twoway (scatter smoothed_yhh`i'_hat age_head if married==0) (scatter smoothed_yhh`i'_hat age_head if married==1), title(Scatter plot for `: variable label yhh`i'') legend(label(1 not married) label(2 married)) 
	 graph export "$output\smoothed_yy`i'_hat.png", replace
}
 
 save "$output\data_project4_clean.dta",replace
 
 
********************************************************************************
* SECOND STAGE ESTIMATION - PRE GOVERNMENT INCOME yhh5
********************************************************************************

clear all  
use "$output\data_project4_clean.dta",replace

* Since the PSID was cast at an annual frequency only until 1997, we want to exclude
* all later years
drop if year > 1997

/* Note: we need to decide on which measure of income to use!
Let's use yy5 i.e. Pre-Government HH Labor Income + .5 Payroll Tax (PRE)
precisely, let's use the residual of the logged and deflated yy5 measure after deducting it's predicted values
*/

* Define income measure we want to use!
gen res_ln_inc = ln_yhh5_tilde

* second stage regression of residual on its lag
sort persnr year
xtset persnr year 
rename res_ln_inc rsd
su rsd

* Generate lead variables:
gen rsd1 = F.rsd 
gen rsd2 = F2.rsd 

* Generate empirical moments.
* Note that in the sample we will consider now, Age = age_head in all cases
gen m00 = rsd^2 
gen m10 = rsd1^2
gen m01 = rsd*rsd1 
gen m02 = rsd*rsd2 
* The next command gives you the mean m00-m02 measures for each age-year pair
collapse (mean) m00 m10 m01 m02 [w=fwgt], by(Age year) 


***************************************************
* Carry out GMM estimation.

	gmm ///
		(mj0: m00 - {rho}^(2*(Age-$MinAge+1))*{var_z} - {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2) - {var_epsi} ) ///
		(mjp10: m10 - {rho}^(2*(Age-$MinAge+2))*{var_z} - {var_nu}*(1-{rho}^(2*(Age-$MinAge+2)))/(1-{rho}^2) - {var_epsi} ) ///
		(mj1: m01 - {rho}*({rho}^(2*(Age-$MinAge+1))*{var_z} + {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2)) ) ///
		(mj2: m02 - {rho}^2*({rho}^(2*(Age-$MinAge+1))*{var_z} + {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2)) ), ///
		instruments(Age, noconstant) instruments(mj0: m00, noconstant) ///
		instruments(mjp10: m10, noconstant) ///
		instruments(mj1: m01, noconstant) ///
		instruments(mj2: m02, noconstant) ///
		winitial(identity) from(rho 0.9 var_z 0.3 var_nu 0.03 var_epsi 0.03) onestep
	mat rt = e(b) 
	

********************************************************************************
* SECOND STAGE ESTIMATION - POST GOVERNMENT INCOME yhh6
********************************************************************************

clear all  
use "$output\data_project4_clean.dta",replace

* Since the PSID was cast at an annual frequency only until 1997, we want to exclude
* all later years
drop if year > 1997

/* Note: we need to decide on which measure of income to use!
Let's use yy6 i.e. Pre-Government HH Labor Income + .5 Payroll Tax +Transfers - Taxes (POST)
precisely, let's use the residual of the logged and deflated yy6 measure after deducting it's predicted values
*/

* Define income measure we want to use!
gen res_ln_inc = ln_yhh6_tilde

* second stage regression of residual on its lag
sort persnr year
xtset persnr year 
rename res_ln_inc rsd
su rsd

* Generate lead variables:
gen rsd1 = F.rsd 
gen rsd2 = F2.rsd 

* Generate empirical moments.
* Note that in the sample we will consider now, Age = age_head in all cases
gen m00 = rsd^2 
gen m10 = rsd1^2
gen m01 = rsd*rsd1 
gen m02 = rsd*rsd2 
* The next command gives you the mean m00-m02 measures for each age-year pair
collapse (mean) m00 m10 m01 m02 [w=fwgt], by(Age year) 


***************************************************
* Carry out GMM estimation.

	gmm ///
		(mj0: m00 - {rho}^(2*(Age-$MinAge+1))*{var_z} - {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2) - {var_epsi} ) ///
		(mjp10: m10 - {rho}^(2*(Age-$MinAge+2))*{var_z} - {var_nu}*(1-{rho}^(2*(Age-$MinAge+2)))/(1-{rho}^2) - {var_epsi} ) ///
		(mj1: m01 - {rho}*({rho}^(2*(Age-$MinAge+1))*{var_z} + {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2)) ) ///
		(mj2: m02 - {rho}^2*({rho}^(2*(Age-$MinAge+1))*{var_z} + {var_nu}*(1-{rho}^(2*(Age-$MinAge+1)))/(1-{rho}^2)) ), ///
		instruments(Age, noconstant) instruments(mj0: m00, noconstant) ///
		instruments(mjp10: m10, noconstant) ///
		instruments(mj1: m01, noconstant) ///
		instruments(mj2: m02, noconstant) ///
		winitial(identity) from(rho 0.9 var_z 0.3 var_nu 0.03 var_epsi 0.03) onestep
	mat rt = e(b) 

