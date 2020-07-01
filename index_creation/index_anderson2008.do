********************************************************************************
*** ANDERSON INDEX CREATION ***
********************************************************************************
/*
notes: inputs for the do file include the following locals:
	variable list: 		`index_varlist' (e.g. local varlist "var1 var2")
	index name:			`index_name' 	(e.g. local index_name "myindex")
	sample restriction: `index_sample' 	(e.g. local index_sample "attrit == 0")
	treatment variable:	`index_treat'   (e.g. local index_treat "")
	normalise to ctrl   `index_norm'	(e.g. local `index_norm' "1"; make non-missing)
	ctrl group must be  `index_treat' == 0

If the normalist to ctrl local is not missing, then demeaning and weights calculation
will be limited to the control group only
*/

** Replace If Missing Locals
if  missing("`index_sample'") 	local index_sample 	"1 == 1"					// sample identifier
if  missing("`index_treat'") 	local index_treat  	"treatment"					// treatment variable
if  missing("`index_name'") 	local index_name   	"index_name"				// variable name
if !missing("`index_norm'")		local index_norm	" & `index_treat' == 0 "	// normalise to control group always
else							local index_norm	"1 == 1"

** Check if Sample Exists
count if `index_sample'
if r(N) > 0 {	// if exists, continue
	** Rename Variables for Index
	local i = 0
	foreach x of local index_varlist {
		local i = `i' + 1
		gen temp_`i' = `x'
	}
	local nvars = `i'

	** Standardize Variables with respect to Control Group
	local index_var_missing ""
	forvalues count = 1/`nvars' {
		su temp_`count'  		if `index_sample' & `index_treat' == 0
		cap local sdev  = r(sd)
		su temp_`count' 		if `index_sample' &	`index_norm'
		cap local mean  = r(mean)
		if r(N) > 0	gen temp_`count'_z  = (temp_`count' - `mean') / `sdev'
		else local index_var_missing `index_var_missing' `count' " " // check if entire variable is missing
	}

	** Check if All Variables Exist
	if missing("`index_var_missing'") {
		** Create Anderson's Covariance Matrix of Variables
		matrix cov = I(`nvars')		// matrix defined in Appendix A
		forvalues i = 1/`nvars' {
			forvalues j = 1/`nvars' {
				if `i' >= `j' {
					egen   temp_cov_`i'`j' = sum(temp_`i'_z * temp_`j'_z) if `index_norm'
					su 	   temp_cov_`i'`j'
					matrix cov[`i',`j'] = r(mean)
					matrix cov[`j',`i'] = r(mean)
				}
			}
		}
		matrix invcov = syminv(cov)
		matrix unity  = J(1, rowsof(invcov), 1)
		matrix weights = unity * invcov	// simple column sum

		** Calculate Weights and Weighted Sum
		svmat weights, names(temp_weight_)
		forvalues count = 1/`nvars' {
			gen temp2_`count'  = temp_`count'_z * temp_weight_`count'[1] if `index_sample'
			gen temp2_weight_`count' = temp_weight_`count'[1] if !missing(temp_`count'_z)
		}

		** Calculate Temporary Index
		egen 	temp_index 			= rowtotal(temp2_* )		// weighted sum from Appendix A
		egen 	temp_index_weight 	= rowtotal(temp2_weight_*)	// W_ij from Appendix A
		replace temp_index			= temp_index / temp_index_weight

		** Normalise Index?
		*su 		temp_index if ( `index_treat' == 0 )
		*replace temp_index = temp_index / r(sd)

		** Create or Replace Index
		cap gen 	`index_name' = temp_index	if `index_sample'
		cap replace `index_name' = temp_index	if `index_sample'
	}

	** Drop Temporary Variables
	cap drop temp_*
	cap drop temp2_*
	cap drop temp_cov_*
	cap drop temp_weight_*
	cap drop temp_index*

	** Warning Message if Missing Variables
	if !missing("`index_var_missing'") di "error: variables `index_var_missing'missing"

}
