********************************************************************************
*** ANDERSON INDEX CREATION ***
********************************************************************************
/*
notes: inputs for the do file include the following locals:
	variable list: 		`index_varlist' (e.g. local index_varlist 	var1 var2)
	index name:			`index_name' 	(e.g. local index_name 		my_index)
	sample restriction: `index_sample' 	(e.g. local index_sample 	condition == "true")
	treatment variable:	`index_treat'   (e.g. local index_treat 	treatment_variable)
	control group value:`index_ctrl' 	(e.g. local index_ctrl 		1)
*/

** Defaults
if  missing("`index_name'") 	local index_name   		index_name	// variable name
if  missing("`index_sample'") 	local index_sample 		1 == 1		// sample identifier
if  missing("`index_treat'") 	local index_treat  		treatment	// treatment variable
if  missing("`index_ctrl'")		local index_ctrl		0			// default

** Check if Sample Exists
count if `index_sample'
if r(N) > 0 {	// if exists, continue
	** Rename Variables for Index
	local i = 0
	foreach x of local index_varlist {
		local i = `i' + 1
		gen temp1_`i' = `x'		if `index_sample'		// limit to relevant sample
	}
	local nvars = `i'

	** Standardize Variables With Respect To Control Group
	local index_var_missing ""
	forvalues count = 1/`nvars' {
		qui su temp1_`count'  		if `index_treat' == `index_ctrl' // normalise by control group
		if r(N) > 0	gen temp1_`count'_z  = (temp1_`count' - r(mean)) / r(sd)
		else local index_var_missing `index_var_missing' `count' " " // check if entire variable is missing for control group
	}

	** Check if All Variables Exist
	if missing("`index_var_missing'") {
		** Create Anderson's Covariance Matrix of Variables
		matrix cov = I(`nvars')				// Sigma matrix defined in Appendix A
		forvalues i = 1/`nvars' {
			forvalues j = 1/`nvars' {
				if `i' >= `j' {
					** Temporary (Non-Sample-Adjusted Pairwise Covariance Matrix)
					egen   temp_cov_`i'`j' = sum(temp1_`i'_z * temp1_`j'_z) if  `index_treat' == `index_ctrl'
					qui su temp_cov_`i'`j'
					matrix cov[`i',`j'] = r(mean)
					matrix cov[`j',`i'] = r(mean)
					* Temporary (Sample-Adjusted Pairwise Covariance Matrix)
					*correl temp1_`i'_z temp1_`j'_z if `index_treat' == `index_ctrl', covariance
					*matrix cov[`i',`j'] = r(cov_12)
					*matrix cov[`j',`i'] = r(cov_12)
				}
			}
		}
		* Temporary (Standard Covariance Matrix)
		*correl temp1_*_z if  `index_treat' == `index_ctrl' , covariance
		*matrix cov = r(C)

		matrix invcov = syminv(cov)			// inverse Sigma matrix defined in Appendix A
		** Calculate Weights
		matrix unity  = J(1, rowsof(invcov), 1)
		matrix weights = unity * invcov		// simple column sum to get w_jk from Appendix A

		** Calculate Weighted Sums
		svmat weights, names(temp_weight_)
		forvalues count = 1/`nvars' {
			gen temp2_`count'  = temp1_`count'_z * temp_weight_`count'[1]
			gen temp3_weight_`count' = temp_weight_`count'[1] if !missing(temp1_`count'_z)
		}

		** Calculate Temporary Index
		egen 	temp_index 			= rowtotal(temp2_* )				// weighted sum from Appendix A
		egen 	temp_index_weight 	= rowtotal(temp3_weight_*)			// W_ij from Appendix A
		replace temp_index			= temp_index / temp_index_weight	// s_ij from Appendix A

		** Normalise Index
		qui su 	temp_index if ( `index_treat' == `index_ctrl' )
		replace temp_index = (temp_index - r(mean)) / r(sd)

		** Create or Replace Index
		cap gen 	`index_name' = temp_index	if `index_sample'
		cap replace `index_name' = temp_index	if `index_sample'

		** Drop Temporary Variables
		cap drop temp2_*
		cap drop temp3_*
		cap drop temp_cov_*
		cap drop temp_weight_*
		cap drop temp_index*
	}
	** Drop Temporary Variables
	cap drop temp1_*

	** Warning Message if Missing Variables
	if !missing("`index_var_missing'") di "error: variables `index_var_missing'missing"
}
else di "error: sample missing"
