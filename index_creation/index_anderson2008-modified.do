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
		gen temp1_`i' = `x'		if `index_sample'		// sample-restricted variables
	}
	local nvars = `i'

	** Standardize Variables With Respect To Control Group
	local index_var_missing ""
	forvalues count = 1/`nvars' {
		qui su temp1_`count'  		if `index_treat' == `index_ctrl' // normalise by control group
		if r(N) > 0	gen temp2_`count'  = (temp1_`count' - r(mean)) / r(sd)
		// normalised sample-restricted variables
		else local index_var_missing `index_var_missing' `count' " " // check if entire variable is missing for control group
	}

	** Check if All Variables Exist
	if missing("`index_var_missing'") {
		** Create Anderson's Covariance Matrix of Variables
		matrix cov = I(`nvars')				// Sigma matrix defined in Appendix A
		forvalues i = 1/`nvars' {
			forvalues j = 1/`nvars' {
				if `i' >= `j' {
					** Sample-Adjusted Pairwise Covariance Matrix)
					correl temp2_`i' temp2_`j' if `index_treat' == `index_ctrl', covariance
					matrix cov[`i',`j'] = r(cov_12)
					matrix cov[`j',`i'] = r(cov_12)
				}
			}
		}
		** Standard Covariance Matrix
		*correl temp2_* 	if  `index_treat' == `index_ctrl' , covariance
		*matrix cov = r(C)

		matrix invcov = syminv(cov)			// inverse Sigma matrix defined in Appendix A
		** Calculate Weights
		matrix unity  = J(1, rowsof(invcov), 1)
		matrix weights = unity * invcov		// simple column sum to get w_jk from Appendix A

		** Calculate Weighted Sums
		svmat weights, names(temp3_) 		// variable weights
		forvalues count = 1/`nvars' {
			gen temp4_`count' = temp2_`count' * temp3_`count'[1]
			// product of weights and normalised sample-restricted variables
			gen temp5_`count' = temp3_`count'[1] if !missing(temp2_`count')
			// observation specific weighted sums
		}

		** Calculate Temporary Index
		egen 	temp_index 			= rowtotal(temp4_* )				// weighted sum from Appendix A
		egen 	temp_index_weight 	= rowtotal(temp5_*)					// W_ij from Appendix A
		replace temp_index			= temp_index / temp_index_weight	// s_ij from Appendix A

		** Normalise Index
		qui su 	temp_index if ( `index_treat' == `index_ctrl' )
		replace temp_index = (temp_index - r(mean)) / r(sd)

		** Create or Replace Index
		cap gen 	`index_name' = temp_index	if `index_sample'
		cap replace `index_name' = temp_index	if `index_sample'

	}
	** Drop Temporary Variables
	cap drop temp1_*		// sample-restricted variables
	cap drop temp2_*		// normalised sample-restricted variables
	cap drop temp3_*		// variable weights
	cap drop temp4_*		// product of weights and normalised sample-restricted variables
	cap drop temp5_*		// observation specific weights
	cap drop temp_index*	// temporary index

	** Warning Message if Missing Variables
	if !missing("`index_var_missing'") di "error: variables `index_var_missing'missing"
}
else di "error: sample missing"
