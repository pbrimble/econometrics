********************************************************************************
*** SUMMARY INDEX ***
********************************************************************************
program drop sumindex
program def sumindex, byable(recall)
	syntax varlist(numeric) [if] [in], Generate(name) ///
		[ Base(string asis) Replace Pairwise Normalise ]

	****************************************************************************
	*** A) PRElIMINARIES ***
	****************************************************************************
	** Marksample
    marksample touse, novarlist strok

	** Check Base Condition
	if "`base'" != "" {
		capture count if `base'
		if _rc {
			display as error ///
				"base condition cannot be evaluated"
			exit 198
		}
	}

	** Create Base Marker
	tempvar tousebase
	if "`base'" != "" {
		gen `tousebase' = (`base' & `touse')
	}
	else {
		gen `tousebase' = `touse'
	}

	** Checks New Variable Name is Appropriate
	if _byindex() == 1 {
	    capture confirm new variable `generate'
		if !_rc {
			qui gen `generate' = .
		}
		else if "`replace'" == "" {
			local lbl_error = "generate() should give new variable name "	///
							+ "or replace option should be selected"
			display as error	"`lbl_error'"
			exit 198
		}
	}


	****************************************************************************
	*** B) INDEX CREATION ***
	****************************************************************************

	** Index Components
	local N: word count `varlist'

	** Create Normalised Variables (Identified with "a")
	tokenize `varlist'
	local nobasevalues = 0
	forvalues n = 1/`N' {
		tempvar a`n'
		qui summarize ``n'' if `tousebase'
		if r(N) > 0 {
			qui gen `a`n'' = (``n'' - r(mean)) / r(sd) if `touse'
		}
		else {
			local nobasevalues = 1
			continue, break
		}
	}

	** Complete Program if No Base Values
	if `nobasevalues' == 1 {
		local lbl_error = "warning: no base group observations; "	///
						+ "programme will return missing index"
		display "`lbl_error'"
	}
	else {
		** Temporary Matrices
		tempname cov invcov unity weights

		** Make Covariance Matrix
		if "`pairwise'" == "" {
			local A = ""
			forvalues n = 1/`N' {
				local A = "`A' `a`n''"
			}
			qui correl `A' if `tousebase'
			matrix `cov' = r(C)
		}
		else {
			matrix `cov' = I(`N')
			forvalues i = 1/`N' {
				forvalues j = 1/`N' {
					if `i' >= `j' {
						qui correl `a`i'' `a`j'' if `tousebase', covariance
						matrix `cov'[`i',`j'] = r(cov_12)
						matrix `cov'[`j',`i'] = r(cov_12)
					}
				}
			}
		}

		** Calculate Weights
		matrix `invcov' 	= syminv(`cov')
		matrix `unity' 		= J(1, rowsof(`invcov'), 1)
		matrix `weights' 	= `unity' * `invcov'

		** Calculate Weighted Sums (Identified with "b" and "c")
		forvalues n = 1/`N' {
			tempvar b`n' c`n'
			qui gen `b`n'' = `a`n'' * `weights'[1,`n'] if `touse'
			qui gen `c`n'' = `weights'[1,`n']	if !mi(`a`n'') & `touse'
		}

		** Calculate Index
		local B = ""
		forvalues n = 1/`N' {
			local B = "`B' `b`n''"
		}
		local C = ""
		forvalues n = 1/`N' {
			local C = "`C' `c`n''"
		}
		tempvar d e f
		qui egen 	`d' = rowtotal(`B')	if `touse'
		qui egen 	`e' = rowtotal(`C') if `touse'
		qui gen 	`f' = `d' / `e'		if `touse'

		** Normalise (Optional)
		if "`normalise'" != "" {
			qui summarize `f' if `tousebase'
			qui replace `f' = (`f' - r(mean)) / r(sd) if `touse'
		}

		** Finalise Index
		replace `generate' = `f' if `touse'
	}

end
