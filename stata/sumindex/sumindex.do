********************************************************************************
*** SUMMARY INDEX ***
********************************************************************************

cap program drop sumindex
program def sumindex,
	syntax namelist [if] [in], Base(str) [Replace Pairwise Normalise]

	** Identify Index Name and Index Components
	local N = 0
	foreach name of local namelist {
		** Counter
		local N = `N' + 1
		**  Index Variable
		if `N' == 1 {
			local index 	= "`name'"
			local `index' 	= ""
		}
		else {
			local `index' = "``index'' `name'"
		}
	}
	** Size of Subcomponents
	local N = `N' - 1

	** Marksample
    marksample touse, strok

	** Create Normalised Variables
	tokenize ``index''
	forvalues n = 1/`N' {
		tempvar a`n'
		su ``n'' if `base' & `touse'
		gen `a`n'' = (``n'' - r(mean)) / r(sd) if `touse'
	}

	** Make Covariance Matrix
	tempname cov invcov unity weights
	matrix `cov' = I(`N')
	if "`pairwise'" == "" {
		local A = ""
		forvalues n = 1/`N' {
			local A = "`A' `a`n''"
		}
		correl `A'
		matrix `cov' = r(C)
	}
	else {
		forvalues i = 1/`N' {
			forvalues j = 1/`N' {
				if `i' >= `j' {
					correl `a`i'' `a`j'' if `base', covariance
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

	** Calculate Weighted Sums
	forvalues n = 1/`N' {
		tempvar b`n' c`n'
		gen `b`n'' = `a`n'' * `weights'[1,`n']
		gen `c`n'' = `weights'[1,`n']	if !mi(`a`n'')
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
	egen 	`d' = rowtotal(`B')
	egen 	`e' = rowtotal(`C')
	gen 	`f' = `e' / `d'

	** Normalise (Optional)
	if "`normalise'" != "" {
		su `f' if `base'
		replace `f' = (`f' - r(mean)) / r(sd)
	}

	** Finalise Index
	if "`replace'" == "" {
		gen `index' = `f'
	}
	else {
		replace `index' = `f' if `touse'
	}

end
