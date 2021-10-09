********************************************************************************
*** PAIRMATCH COMMAND ***
********************************************************************************
/*
pairmatch.ado:
	performs nearest pairwise match merges of var_near within exact matches of
	optional varlist. Uses a similar syntax as the nearmrg code written by
	Booth, M Blasnik and K Smith.

Version 1.0 (2020-02-14) by Paul Brimble

This code achieves the following algorithm. For each group (identified with
varlist), execute the following:
	1) Calculate the differences between var_near for all pairwise
	   combinations of observations from the main and using data and denote
	   these with {Ax, By} where Ax denotes observation x from the main
	   data and By denotes observation y from the using data.

	2) Identify the minimum number of distinct observations in the main and
	   using data and denote this with N. If there are equal observations
	   in both the main and using data, there will be no unmatched data.

	3) For each integer n in {1,...,N}":
		3.1) Identify the pairwise combination with the n^th smallest
		     difference and label the elements of this combination {An, Bn}.
		3.2) Omit all other pairwise combinations that contain either An or
		     Bn. Specifically omit all pairwise combinations {An, 'Bn} and
			 {'An, Bn} where ' denotes "not".

	4) For each group, you are left with a set of unique nearest pairwise
	   combinations {A1,B1}, ... ,{An,Bn}.

A note on tiebreaks: in cases where there is a tie for minimum differences, this
is only problematic if there is overlap with at least one of the elements. For
example, suppose {Ax,Bx} and {Ay,By} have the same difference, then the order
of the algorithm does not make a difference. However, if we have the same
difference for {Ax,Bx} and {Ax,By}, then the arbitrary ordering from the
'stable' option will make a difference. This is because if {Ax, Bx} is chosen,
the eventual pairing of By will be affected. The intuition is further
illustrated if we have N = 1 with one main observation and two using
observations. If the main observation is equally close to the using observation,
then an arbitrary decision is made.
*/

program define pairmatch
	syntax [varlist(default=none)] using ,  ///
		var_near(varname numeric) ///
		[LIMit(str asis) SIGNed(str asis) ORDERing(str asis) GENMATCH(str asis)   * ]
	****************************************************************************
	*** A) PRElIMINARIES ***
	****************************************************************************
	** Check Options Are Consistent
	if "`ordering'" != "" & !inlist("`ordering'","strict","weak") {
		di as err "Ordering option mispecified; it must be strict or weak"
		error 461
	}
	if "`signed'" != "" & !inlist("`signed'","positive","negative") {
		di as err "Signed option mispecified; it must be positive or negative"
		error 461
	}

	** Default Optoins
	if "`ordering'" == "" {
		local ordering weak
	}

	if `"`genmatch'"'!=`""' confirm new var `genmatch'

	** Identify All Variables
	local vars_full `"`varlist' `var_near'"'

	** Preserve Current Main File
	sort 	 `vars_full', stable
	cap isid `vars_full' // check uniqueness
	if _rc {
		di as err "Variables: `vars_full' not unique in main dataset"
		error 459
	}
	tempfile  dta_main_full
	qui save `dta_main_full'

	** Identify Relevant Variables in Using Data
	use `vars_full' `using', clear

	** Check Near Variable
	if mi(`var_near') {
		local lbl_error = "`var_near' contains non-numeric chars, cannot "	///
			+ "destring (if `var_near' is a date, convert using time-date " ///
			+ "functions)"
		di as err "`lbl_error'"
		exit 198
	}
	sort 	 `vars_full', stable
	cap isid `vars_full' // check uniqueness
	if _rc {
		di as err "Variables: `vars_full' not unique in using dataset"
		error 459
	}

	** Create Using Copy of var_near (Labelled with Suffix "_2")
	tempvar   var_near_2
	clonevar `var_near_2' = `var_near'
	drop 	 `var_near'

	** Temporary Save Relevant Using Data
	tempfile  dta_using
	qui save `dta_using'

	** Identify Relevant Variables in Main Data
	use 	 `vars_full' using "`dta_main_full'", clear
	tempvar   var_near_1
	clonevar `var_near_1' = `var_near'
	drop 	 `var_near'

	****************************************************************************
	*** B) NEAREST MATCHING ***
	****************************************************************************
	** Joinby Data to Form All Pairwise Comparisons
	joinby 	 `varlist' 	using "`dta_using'"

	** Generate Absolute Differences
	tempvar dif
	gen double `dif' = abs(`var_near_1' - `var_near_2')

	** Signed and Limit Options
	if "`signed'" == "" {
		** Omit Pairs with Differences Exceeding Limit (If Specified)
		if 		`"`limit'"' != "" & "`ordering'" == "weak" {
			keep if `dif' <= `limit'
		}
		else if `"`limit'"' != "" & "`ordering'" == "strict" {
			keep if `dif' < `limit'
		}
	}
	else {
		** Create Limit of 0 If No Limit
		if "`limit'" == "" {
			local limit = 0	// assign limit of 0
		}
		** Generate Signed Differences
		tempvar dif_sign
		gen double `dif_sign' = `var_near_1' - `var_near_2'
		** Omit Pairs
		if 		"`signed'" == "positive" & "`ordering'" == "weak" {
			keep if `dif_sign' >= `limit'
		}
		else if "`signed'" == "positive" & "`ordering'" == "strict" {
			keep if `dif_sign' > `limit'
		}
		else if "`signed'" == "negative" & "`ordering'" == "weak" {
			keep if `dif_sign' <= `limit'
		}
		else if "`signed'" == "negative" & "`ordering'" == "strict" {
			keep if `dif_sign' < `limit'
		}
	}

	** Sort in Ascending Order of Differences
	sort `varlist' `dif', stable

	** Omit Non-nearest Pairs
	qui duplicates drop `varlist' `var_near_1', force
	qui duplicates drop `varlist' `var_near_2', force

	****************************************************************************
	*** C) MERGE DATA ***
	****************************************************************************
	** Rename Variables to Merge Original Using Data
	rename `var_near_2' `var_near'
	keep   `vars_full'  `var_near_1'
	/*
	note: `var_near' is the merging variable for the using data. `var_near_1'
	will be the eventual merging variable for the main data.
	*/

	** Merge Original Using Data
	qui merge 1:1 `vars_full'  `using', nogen

	** Keep a Copy of Nearest Variable from Using Data
	if "`genmatch'"!="" {
		clonevar `genmatch'=`var_near'
		if `"`vars_full'"' !="" local text ", matched on `varlist'"
		label var `genmatch' `"nearest match to `var_near' `text'"'
	}

	** Create Main Data Nearest Variable for Matching
	replace `var_near' = `var_near_1' if !mi(`var_near_1')

	** Save Modified Using Data
	tempfile  dta_match_using
	qui save `dta_match_using'

	** Merge Modified Using into Main Data
	use "`dta_main_full'", clear
	merge 1:1 `vars_full' using "`dta_match_using'", `options'

end
