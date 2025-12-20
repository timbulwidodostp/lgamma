*! version 1.0.3  17jun1999  (STB-53: sg126)
program define gammalog, eclass
	version 6
	if replay() {
		if "`e(cmd)'"!="gammalog" {
			error 301
		}
		Display `0'
		error `e(rc)'
		exit
	}

/* Parse. */

	syntax varlist(numeric) [fw pw iw] [if] [in] [, IRr EForm /*
	*/ FROM(string) Level(integer $S_level) Offset(varname numeric) /*
	*/ Exposure(varname numeric) noCONstant Robust CLuster(varname) /*
	*/ SCore(string) noLOg * ]

	mlopts mlopts, `options'

/* Check syntax. */

	if `level' < 10 | `level' > 99 {
		di in red "level() must be between 10 and 99"
		exit 198
	}
	if `"`score'"'!="" {
		confirm new variable `score'
		local nword : word count `score'
		if `nword' > 2 {
			di in red "score() must contain the name of " /*
			*/ "one or two new variables"
			exit 198
		}
		local scname1 : word 1 of `score'
		local scname2 : word 2 of `score'
		tempvar scvar1 scvar2
		local scopt "score(`scvar1' `scvar2')"
	}
	if "`offset'"!="" & "`exposur'"!="" {
		di in red "only one of offset() or exposure() can be specified"
		exit 198
	}
	if "`constan'"!="" {
		local nvar : word count `varlist'"
		if `nvar' == 1 {
			di in red "independent variables required with " /*
			*/ "noconstant option"
			exit 102
		}
	}

/* Mark sample except for offset/exposure. */

	marksample touse

	if `"`cluster'"'!="" {
		markout `touse' `cluster', strok
		local clopt cluster(`cluster')
	}

/* Process offset/exposure. */

	if "`exposur'"!="" {
		capture assert `exposur' > 0 if `touse'
		if _rc {
			di in red "exposure() must be greater than zero"
			exit 459
		}
		tempvar offset
		qui gen double `offset' = ln(`exposur')
		local offvar "ln(`exposur')"
	}

	if "`offset'"!="" {
		markout `touse' `offset'
		local offopt "offset(`offset')"
		if "`offvar'"=="" {
			local offvar "`offset'"
		}
	}

/* Count obs and check for negative values of `y'. */

	gettoken y xvars : varlist

	if "`weight'"!="" {
		if "`weight'"=="fweight" {
			local wt `"[fw`exp']"'
		}
		else	local wt `"[aw`exp']"'
	}

	summarize `y' `wt' if `touse', meanonly

	if r(N) == 0 { error 2000 }
	if r(N) == 1 { error 2001 }

	if r(min) <= 0 {
		di in red "`y' must be greater than zero"
		exit 415
	}

	tempname mean
	scalar `mean' = r(mean)

/* Remove collinearity. */

	_rmcoll `xvars' [`weight'`exp'] if `touse', `constan'
	local xvars `r(varlist)'

/* Estimate constant-only model. */

	if "`constan'"=="" & `"`from'"'=="" {
		tempname b00

	/* Get starting values for full model. */

		if "`xvars'"!="" {
			tempname b0
			tempvar z
			if "`offset'"=="" {
				qui gen double `z' = ln(`y') if `touse'
			}
			else	qui gen double `z' = ln(`y')-`offset' if `touse'

			qui reg `z' `xvars' `wt' if `touse'
			drop `z'
			mat `b0' = e(b)
			local dim = colsof(`b0')
			mat `b0' = `b0'[1,1..`dim'-1]
			mat coleq `b0' = `y'

			qui mat score double `z' = `b0' if `touse'
			if "`offset'"!="" {
				qui replace `z' = `z' + `offset'
			}

			StartVal `touse' `y' `mean' `"`wt'"' `z' `b00'
			drop `z'
			mat `b0' = (`b0',`b00')
			local initopt "init(`b0')"
		}

	/* Fit constant-only model. */

		if "`log'"=="" {
			di in gr _n "Fitting constant-only model:"
		}

	/* Get starting values for constant-only model. */

		StartVal `touse' `y' `mean' `"`wt'"' `"`offset'"' `b00'

		if "`weight'"!="" {
			local iwt `"[iw`exp']"'
		}

		ml model d2 gamln_lf (`y': `y'=, `offopt') /ln_phi /*
		*/ `iwt' if `touse', collinear missing max nopreserve /*
		*/ wald(0) init(`b00') search(off) `mlopts' `log'

		local continu "continue"

		if "`log'"=="" {
			di in gr _n "Fitting full model:"
		}
	}

/* Estimate full model. */

	if `"`from'"'!="" { local initopt `"init(`from')"' }

	ml model d2 gamln_lf (`y':`y'=`xvars', `constan' `offopt') /ln_phi /*
	*/ [`weight'`exp'] if `touse', collinear missing max nopreserve /*
	*/ `initopt' search(off) `mlopts' `log' `scopt' `robust' `clopt' /*
	*/ `continu' title("Gamma distribution-log link model")

	est local cmd

/* Make score variables real. */

        if "`score'"!="" {
		label var `scvar1' "Score index for x*b from gammalog"
		rename `scvar1' `scname1'
		if "`scname2'"!="" {
			label var `scvar2' /*
			*/ "Score index for /ln_phi from gammalog"
			rename `scvar2' `scname2'
		}
	}

/* Fill in e(). */

	est scalar phi = exp(_b[/ln_phi])
	est local offset  "`offvar'"
	est local offset1 /* erase; set by -ml- */
	est local predict "poisso_p"
	est local cmd "gammalog"

/* Display final results. */

	Display, `irr' `eform' level(`level')

	error `e(rc)'
end

program define Display
	syntax [, IRr EForm Level(int $S_level)]
	if `level' < 10 | `level' > 99 {
		di in red "level() must be between 10 and 99"
		exit 198
	}
	if "`irr'"!="" | "`eform'"!="" {
		local eopt "eform(IRR)"
	}

	ml display, level(`level') `eopt' plus first
	_diparm ln_phi, level(`level')
	di in gr _dup(9) "-" "+" _dup(68) "-"
	_diparm ln_phi, level(`level') exp label("phi")
	di in gr _dup(78) "-"
end

program define StartVal
	args touse y mean wt offset b0

	tempname c f0 lnphi olnphi
	tempvar z

	if "`offset'"=="" {
		scalar `c' = ln(`mean')
		local offset 0
	}
	else {
		qui gen double `z' = `y'*exp(-`offset') if `touse'
		summarize `z' `wt' if `touse', meanonly
		scalar `c' = ln(r(mean))
		drop `z'
	}

	qui gen double `z' = `c' + `offset' - ln(`y') if `touse'
	summarize `z' `wt' if `touse', meanonly
	scalar `f0' = r(mean)

	scalar `lnphi' = (ln(`f0') + 0.3559669)/1.047473
	local converg 0
	local i 1
	while `i'<=10 & !`converg' {
		scalar `olnphi' = `lnphi'
		scalar `lnphi' = (digamma(exp(-`lnphi'))+`lnphi'+`f0') /*
		*/ /(exp(-`lnphi')*trigamma(exp(-`lnphi'))-1) + `lnphi'
		local converg = (reldif(`lnphi',`olnphi')<1e-8)
		local i = `i' + 1
	}

	mat `b0' = (0, 0)
	mat colnames `b0' = `y':_cons ln_phi:_cons
	mat `b0'[1,1] = cond(`c'!=.,`c',0)
	mat `b0'[1,2] = cond(`lnphi'!=.,`lnphi',0)
end

