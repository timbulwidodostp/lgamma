*! version 1.0.1  20apr1999  (STB-53: sg126)
/* log-likelihood program to accompany gammalog.ado */
program define gamln_lf
	version 6
	args todo b lnf g H s1 s2

/* Calculate the log-likelihood. */

	tempvar xb
	tempname lnphi iphi
	mleval `xb' = `b', eq(1)
	mleval `lnphi' = `b', eq(2) scalar
	scalar `iphi' = exp(-`lnphi')
	local y $ML_y1

	mlsum `lnf' = (`iphi'-1)*ln(`y') - lngamma(`iphi') /*
	*/ - `iphi'*(`y'*exp(-`xb') + `xb' + `lnphi')

	if `todo' == 0 | `lnf'==. { exit }

/* Calculate the scores and gradient. */

	tempname g1 g2
	qui replace `s1' = -`iphi'*(1 - `y'*exp(-`xb')) if $ML_samp
	qui replace `s2' = `iphi'*(-ln(`y') + digamma(`iphi') /*
	*/ + `y'*exp(-`xb') + `xb' + `lnphi' - 1) if $ML_samp

	mlvecsum `lnf' `g1' = `s1', eq(1)
	mlvecsum `lnf' `g2' = `s2', eq(2)
	matrix `g' = (`g1',`g2')

	if `todo' == 1 | `lnf'==. { exit }

/* Calculate the negative hessian. */

	tempname d11 d12 d22
	mlmatsum `lnf' `d11' = `iphi'*`y'*exp(-`xb'), eq(1)
	mlmatsum `lnf' `d12' = `s1', eq(1,2)
	mlmatsum `lnf' `d22' = `s2' - `iphi'*(1 - `iphi'*trigamma(`iphi')), /*
	*/ eq(2)

	matrix `H' = (`d11',`d12' \ `d12'',`d22')
end
