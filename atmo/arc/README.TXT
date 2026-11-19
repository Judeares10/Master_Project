ARC RETRIEVAL

The ATMO Retrieval Code (ARC)

Purpose is to find a best-fitting solution to eclipse or transit data (CHI^2 minimizer), measure the evidence of the model (nested sampling), and measure the posterior (MCMC). ARC is written and run in IDL.

IDL libraries are included so no user IDL setup should be required.

Outputs include various diagnostic plots.

ARC Programs are contained in a directory in pacode/atmo/arc/ and executed from a directory in /atmo/run/ as with a normal atmo run e.g. pacode/atmo/run/ARC_Wasp107/. The outputs are saved in a separate directory, suggested on /data/ as many of the outputs are very large (i.e. MCMC chains and spectra at each chain).

To run:

Follow normal ATMO install procedures, compile ATMO into a new directory within /atmo/run/ to execute the retrieval (e.g. atmo/run/ARC_Wasp107/). 1) Copy a setup file pacode/atmo/arc/ARC_retrieval_input_tr_eq.pro into the run directory (e.g pacode/atmo/run/ARC_Wasp107/)
Copy data to be fit into pacode/atmo/run/ARC_Wasp107/
EDIT arc_retrieval_input.pro to setup the retrieval (See QUICKGUIDE) and also enter basic planet and stellar data. An appropriate Stellar file will be needed.

Compile and Run in IDL
	IDL> .r arc_retrieval_input.pro
	IDL> arc_retrieval_input

An examples is included (/atmo/arc/arc_retrieval_input_Wasp121example.pro) which can be copied into a run directory /atmo/run/ARC_Wasp121_Trans_Fe/ 
along with the data /atmo/arc/w121_trspec.g430l_n30.g750l_n30.ldatlas_free.G430Lshifted.txt.  The example is a transmission spectra of Wasp-121b.
With atmo.x in the same run directory, there should be at minimum the three mentioned files in /atmo/run/ARC_Wasp121_Trans_Fe/ needed to run a retrieval.
Compile and Run Example in IDL
	IDL> .r arc_retrieval_input_Wasp121example.pro
	IDL> arc_retrieval_input_Wasp121example


Quick Guides - See Attachment, link at top of page

Suggested Fitting Sequence:

1) Start with CHI^2 minimizer (turn off nested sampling & mcmc)
2) Find CHI^2 minimum for a desired ATMO fit.
3) Adjust the number of parameters and complexity of the model as desired, making a new arc_retrieval_input_EXONAMEHERE.pro for each model. Saving each model fit in it's own unique directory.
3a) Approximate evidence with BIC using MPFIT and/or MCMC
3b) Find evidence of each model tried with Nested sampling. Make course nested sampling fit. Can then compare likelihood of each model.
4) Measure the posterior distribution with MCMC for the best model(s).
5) Publish great results.

Known Issue:
What to do when MCMC crashes when looking for the step scale?? This can easily happen after an MPFIT minimization finnishes but before full MCMC is started. MCMC first calculates parameter change that results in change of CHI^2 by 1, which informs the step size when the MCMC is first started. This value only has to be close, as the DEMC does the rest after the MCMC chains have started. Crash happens when a parameter is unconstrained and it can't find Delta CHI^2=1.

Solution: restore file "dirs+MCMCscale_seed.sav" and manually set the scale for each free parameter, saving as an array called "scale" in the file "dirs+MCMC.scale.sav "

Example:

Restore what scale information has been found thus far

IDL> dirs='/name/of/output/directory/here/it's/in/the/first/line/of/the/arc_retrieval_input.pro/file/'
IDL> restore,dirs+'MCMCscale_seed.sav'
Array with 2 elements for each free parameter (and only the free parameters, fixed parameters are not included).
IDL> print,MCMCSCALE
0.12458936 0.027898464 3.2157411 0.0086544147 0.81642712 5.9506114 23.681758 0.57607244 10.159983 12.227367 0.86924807 1.0924662
0.12742088 0.028737486 NaN 0.0087576922 1.4341310 3.7813923 NaN 0.71345100 NaN NaN 0.99088580 23.384344
Top row is the scale when increasing the free parameter by that amount, bottom row is the negative change in the parameter. Each column corresponds to a free parameter. Average of the two + and - values is recomended for the size of the scale. (fixed parameters are not included). 
Default value is 0.001, thus if MCMCSCALE has this value the crash happened before a measurement could be made, it's likely better to set the value to a more sensible one in this case.
IDL> scale=(MCMCSCALE[*,0]+MCMCSCALE[*,1])/2. 
IDL> print,scale
0.12600512 0.028317975 NaN 0.0087060534 1.1252790 4.8660019 NaN 0.64476172 NaN NaN 0.93006693 12.238405
Replace NaN's with either a value found above in MCMCSCALE (higlighted in greed above), or your best guess if both values are NaN's which is likely the case.
IDL> scale=[0.12600512, 0.028317975, 3.2157411, 0.0087060534 , 1.1252790, 4.8660019, 23.681758, 0.64476172 , 10.159983 , 12.227367 , 0.93006693, 12.238405] 
IDL> save,filename=dirs+'MCMC.scale.sav',scale
Now restart ARC, it should find the scale, and skip straight to starting the MCMC.

Key References:

Line et al. 2013, ApJ, 775, 137 - retrieval techniques
Guillot, T. 2010, A&A, 520, A27 - T-P profile parameterization
Marquardt 2009 - Levenberg-Markquart least-squares fit 
Eastman et al. 2014 - MCMC. Differential Evolution Chain Monte Carlo
Feroz et al. (2009, MNRAS, 398, 1601; 2009, MNRAS, 384, 449) - Multinested sampling
Example ARC Publications:

Evans et al. 2017 - WASP-121b eclipse spectra with inversion.
Wakeford et al. 2017 - HAT-P-26b transmission spectra with H2O.

History

Version 1.0 - 22 June 2017 - David Sing