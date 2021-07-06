*************************************************************
*E-Flux algorithm for integrating transciptomics data
*************************************************************
*****************Niaz Bahar Chowdhury************************
*************************************************************
$INLINECOM /*  */

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "mets_1_filled.txt"	

	j					set of reactions

$include "rxns_1_filled.txt"

;

alias(j,j1);
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "Sij_1_filled.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"

	c(j)					choose set for FVA
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance(i)				steady state mass balance
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
*	biomass_enforcement			enforcing biomass
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	Z					objective value
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= sum(j, c(j) * v(j));

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

*biomass_enforcement..		v('Root_Biomass[GS1]') =e= 1936.542316;

Model maize_root /all/;
******************************************************************

****************Flux Variability Analysis*************************

maize_root.optfile = 1;

maize_root.holdfixed = 1;

FILE RESULTS /WT_FVA.txt/;

PUT RESULTS;

PUT "reaction    WT_min_rate     WT_max rate"/;

LOOP(j1,

	c(j) = 0;
	c(j1) = 1;
	
	PUT j1.tl:0:100;

	solve maize_root using lp MINIMIZING Z;
	
	PUT Z.l:20:5;
	
	solve maize_root using lp MAXIMIZING Z;
	
	PUT Z.l:20:5/;
	
);

PUTCLOSE;
