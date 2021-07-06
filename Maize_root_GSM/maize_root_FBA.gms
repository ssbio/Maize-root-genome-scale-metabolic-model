*Flux Balance Analysis Program

OPTIONS

reslim=100000000000000000
decimals=8
optca=1E-9
optcr=1E-9

;

$INLINECOM /*  */
$onlisting
$onempty

SETS

	i		set of metabolites
$include "mets_1_filled.txt"

	j		set of reactions
$include "rxns_1_filled.txt"

;

PARAMETERS

	S(i,j)	stoichiometric matrix
$include "Sij_1_filled.txt"

	rxn_type(j)		reaction type information
$include "rxntype_1_filled.txt"

;

EQUATIONS

	objective
	
	mass_balance(i)
;

VARIABLES

	Z		objective value
	v(j)		reaction rate

;

objective..		Z =e= v('Root_Biomass[GS1]');

*objective..		Z =e= v('R01281[R,c,GS1]');

mass_balance(i)..	sum(j, S(i,j) * v(j)) =e= 0;

MODEL FBA
/

	objective
	mass_balance

/
;

v.lo(j)$(rxn_type(j) eq 1) = 0;
v.up(j)$(rxn_type(j) eq 1) =  30000;

v.lo(j)$(rxn_type(j) eq 2) = -30000;
v.up(j)$(rxn_type(j) eq 2) = 0;

v.lo(j)$(rxn_type(j) eq 3) = -30000;
v.up(j)$(rxn_type(j) eq 3) =  30000;

v.lo(j)$(rxn_type(j) eq 4) = 0;
v.up(j)$(rxn_type(j) eq 4) = 30000;

v.lo(j)$(rxn_type(j) eq 5) = 0;
v.up(j)$(rxn_type(j) eq 5) = 0;

**********Nutrients********************
*Phosphate
v.lo('Exe3[R,GS1]') = -600;
*Ammonia
v.lo('Exe5[R,GS1]') = -30000;
*Nitrate
v.lo('Exe10[R,GS1]') = -30000;
*Sulfate
v.lo('Exe6[R,GS1]') = -30000;
*Manganese
v.lo('Exe_C00034[R,GS1]') = -30000;
*Zinc
v.lo('Exe_C00038[R,GS1]') = -30000;
*Copper
v.lo('Exe_C00070[R,GS1]') = -30000;
*Calcium
v.lo('Exe50[R,GS1]') = -30000;
*Pottasium
v.lo('Exe9[R,GS1]') = -30000;
*Magnesium
v.lo('Exe19[R,GS1]') = -30000;
*Chloride
v.lo('Exe13[R,GS1]') = -30000;
*Sodium
v.lo('Exe18[R,GS1]') = -30000;
*Ferrous
v.lo('Exe21[R,GS1]') = -30000;
****************************************

*********Phosphorus Phloem Transport***************
*Uridine diphosphate
v.lo('PhloemLoad_C00015_I[P,GS1]') = -30000;
*Cytidine triphosphate
v.lo('PhloemLoad_C00063_I[P,GS1]') = -30000;
*Uridine triphosphate
v.lo('PhloemLoad_C00075_I[P,GS1]') = -30000;
*Guanosine triphosphate
v.lo('PhloemLoad_C00044_I[P,GS1]') = -30000;
**************************************************

*********Nitrogen Phloem Transport*****************
*Threonine
*v.lo('PhloemLoad_C00188[P,GS1]') = -30000;
*Serine
*v.lo('PhloemLoad_C00065[P,GS1]') = -30000;
*Glycine
*v.lo('PhloemLoad_C00037[P,GS1]') = -30000;
*Alanine
*v.lo('PhloemLoad_C00041[P,GS1]') = -30000;
*Cysteine
*v.lo('PhloemLoad_C00097[P,GS1]') = -30000;
*Methionine
*v.lo('PhloemLoad_C00073[P,GS1]') = -30000;
*Isoleucine
*v.lo('PhloemLoad_C00407[P,GS1]') = -30000;
*Tyrosine
*v.lo('PhloemLoad_C00082[P,GS1]') = -30000;
*Phenylalanine
*v.lo('PhloemLoad_C00079[P,GS1]') = -30000;
*Leucine
*v.lo('PhloemLoad_C00123[P,GS1]') = -30000;
*Histidine
*v.lo('PhloemLoad_C00135[P,GS1]') = -30000;
*Arginine
*v.lo('PhloemLoad_C00062[P,GS1]') = -30000;
**************************************************

*********N-P Phloem Transport*****************
*Uridine diphosphate
v.lo('PhloemLoad_C00105_I[P,GS1]') =-30000;
**************************************************


SOLVE FBA USING LP MAXIMIZING Z;
*SOLVE FBA USING LP MINIMIZING Z;

FILE RESULTS /FBA_result.txt/;

PUT RESULTS;

PUT "reaction      FLUX"/;

LOOP(j,
	
	PUT j.tl:0:100,"    ", v.l(j):20:5/;
		
);

PUTCLOSE;