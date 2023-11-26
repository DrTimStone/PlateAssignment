# PlateAssignment

This is a script to assign 192 cancer and non-cancer samples onto two Illumina Infinium MethylationEPIC plates (epigenetic array plates) in a way that optimises the balance between several clinical and lifestyle covariates.

This script is used after executing the CaseControlMatching.R script, that finds 192 matched case-controls from each of the groups in our cohort: 3 cancer groups (adenocarcinoma, intramucosal and high-grade dysplasia, and 3 controls ("normals", healthy volunteers, non-dysplastic Barrett's esophagus).

The motivation for this analysis was to disassociate the huge batch effect of the epigenetic plates from the clinical covariates. We would be using stringent batch correction for the technical variables, which are the two plates, plus the individual rows of the plate. This removal is achieved by means of the residuals of a fitted linear model, we would preserve the clinical covariates so that they could be investigated controlled for and investigated seperately. This is therefore a script that seperates out the main (uninteresting) technical covariate from the more interesting clinical and lifestyle covariates. 

The covariates are: Case-Control (2-group), then the sub-groups, age, sex, bmi, alcohol consumption, cigarette consumption, PPI consumption and accumulated instance of heartburn. 

The most original part of this script was to create a measure of homogeneity based on weighted p-values, used in the reverse way a p-value is normally used. P-values are normally used to indicate differences, but it follows that the more similar and homogeneous two samples are, the higher the p-value. We can therefore perform a seires of statistical tests and weight them according to the covariates that we consider more important, add them all up and search random permutations for a maximum. 

It is usually assumed that balance between the clinical groups has the highest weighting, followed by sex differences, and maybe heartburn and PPI differences come last, but the parameters of these can be set by the user.

I assume that a near-even distribution of cases and controls on each plate is the optimal solution and have set up the script not to deviate too far from this situation, however, an increase of a very small number number of cases on one plate (and a corresponding decrease in the other) COULD produce learge dividends in balancing out the other main covariates, but the script has a limit on the extent of this "jittering".

Later versions of this script incorporate the individual rows of the array plate for balance. 
