# Multi-state model solution to the leukemia data

Treatment: transplant modality (MSDT, Haplo-SCT) 

Terminal (primary) event: death

Intermediate event: relapse

Baseline confounders: age, sex, diagnosis, disease status

Time-varying confounder: GVHD

## leukemia_mediation.R
Under the separable effects framework, to estimate the counterfactual cumulative indicences (CIFs) and treatment effects.
This is included in "leukemia_interventional.R" as a special case.

## leukemia_interventional.R
Under the separable effects, interventional effects and conditional interventional effects framework, to estimate the counterfactual cumulative indicences (CIFs) and treatment effects (bootstrap use "leukemia_bootstrap.R").
