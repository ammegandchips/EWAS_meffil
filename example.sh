#!/bin/bash
#
#
#PBS -l nodes=1:ppn=1,walltime=12:00:00

WORK_DIR="/panfs/panasas01/sscm/gs8094/EWAS/example_project"
module add languages/R-3.0.2

R CMD BATCH --no-save --no-restore '--args Trait CellData CellAdj Phenofile Method Removal BorM TP PACE Covariates Crude_or_Adj WD' /panfs/panasas01/sscm/gs8094/Common_files/meffil_EWAS_script.r /panfs/panasas01/sscm/gs8094/EWAS/example_project/example.out