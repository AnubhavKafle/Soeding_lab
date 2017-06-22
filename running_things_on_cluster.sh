#!/bin/sh
#BSUB -J Gxxx
#BSUB -q mpi
#BSUB -W 48:00
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -n 16
#BSUB -R scratch
#BSUB -R span[hosts=1]
#BSUB -a openmp

STUDY=Gxxx
BASEDIR="/scratch/akaphle/New_Genotype/${STUDY}"
SNPTEST="${HOME}/anubhav_GWAS/scripts/snptest_v2.5.2"
OUTDIR="${HOME}/anubhav_GWAS/snptest_results/${STUDY}"
GENO="${BASEDIR}/${STUDY}_QC_Chr6_imputed_modified"
PHENO="${BASEDIR}/${STUDY}_age_at_onset.sample"

if [ ! -d ${OUTDIR} ]; then
    mkdir -p ${OUTDIR}
fi

OUTFILE="`basename ${GENO}`.out"
$SNPTEST -data $GENO $PHENO -o ${OUTDIR}/${OUTFILE} -frequentist 1 -pheno pheno -method score
