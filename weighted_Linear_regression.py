
# Weighted Linear regression for duplicated data

# In[ ]:


import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import collections
import os, csv
import glob

SAMPLE_FIELDS = ['sid', 'sex', 'pheno', 'weight']
class SampleInfo(collections.namedtuple('_SampleInfo', SAMPLE_FIELDS)):
    __slots__ = ()
    
SNPINFO_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()


# In[ ]:


def dosage_to_genotype(this_dosage):
    genotype = np.array(([this_dosage[i+1] + 2 * this_dosage[i + 2] for i in range(0, len(this_dosage),3)])).astype(float)
    return genotype

def norm_binom(gt, f):
    gt = (gt - ( 2* f)) / np.sqrt( 2* f * (1 - f))
    return gt

def read_dosage(filepath):
    snpinfo = list()
    norm_genotype = list()
    with open(filepath, 'r') as genfile:
        for mline in genfile:
            mlinesplit = mline.split()
            this_snp = SnpInfo(rsid = mlinesplit[1],
                               bp_location = int(mlinesplit[2]),
                               alt_allele = mlinesplit[3],
                               ref_allele = mlinesplit[4])
            this_snp_dosage = np.array(mlinesplit[5:], dtype=float)
            genoss = dosage_to_genotype(this_snp_dosage)
            snpinfo.append(this_snp)
            freq = (np.sum(genoss) / (len(genoss) * 2))
            normalized = norm_binom(genoss, freq)
            norm_genotype.append(normalized)
    return snpinfo, norm_genotype




def read_expected_phenotype(filepath):
    samples = list()
    with open(filepath, 'r') as samfile:
        next(samfile)
        next(samfile)
        for mline in samfile:
            mlinesplit = mline.split()
            this_sample = SampleInfo(sid = mlinesplit[0],
                                     sex = int(mlinesplit[5]), # careful, assuming no NA
                                     pheno = float(mlinesplit[7]),
                                     weight = 1)
            samples.append(this_sample)
    return(samples)

def read_weighted_phenotype(filepath):
    samples = list()
    with open(filepath, 'r') as samfile:
        next(samfile)
        next(samfile)
        for mline in samfile:
            mlinesplit = mline.split()
            this_sample = SampleInfo(sid = mlinesplit[1],
                                     sex = int(mlinesplit[5]), # careful, assuming no NA
                                     pheno = float(mlinesplit[6]),
                                     weight = float(mlinesplit[7])) # careful, assuming no NA 
            samples.append(this_sample)
    return(samples)

def get_pval(beta, var_beta):
    stat = np.square(np.array(beta))/var_beta
    pval = 1 - stats.chi2.cdf(stat, 1)
    return pval
    

def snptest(x, y):
    xtx = np.einsum('i, i', x, x)
    xty = np.einsum('i, i', x, y)
    beta = xty / xtx
    yres = y - x * beta
    var_beta = np.einsum('i, i', yres, yres) / (x.shape[0] * xtx)
    pval = get_pval(beta, var_beta)
    return beta, var_beta, pval
    
def weighted_snptest(gt, pheno, wt, indices):
    xnew = np.array([gt[indices[i]] for i in range(len(pheno))])
    xtwx = np.einsum('i, i, i', xnew, xnew, wt)
    xtwy = np.einsum('i, i, i', xnew, wt, pheno)
    beta = xtwy / xtwx
    yres = pheno - xnew * beta
    var_beta = np.einsum('i, i', yres, yres) / (xnew.shape[0] * xtwx)
    pval = get_pval(beta, var_beta)
    return beta, var_beta, pval


# In[ ]:


dosagefilename = '/home/anubhavk/Desktop/GWAS_Anubhav/expected_phenotype/G1_QC_Chr1.imputed_modified.top10lines'
phenofilename = '/home/anubhavk/Desktop/GWAS_Anubhav/expected_phenotype/G1_age_at_onset.sample_1'
weighted_phenofilename = '/home/anubhavk/Desktop/GWAS_Anubhav/distributed_phenotype/G1_age_duplicated_sample_diffid'
outputfile = '/home/anubhavk/Desktop/GWAS_Anubhav/distributed_phenotype/pppp.txt'

snpinfo, gt = read_dosage(dosagefilename)
pheno_expected = read_expected_phenotype(phenofilename)
pheno_weighted = read_weighted_phenotype(weighted_phenofilename)


# In[ ]:


pheno = np.array([x.pheno for x in pheno_weighted])
norm_pheno = (pheno - np.mean(pheno))/np.std(pheno)
weights = np.array([x.weight for x in pheno_weighted])
_sids = [x.sid for x in pheno_expected]
indices = np.array([_sids.index(x.sid) for x in pheno_weighted])

with open(outputfile, "w") as output:
    output.write("rsid pos allele_A allele_B info P_value beta se\n")
    for i, snp in enumerate(snpinfo):
        beta, var_beta, pval = weighted_snptest(gt[i], norm_pheno, weights, indices)
        se = np.sqrt(var_beta)
        output.write('{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(snp.rsid, snp.bp_location,snp.alt_allele, snp.ref_allele,'1',pval, beta, se))
