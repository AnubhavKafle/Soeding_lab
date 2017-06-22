
# coding: utf-8

# In[138]:

import numpy as np
import collections
import csv
import glob

SAMPLE_FIELDS = ['sid', 'sex', 'pheno']
class SampleInfo(collections.namedtuple('_SampleInfo', SAMPLE_FIELDS)):
    __slots__ = ()
    
SNPINFO_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

samplefilename = 'clinical_pheno/G1_QC.sample'
phenofilename = 'clinical_pheno/1_GerMIFS_I_formated_pheno.txt'
expctfilename = 'expected_age_at_onset_male.dat'
expctfilename_female = 'expected_age_at_onset_female.dat'
sampleoutfilename = 'G1_age_at_onset.sample'
gtfpath = r'/home/anubhavk/Desktop/GWAS_Anubhav/Age_quant/test/genotype'
gtffilenames = glob.glob(gtfpath + "/*.imputed.test")
gtoutfilename = [ i.split("/")[8] for i in gtffilenames] 

# In[114]:

# Read the samplefile
samples = list()
with open(samplefilename, 'r') as samfile:
    next(samfile)
    next(samfile)
    for mline in samfile:
        mlinesplit = mline.split()
        this_sample = SampleInfo(sid = mlinesplit[0],
                                 sex = int(mlinesplit[4]), # careful, assuming no NA
                                 pheno = int(mlinesplit[5])) # careful, assuming no NA 
        samples.append(this_sample)
samples


# In[115]:

#nsample = len(samples)
# Read the full genotype
#ncol = nsample * 3 + 5
#dosage = np.loadtxt(gtfilename, usecols=range(5, ncol)).T

'''dosage = list()
snpinfo = list()
# Read the SNP information
with open(gtfilename, 'r') as genfile:
    for mline in genfile:
        mlinesplit = mline.split()
        this_snp = SnpInfo(rsid = mlinesplit[1],
                           bp_location = mlinesplit[2],
                           alt_allele = mlinesplit[3],
                           ref_allele = mlinesplit[4])
        this_snp_dosage = np.array(mlinesplit[5:], dtype=float)
        snpinfo.append(this_snp)
        dosage.append(this_snp_dosage)
dosage = np.vstack(dosage).T
dosage'''


# In[54]:

# Read the age-of-interview
ageinfo = collections.defaultdict(lambda:0)
sexinfo = collections.defaultdict(lambda:0)
samples_without_ages = list()
samples_without_sex = list()
with open(phenofilename, 'r') as phenofile:
    next(phenofile)
    for mline in phenofile:
        mlinesplit = mline.split()
        sid = mlinesplit[0].strip('"')
        age = mlinesplit[6]
        gender = mlinesplit[4]
        sex = mlinesplit[3]
        mage = -1
        msex = -1
        if age == 'NA':
            samples_without_ages.append(sid)
        else:
            if gender == 'NA' and sex == 'NA':
                samples_without_sex.append(sid)
            else:
                if gender == 'NA':
                    msex = 2 if sex == 0 else 1
                else:
                    msex = gender
                mage = int(age) # we have checked age is not NA
        ageinfo[sid] = mage
        sexinfo[sid] = msex


# In[139]:

# Calculate age at onset
with open(expctfilename, 'r') as mfile:
    mline = list(csv.reader(mfile, delimiter='\t'))

expct_age_dict = {}
for x in mline:
    expct_age_dict[int(x[0])] = float(x[1])

# expct_age_dict_female
    
with open(expctfilename_female, 'r') as ffile:
    fline = list(csv.reader(ffile, delimiter='\t'))

expct_age_dict_female = {}
for x in fline:
    expct_age_dict_female[int(x[0])] = float(x[1])
    


# In[85]:

# Write new sample file
with open(sampleoutfilename, 'w') as mfile:
    mfile.write("ID_1 ID_2 missing father mother sex pheno age-at-onset\n")
    mfile.write("0 0 0 D D D B P\n")
    for sample in samples:
        this_sid = sample.sid
        if ageinfo[this_sid] > 0:
            if sample.pheno == 1:
                newpheno = ageinfo[this_sid]
            else:
                if sample.sex == 1:
                    newpheno = expct_age_dict[ageinfo[this_sid]]
                else:
                    newpheno = expct_age_dict_female[ageinfo[this_sid]]
            mstring = '{0} {0} {1} {1} {1} {2} {3} {4}\n'.format(this_sid, 0, sample.sex, sample.pheno, newpheno)
            mfile.write(mstring)


# In[124]:

def genofile_write(dosage,ageinfo, samples, snpinfo,j):
    #all_geno_index = [i for i in range(nsample * 3)]
    selected_samples = [range(i*3, (i+1)*3) for i, sample in enumerate(samples) if ageinfo[sample.sid] > 0]
    selected_samples = np.array(selected_samples)
    selected_samples = selected_samples.reshape(selected_samples.shape[0]*3, )
    dosage_qc = dosage[selected_samples,:]
    dosage_qc_t = dosage_qc.T

    with open(gtoutfilename[j], 'w') as mfile:
        for i, snp in enumerate(snpinfo):
            snpmetainfo = '--- {0} {1} {2} {3}'.format(snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele)
            gtinfo = ' '.join(['%g'% x for x in list(dosage_qc_t[i,:])])
            mfile.write('%s %s\n' % (snpmetainfo, gtinfo))
    return(None)


for j, gtfilename in enumerate(gtffilenames):
    dosage = list()
    snpinfo = list()
    with open(gtfilename, 'r') as genfile:
        for mline in genfile:
            mlinesplit = mline.split()
            this_snp = SnpInfo(rsid = mlinesplit[1],
                               bp_location = mlinesplit[2],
                               alt_allele = mlinesplit[3],
                               ref_allele = mlinesplit[4])
            this_snp_dosage = np.array(mlinesplit[5:], dtype=float)
            snpinfo.append(this_snp)
            dosage.append(this_snp_dosage)
    dosage = np.vstack(dosage).T
    genofile_write(dosage,ageinfo, samples, snpinfo,j)





