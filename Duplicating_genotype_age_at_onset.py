
# coding: utf-8

# In[ ]:

import numpy as np
import collections
import csv
import glob
SAMPLE_FIELDS = ['sid', 'old_id', 'sex', 'pheno', 'weight']
class SampleInfo(collections.namedtuple('_SampleInfo', SAMPLE_FIELDS)):
    __slots__ = ()
    
SNPINFO_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

home = os.getenv("HOME")
samplefilename = '/scratch/akaphle/Pheno/Gxxx/Gxxx_QC.sample'
phenofilename = '/scratch/akaphle/clinical/Gxxx/Gxxx_GerMIFS_formated_pheno'
expctfilename_male = home+'/anubhav_GWAS/scripts/quantised_age_at_onset_male.dat'
expctfilename_female = home+'/anubhav_GWAS/scripts/quantised_age_at_onset_female.dat'
sampleoutfilename = '/scratch/akaphle/New_pheno_duplicated/Gxxx/Gxxx_age_duplicated_sample'
gtfpath = r'/scratch/akaphle/Genotype/Gxxx'
gtffilename = glob.glob(gtfpath + "/*.imputed")
gtoutfilename = [i.split("/")[5] for i in gtffilename]


# In[ ]:

# Read the samplefile
samples = list()
with open(samplefilename, 'r') as samfile:
    next(samfile)
    next(samfile)
    for mline in samfile:
        mlinesplit = mline.split()
        this_sample = SampleInfo(sid = mlinesplit[0],
                                 old_id = '',
                                 sex = int(mlinesplit[4]), # careful, assuming no NA
                                 pheno = int(mlinesplit[5]),
                                 weight = 1) # careful, assuming no NA 
        samples.append(this_sample)


# In[ ]:

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


# In[ ]:

# Calculate age at onset
with open(expctfilename_male, 'r') as mfile:
    mline = list(csv.reader(mfile, delimiter='\t'))

expct_age_dict_male = {}
for x in mline:
    if int(x[0]) in expct_age_dict_male:
        expct_age_dict_male[int(x[0])].append([float(x[1]),float(x[2])])
    else:
        expct_age_dict_male[int(x[0])] = [[float(x[1]),float(x[2])]]
# expct_age_dict_female
    
with open(expctfilename_female, 'r') as ffile:
    fline = list(csv.reader(ffile, delimiter='\t'))

expct_age_dict_female = {}
for x in fline:
    if int(x[0]) in expct_age_dict_female:
    	expct_age_dict_female[int(x[0])].append([float(x[1]),float(x[2])])
    else:
        expct_age_dict_female[int(x[0])] = [[float(x[1]),float(x[2])]]


# In[ ]:

# Create new samples
newsamples = list()
for ii, sample in enumerate(samples):
    this_sid = sample.sid
    if ageinfo[this_sid] > 0:
        if sample.pheno == 1:
            this_sample = SampleInfo(sid = this_sid,
                                     old_id = this_sid,
                                     sex = sample.sex,
                                     pheno = ageinfo[this_sid],
                                     weight = 1)
            newsamples.append(this_sample)
        else:
            if sample.sex == 1:
                newphenos = expct_age_dict_male[ageinfo[this_sid]]
                nnew = len(newphenos)
                for inew in range(nnew):
                    new_sid = '%s_%02i' % (this_sid, inew)
                    this_sample = SampleInfo(sid = new_sid,
                                             old_id = this_sid,
                                             sex = sample.sex,
                                             pheno = newphenos[inew][0],
                                             weight = newphenos[inew][1])
                    newsamples.append(this_sample)
            elif sample.sex == 2:
                newphenos = expct_age_dict_female[ageinfo[this_sid]]
                nnew = len(newphenos)
                for inew in range(nnew):
                    new_sid = '%s_%02i' % (this_sid, inew)
                    this_sample = SampleInfo(sid = new_sid,
                                             old_id = this_sid,
                                             sex = sample.sex,
                                             pheno = newphenos[inew][0],
                                             weight = newphenos[inew][1])
                    newsamples.append(this_sample)


# In[ ]:

# Write new sample file
with open(sampleoutfilename, 'w') as mfile:
    mfile.write("ID_1 ID_2 missing father mother sex age_at_onset weight\n")
    mfile.write("0 0 0 D D D P P\n")
    for sample in newsamples:
            mstring = '{0} {5} {1} {1} {1} {2} {3} {4} \n'.format(sample.sid, 0, sample.sex, sample.pheno, sample.weight, sample.old_id)
            mfile.write(mstring)


# In[ ]:
# Create new genotype
def genofile_write(dosage, ageinfo, newsamples, snpinfo,p):
    newdosage = list()
    file_name = '/scratch/akaphle/New_Genotype_duplicated/Gxxx/' + gtoutfilename[p]+'_duplicated'
    with open(file_name, 'w') as mfile:
        for ii, snp in enumerate(snpinfo):
            this_dosage = list()
            for sample in newsamples:
                oid = sample.old_id
                sample_dosage = dosage[(oid * 3) : ((oid + 1) * 3), ii]
                this_dosage.append(sample_dosage)
            snpmetainfo = '--- {0} {1} {2} {3}'.format(snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele)
            gtinfo = ' '.join(['%g'% x for x in this_dosage])
            mfile.write('%s %s\n' % (snpmetainfo, gtinfo))
            return(None)

# Read the SNP information
#for i, gtfilename in enumerate(gtffilenames):
#    dosage = list()
#    snpinfo = list()
#    with open(gtfilename, 'r') as genfile:
#        for mline in genfile:
#            mlinesplit = mline.split()
#            this_snp = SnpInfo(rsid = mlinesplit[1],
#                               bp_location = mlinesplit[2],
#                               alt_allele = mlinesplit[3],
#                               ref_allele = mlinesplit[4])
#            this_snp_dosage = np.array(mlinesplit[5:], dtype=float)
#            snpinfo.append(this_snp)
#            dosage.append(this_snp_dosage)
#    genofile_write(dosage,ageinfo, newsamples, snpinfo,i)

