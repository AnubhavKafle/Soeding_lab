
# coding: utf-8

# In[2]:


# All files

#Fitting age - Cardiovascular relation as an exponential function
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from scipy import stats, optimize
from pylab import rc
import math
#import seaborn as sns


# In[3]:


# Defining function for incidence of CVD using population data
###############################################################################################
def y_male_age_incidence(age):
    #incidence = -1.782385 - (-0.02748726/-0.03572809)*(1 - np.exp(+0.03572809*age))
    incidence_power = 0.000001780169*math.pow(age,3.566421)
    return incidence_power

def y_female_age_incidence(age):
    incidence = -0.3940427 - (-0.003166825/-0.06495186)*(1 - np.exp(+0.06495186*age))
    return incidence

def likelihood_func():
    phenotype = 1
    return phenotype
###############################################################################################

random_age_male = np.array(range(17,100,5))
incidence=[]
for i in random_age_male:
    incidence.append(y_male_age_incidence(i))

fig_fit=plt.figure()
fig_fit = plt.plot(random_age_male, incidence)
fig_fit = plt.show()

CVD_age_incidence_normalized_male = pd.DataFrame({"Age":random_age_male,"Incidence(per 1000)":incidence})
#print(CVD_age_incidence_normalized_male)
population_age = random_age_male
population_size = np.array([14,15.9,17,18,18.5,19,19.5,18.8,17,15.0,12.5,10,8,5,3,1,0.3])
EU_population_ageWise = pd.DataFrame({"Age":population_age,"Size(X 10*6)":population_size})
Unnormalized_incidence = 10**3*population_size*incidence
CVD_age_incidence_unnormalized_male_PRIOR = pd.DataFrame({"Age":random_age_male,"Incidence":Unnormalized_incidence})

#print(sum(population_size)) # 168 Millions male population of europe

############## Distribution plot for priori################################# 
fig = plt.figure()
fig = plt.plot(CVD_age_incidence_unnormalized_male_PRIOR["Age"],CVD_age_incidence_unnormalized_male_PRIOR["Incidence"])
fig = plt.show()
############################################################################
############################ Posterior Probability distribution ###########

#posterior   = []

############################### Importing GERMifs files #####################################
#GermI = 


print(CVD_age_incidence_unnormalized_male_PRIOR)


# In[8]:


################# smoothing the incidence distribution graph using spline ############################
from scipy import interpolate
age=random_age_male
incidence = Unnormalized_incidence
age_new = np.arange(17,98,1)
print(age_new)
tck = interpolate.splrep(age, incidence, s=0)
Incidence_new = interpolate.splev(age_new, tck, der=0)
print(len(Incidence_new), len(age), len(age_new))
interpolated = plt.figure()
interpolated = plt.plot(age_new,Incidence_new)
plt.title("Male CVD incidence")
plt.savefig("/home/anubhavk/Desktop/Male_CVD_Incidence.pdf")
interpolated = plt.show()

posterior_distribution = pd.DataFrame({"Age":age_new, "Incidence":Incidence_new})
p=posterior_distribution.loc[posterior_distribution["Age"] >= 17]
#######################################################################################################


# In[4]:


##################### Function for calculating area under the curve using simpson method ##################

def integrate(y_vals, h):
    i=1
    total=y_vals[0]+y_vals[-1]
    for y in y_vals[1:-1]:
        if i%2 == 0:
            total+=2*y
        else:
            total+=4*y
        i+=1
    return total*(h/3.0)
##############################################################################################################


# In[5]:


############# Function for calculating onset age expected value from age at interview ########################

def expected_onset_age(age_at_interview):
    Expected_age={}
    for i in age_at_interview:
        if i not in Expected_age:
            Trimmed_distribution = posterior_distribution.loc[posterior_distribution["Age"] >= i]
            trimmed_age = np.array(Trimmed_distribution["Age"])
            trimmed_incidence = np.array(Trimmed_distribution["Incidence"])
            area = integrate(trimmed_incidence,1)
            probability = trimmed_incidence/area
            #print(sum(probability))
            expected_onset = sum(trimmed_age * probability)
            Expected_age[i]=expected_onset
        
        else: 
            pass
    return(Expected_age)
##############################################################################################################


# In[6]:


banskt_age_at_interview = [x for x in range(17, 98)]
test_expct = expected_onset_age(banskt_age_at_interview)
with open('New_expected_age_at_onset.dat', 'w') as mfile:
    for key, val in test_expct.items():
        mfile.write("%i\t%g\n" % (key, val))

