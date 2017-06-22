
# coding: utf-8

# In[2]:


#Fitting age - Cardiovascular relation as an exponential function for female
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from scipy import stats, optimize
from pylab import rc
import math
#import seaborn as sns


# In[3]:


def y_female_age_incidence(age):
    incidence = 1.646609*math.pow(10,-9)*math.pow(age,5.099324) #0.000001780169*math.pow(age,3.566421)
    #incidence = 43.29204 - (-1.505984/-0.1522527)*(1 - np.exp^(+0.1522527*age))
    return incidence

random_age_female = np.array(range(17,100,5))
incidence=[]

for i in random_age_female:
    incidence.append(y_female_age_incidence(i))
    
CVD_age_incidence_normalized_female = pd.DataFrame({"Age":random_age_female,"Incidence(per 1000)":incidence})
#print(CVD_age_incidence_normalized_male)
population_age = random_age_female
population_size = np.array([13.3,15,15.7,17,17.5,18,19,18,17,16,14,12,11,8.8,5.5,2.8,0.5])
EU_population_ageWise = pd.DataFrame({"Age":population_age,"Size(X 10*6)":population_size})
Unnormalized_incidence = 10**3*population_size*incidence
CVD_age_incidence_unnormalized_female_PRIOR = pd.DataFrame({"Age":random_age_female,"Incidence":Unnormalized_incidence})
print(CVD_age_incidence_unnormalized_female_PRIOR)

########################


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


# In[5]:



fig_fit=plt.figure()
fig_fit = plt.plot(random_age_female, incidence)
fig_fit = plt.show()



fig = plt.figure()
fig = plt.plot(CVD_age_incidence_unnormalized_female_PRIOR["Age"],CVD_age_incidence_unnormalized_female_PRIOR["Incidence"])
fig = plt.show()


# In[12]:


from scipy import interpolate
age=random_age_female
incidence = Unnormalized_incidence
age_new = np.arange(17,98,1)
print(age_new)
tck = interpolate.splrep(age, incidence, s=0)
Incidence_new = interpolate.splev(age_new, tck, der=0)
print(len(Incidence_new), len(age), len(age_new))
interpolated = plt.figure()
interpolated = plt.plot(age_new,Incidence_new)
plt.title("Female CVD incidence")
plt.savefig("/home/anubhavk/Desktop/Female_CVD_Incidence.pdf")
#interpolated = plt.show()
posterior_distribution = pd.DataFrame({"Age":age_new, "Incidence":Incidence_new})
p=posterior_distribution.loc[posterior_distribution["Age"] >= 17]
#######################################################################################################


# In[17]:


banskt_age_at_interview = [x for x in range(17, 98)]
test_expct = expected_onset_age(banskt_age_at_interview)
with open('New_expected_age_at_onset_female.dat', 'w') as mfile:
    for key, val in test_expct.items():
        mfile.write("%i\t%g\n" % (key, val))


# In[33]:


############### Function to take distribution information to create pseudo age of onset from age at interview ##

arbitrary_onset_age = np.array([50,60,70,80,90]) #Random sampling from the distribution

#dict_age_onset_age_weight = {}

def onset_age_distribtuion(age_at_interview):
    Trimmed_distribution = posterior_distribution.loc[posterior_distribution["Age"] >= age_at_interview]
    trimmed_age = np.array(Trimmed_distribution["Age"])
    trimmed_incidence = np.array(Trimmed_distribution["Incidence"])
    area = integrate(trimmed_incidence,1)
    probability = trimmed_incidence/area
    #print(probability)
    PDF = pd.DataFrame({"Age":trimmed_age, "weight":probability})
    #print(PDF)
    #print(Trimmed_distribution)
    return(PDF)

def dic_expectedOnsetAges(age_at_interview):
    dict_age_onset_age_weight = {}
    for i in age_at_interview: #32 is the least age in the distribution and it grows in 1 interval
        PDF = onset_age_distribtuion(i)
        weights=[]
        age_trimmed = []
        #ranges=[age_at_interview[0],65,75,85,age_at_interview[-1]]
        p=i
        for j in arbitrary_onset_age:
            if j >= i:
                age_trimmed.append(j)
                area = integrate(np.array(PDF.loc[(PDF["Age"] >= p) & (PDF["Age"]<= j+5 ) ]["weight"]),1)
                p=j+5
                weights.append(area)
            else:
                pass
        #print(age_trimmed)
        distribution_points = pd.DataFrame({"onset_age":age_trimmed,"weight":weights})
        if i not in dict_age_onset_age_weight:
            dict_age_onset_age_weight[i] = distribution_points
        else:
            pass
    return(dict_age_onset_age_weight)


# In[49]:


banskt_age_at_interview = [x for x in range(32, 98)]
test_expct = dic_expectedOnsetAges(banskt_age_at_interview)
with open('quantised_age_at_onset_female.dat', 'w') as mfile:
    for key,val in test_expct.items():
        for j in range(val.shape[0]):
            mfile.write("%i\t%g\t%g\n" % (key, val.onset_age[j],val.weight[j] ))


# In[48]:


p = test_expct[32]
#print(p)
#print(p.onset_age[1])

for i in range(p.shape[0]):
    print("%g\t%g\n" %(p.onset_age[i],p.weight[i]))
    
    
for key,val in test_expct.items():
    for j in range(val.shape[0]):
        print("%i\t%g\t%g\n" % (key, val.onset_age[j],val.weight[j] ))
