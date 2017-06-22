def y_male_age_incidence(age):
    #incidence = -1.782385 - (-0.02748726/-0.03572809)*(1 - np.exp(+0.03572809*age))
    incidence_power = 0.000001780169*math.pow(age,3.566421)
    return incidence_power

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

############### Function to take distribution information to create pseudo age of onset from age at interview ##

arbitrary_onset_age = np.array([60,70,80,90]) #Random sampling from the distribution

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
##############################
    
