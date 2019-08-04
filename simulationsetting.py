#Bibliotheken
import numpy as np
import pandas as pd
import random
#Funktionen
def cv(Input):
    cv = pd.DataFrame(columns=np.arange(1,51), index=Input.index.values)
    one = int(cv.shape[0]*0.2)
    null = int(cv.shape[0] - one)
    sample = np.array([0] * null + [1] * one)
    x = np.arange(0, 50)
    for j in x:
        np.random.shuffle(sample)
        cv.iloc[:,j] = sample
    z = np.arange(1, 51)
    names = np.array([])
    for j in x:
        name = "cv_" + str(z[j])
        names = np.append(names, name)
    cv.columns = names
    return cv
def addNA(accessions, Marker, x, geno, loop):
    random.seed(loop)
    for j in x:
        colum = random.randint(0, Marker)
        index = random.randint(0, accessions)
        geno.iloc[index, colum] = np.nan
    return(geno)
def sim_phe(Input):
    lat = pd.read_csv("../BachelorThesis/Data/Simulation/Latitude.csv", index_col=0)
    accessions = genos.index.values
    ecotypeid = lat.iloc[:,0]
    ecotypeid = ecotypeid.astype("int")
    filter = np.isin(ecotypeid, accessions)
    lat.index = ecotypeid
    lat2 = lat[filter]
    pheno = pd.DataFrame(index=lat2.index)
    latitude = lat.iloc[:,2]
    pheno["Latitude"] = latitude
    mean = np.mean(latitude)
    sample = np.random.choice(lat2.index, 50)
    for i in sample:
        pheno.loc[i,"Latitude"] = mean
    value = []
    values = Input.sum(axis=1)
    pheno["Sum_NAs"] = values
    pheno["Phenotype_values"] = pheno["Latitude"] + 100 * pheno["Sum_NAs"]
    return pheno
#Skript
x = np.arange(0, 100)
x = x.astype("str")
for i in x:
    name = "S" + i + ".csv"
    path = "../BachelorThesis/Data/Simulation/genos_f/" + name
    genotyp = pd.read_csv(path, index_col=0)
    select = genotyp.columns.values
    index = genotyp.index
    genos = addNA(accessions=199, Marker=19, x = np.arange(0, 800), geno = genotyp, loop = i)
    genos = genos.astype("Int64")
    path2 = "../BachelorThesis/Data/Simulation/genos_sim25/" + name
    genos.to_csv(path2)
    #CV-File
    #cv = cv(Input=genotyp)
    #path3 = "../../Data/Simulation/CV/" + name
    #cv.to_csv(path3)
    #Coding imputed
    genos = genos.replace(1, 0, regex = True)
    genos = genos.fillna(1)
    genos = genos.astype("Int64")
    path5 = "../BachelorThesis/Data/Simulation/imputed_sim25/" + name
    genos.to_csv(path5)
    #Simulierung Phänotyp
    index = genos.index
    pheno = sim_phe(Input=genos)
    path4 = "../BachelorThesis/Data/Simulation/phenos_20/" + name
    pheno.to_csv(path4)
