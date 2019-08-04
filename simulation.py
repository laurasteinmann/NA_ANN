#Bibliotheken
import numpy as np
import pandas as pd
import random
#Funktionen
#Funktion die Genotyp-Matrix mit zufälligen 20 Markern und 200 Accessions erstellt, ebenso wie die Auswahl der Marker mit Allelefrequency
def set_geno(Input, Input_filter, number):
    marker = np.zeros(number)
    for j in np.arange(0, number) :
        i = np.random.choice(Input_filter.shape[1],1)
        x = np.sum(Input.iloc[:,i])[0]      #Auswahl von Markern
        while (x < 50) | (x > 150):         #Filter der Allelfrequency
            i = np.random.choice(Input_filter.shape[1],1)
            x = np.sum(Input.iloc[:,i])[0]
        marker[j] = i
    Input = Input.transpose()
    marker = marker.astype("Int64")
    genotyp2 = Input.iloc[marker]
    genotyp2 = genotyp2.transpose()
    genotyp2 = genotyp2.astype("Int64") #Dieser Datentyp wird für GBLUP und Netzwerke benötigt
    genotyp2 = Input.iloc[marker]
    return genotyp2
#Skript
FT10 = pd.read_csv("../BachelorThesis/Data/imputed_full/genos_f/FT10.csv", index_col=0)
pheno = pd.read_csv("../BachelorThesis/Data/imputed_full/phenos_f/FT10.csv", index_col=0)
index = pheno.iloc[:,0]
FT10.index = index
lat = pd.read_csv("../BachelorThesis/Data/Simulation/Latitude.csv", index_col=0)
accessions_all = lat.iloc[:,0]
filter = np.isin(accessions_all,index)
possible_acc = accessions_all[filter]
x = np.arange(0, 100)
x = x.astype("str")
for k in x:
    sample = np.random.choice(possible_acc,300) #Auswahl von Accessions
    sample = np.unique(sample)
    sample = sample[:200]
    sample = np.asarray(sample)
    sample = np.astype("int")
    genos = FT10.loc[sample]
    genotyp = set_geno(Input=genos, Input_filter=genos, number=20)  #Auswahl der Marker
    del genotyp.index.name
    genotyp = genotyp.astype("Int64")
    name = "S" + k + ".csv"
    path = "../BachelorThesis/Data/Simulation/genos_f/" + name
    genotyp.to_csv(path)
