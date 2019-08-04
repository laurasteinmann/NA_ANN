#Bibliographies
import numpy as np
import pandas as pd
import os
#Design Filter for filtering vcf_file
ordner = "../BachelorThesis/Data/imputed_full/genos_f/" #Imputed-Data
ordner2 = "../BachelorThesis/Data/Filter_pos/" #Speicherort
dateien = os.listdir(ordner) #Liste von allen Phänotypen
z = np.arange(0,146)  #Anwendung auf alle Phänotypen
for i in z:
    load = ordner + dateien[i]
    save = ordner2 + "pos_" + dateien[i]
    geno = pd.read_csv(load, index_col=0) #Einlesen der Daten 
    positions = list(geno.columns.values) #Extraktion von Markerpositionen
    marker = pd.DataFrame(positions, columns=["chromosom"]) #Erstellung der neuen Filtermatrix
    marker["position"] = positions
    marker.iloc[:,0] = marker.iloc[:,0].map(lambda x: str(x)[:1]) #Chromsomnummer
    marker.iloc[:,1] = marker.iloc[:,1].map(lambda x: str(x)[2:]) # Position auf dem Chromosom
    marker.to_csv(save, header=False, index=False, sep="\t")
