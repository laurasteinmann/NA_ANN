#Bibliotheken
import numpy as np
import pandas as pd

#Data Setting
dateien = ("Seed_bank_133-91", "Seed_Dormancy", "Seedling_Growth", "Silique_16", "Silique_22",
	   "SizeLocSweden2009", "SizeMainEffect2009", "SizePlantingSummer2009", "SizePlantingSummerLocSweden2009", "Size_spain_2009_1st_experiment", 
	   "Size_spain_2009_2nd_experiment", "Size_sweden_2009_1st_experiment", "Size_sweden_2009_2nd_experiment", "Storage_28_days", "Storage_56_days",
	   "Storage_7_days", "Trichome_avg_C", "Trichome_avg_JA", "Vern_Growth", "Width_10", "Width_16", "Width_22", "YEL_", "YieldLocSweden2009", 
	   "YieldMainEffect2009", "YieldPlantingSummer2009", "YieldPlantingSummerLocSweden2009", "Yield_spain_2009_1st_experiment",
	   "Yield_spain_2009_2nd_experiment", "Yield_sweden_2009_1st_experiment", "Yield_sweden_2009_2nd_experiment", "Zn66" )

z = np.arange(0,32)
##Make cv file and cut xset as phenotype accessions
#Set x and phenotype accessions
for i in z:
    datei = "../BachelorThesis/Data/genos_f/pos_" + dateien[i] + ".csv"
    genos = pd.read_csv(datei, sep="\t", index_col=0)
    datei2 = "../BachelorThesis/Data/imputed_full/phenos_f/" + dateien[i] + ".csv"
    phenos = pd.read_csv(datei2, index_col=0)
    phenos2 = pd.DataFrame(columns=["accession_id","phenotype_value"])
    a = phenos.iloc[:,0]
    phenos2.iloc[:,0] = a
    b = phenos.iloc[:,1]
    phenos2.iloc[:,1] = b
    c = genos.index.values
    d = np.isin(c,a)
    genos2 = genos[d]
    e = genos2.index.values
    f = np.isin(a,e)
    phenos2 = phenos2[f]
    genos2 = genos2.astype("Int64")
    save1 = "../BachelorThesis/Data/phenos_y/" + dateien[i] + ".csv"
    phenos2.to_csv(save1,index=False)
    save2 = "../BachelorThesis/Data/genos_x/" + dateien[i] + ".csv"
    genos2.to_csv(save2)
#CV_File
    cv = pd.DataFrame(columns=np.arange(1,51), index=e)
    one =int(cv.shape[0]*0.2)
    null = int(cv.shape[0] - one)
    sample = np.array([0] * null + [1] * one)
    x = np.arange(0,50)
    for j in x:
        np.random.shuffle(sample)
        cv.iloc[:,j] = sample
    z = np.arange(1,51)
    names=np.array([])
    for j in x:
        name = "cv_" + str(z[j])
        names = np.append(names,name)
    cv.columns = names
    save3 = "../BachelorThesis/Data/CV/" + dateien[i] + ".csv"
    cv.to_csv(save3)
#Imputed x
    datei3 = "../BachelorThesis/Data/imputed_full/genos_f/" + dateien[i] + ".csv"
    imputed = pd.read_csv(datei3, index_col=0)
    datei4 = "../BachelorThesis/Data/phenos_y/" + dateien[i] + ".csv"
    pheno = pd.read_csv(datei4)
    datei5 = "../BachelorThesis/Data/imputed_full/phenos_f/" + dateien[i] + ".csv"
    imputed2 = pd.read_csv(datei5, index_col=0)
    a = imputed2.iloc[:,0]
    a = a.values
    imputed.index = a
    b = pheno.iloc[:,0]
    c = np.isin(a,b)
    imputed3 = imputed[c]
    imputed3 = imputed3.astype("int64")
    imputed4 = imputed3.transpose()
    d = np.asarray(imputed3.sum(axis=0)!=0)
    imputed5 = imputed4[d]
    imputed6 = imputed5.transpose()
    save4 = "../BachelorThesis/Data/imputed_x/" + dateien[i] + ".csv"
    imputed6.to_csv(save4)






