source("/storage/evolgen/R/scripts/GWAS/emma.r")
source("/storage/evolgen/R/scripts/GWAS/gwas.r")
load("GWAS/K_2029.rda") #Kann durch eine eigene emma.kinship.matrix ersetzt werden
results = []
for (i  in 1:100){
    datei = "S" + i + ".csv"
    pfadx = "Simulation/genos_f/" + datei
    pfady = "Simulation/phenos_lat80" + datei
    X <- read.csv(pfadx, row.names=1) #Genotyp-matrix
    Y <- read.csv(pfady)  #Phenotyp-matrix
    h <- amm_gwas(Y, X, K, run=FALSE) #Berechnung der heritability
    results[i] = h 
  }
write.csv(results, "results_lat90.csv")
