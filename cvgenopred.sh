#!/bin/bash 
#echo   $SLURM_ARRAY_TASK_ID

i=$SLURM_ARRAY_TASK_ID
echo $i ;
A=($(cat phenos)) ;

phe=$(echo "${A[@]:$i:1}")

 
g='~/Simulation/genos_f/'$phe'.csv' 
p='~/Simulation/phenos_sim85/'$phe'.csv' 
c='~/Simulation/CV/'$phe'.csv' 


Rscript ex.geno.pred.r -x $g -y $p -cv $c -phe $phe -label Phenotype_values
