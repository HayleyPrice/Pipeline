#!/usr/bin/env bash

# srun -c 40 --pty bash

# Provide file name and dataset as args 1/2
# .txt file with line above heading in format  0\t 0\t 0\t 1\t 1\t 1\t 2\t 2\t 2\t
# protein ids and all columns with non-quantitative data '0'
#    Replicates should be grouped numerically.
# Protein ID format 'sp|O75822|EIF3J_HUMAN'

# # # Performs normalisation of raw data
# echo '################### Reading file: ' $1 ' and normalising data'
# cd Normalisation
# Rscript Normalyser_Manual.R $1 $2
 
# # runs R scripts to 
# 1. Perform t-test analysis on normalised data 
# 2. Separate each of the normalised data into subsets ready for QPROT model analysis
# cd Normalisation
# norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
    # echo '################### Perfoming t-Test: ' $n
    # Rscript ServerT.r $n $2
    # #Rscript QPROTmodel_dataSplit_Server.R $n $2 25 &
# done
# wait

# # # # Performs DE analysis on each data subset using QPROT model
# path="/mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test/"$2"_subs"
# cd $path
# for f in *csv
# do
	# echo '################### Performing QPROT: ' $f
	# Rscript /mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test/subServer.r $f &
# done
# wait

# # # Recombines the subset results and calculates adjusted p value from z score
# cd $path
# cd out
# echo '################### Combing output'
# Rscript QPROTmodel_dataComb.R

# # runs R script to separate each of the normalisaed data into subsets
cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T
norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
for n in ${norms[@]};
do
	echo '################### Processing normalisation: ' $n
	Rscript T_PAnoKegg.R $n $2 &
done
wait

# # # runs R script to separate each of the normalisaed data into subsets
# cd /mnt/hc-storage/users/hprice/Pipeline/FDR
# norms=(Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
	# echo '################### Processing normalisation: ' $n
	# Rscript QPROTmodel_FDR.R PXD004682 $n &
# done
# wait

# # # runs R script to separate each of the normalisaed data into subsets
# cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis
# norms=(Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
	# echo '################### Processing normalisation: ' $n
	# Rscript QPROTmodel_PAnoKegg.R $ &
# done
# wait

# # # Evaluates the best  result
# #echo '################### Analysing results'
# cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T/results
# Rscript T_results.r $2

# # Evaluates the best  result
# echo '################### Analysing results'
# cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/results
# Rscript ResultsSummary.r $2

# echo '#################'