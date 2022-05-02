#!/usr/bin/env bash

# # Performs normalisation of raw data
# echo '################### Normalising data'
# cd Normalisation
# Rscript Normalyser_Manual.R
# 
# # # runs R script to separate each of the normalisaed data into subsets
# for n in *normalized_Log2
# do
# 	echo '################### Processing file: ' $n
# 	Rscript QPROTmodel_dataSplit_Server.R $n PXD004682 25 &
# done
# wait
# 
# # # Performs DE analysis on each data subset using QPROT model
# cd /mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test
# for f in *.csv
# do
# 	echo '################### Processing file: ' $f
# 	Rscript subServer.r $f PXD004682 &
# done
# wait
# 
# # Recombines the subset results and calculates adjusted p value from z score
# cd out
# echo '################### Combing output'
# Rscript QPROTmodel_dataComb.R

# # # FDRs calculated
# #norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# #for n in ${norms[@]};
# do
	# echo '################### Calculating FDR: ' $n
	# Rscript QPROTtest_FDR.R $1 $2 &
# done
# wait

# runs R script to separate each of the normalisaed data into subsets
#cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis
thresholds=(0.000001 0.000002 0.000003 0.000004 0.000005 0.000006 0.000007 0.000008 0.000009 \
0.00001 0.00002 0.00003 0.00004 0.00005 0.00006 0.00007 0.00008 0.00009 \
0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 \
0.001  0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 \
0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)

for t in ${thresholds[@]};
do
	echo '################### Processing normalisation: ' $t
	Rscript QPROTtest_PAnoKegg.R $1 $2 $t &
done
wait

# # Evaluates the best  result
# #echo '################### Analysing results'
# #cd $2
# Rscript QPROTmodel_resultsQPROT.r $2

# echo '################### Complete'_QPROT

#srun -c 40 --pty bash