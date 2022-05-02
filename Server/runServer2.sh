#!/usr/bin/env bash

# srun -c 40 --pty bash

# Provide file name and dataset as args 1/2
# .txt file with line above heading in format  0\t 0\t 0\t 1\t 1\t 1\t 2\t 2\t 2\t
# protein ids and all columns with non-quantitative data '0'
#    Replicates should be grouped numerically.
# Protein ID format 'sp|O75822|EIF3J_HUMAN'

# # # # Performs normalisation of raw data
# echo '################### Reading file: ' $1 ' and normalising data'
# Rscript Normalyser_Manual.R $1 $2
 
# # runs R scripts to 
# #1. Perform t-test analysis on normalised data 
# #2. Separate each of the normalised data into subsets ready for QPROT model analysis
# norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
    # echo '################### Perfoming t-Test: ' $n
    # Rscript ServerT.r $n $2
    # Rscript QPROTmodel_dataSplit_Server.R $n $2 25 &
	# echo '################### Data seperated into ' 25 ' subsets'
# done
# wait

# # # # Performs DE analysis on each data subset using QPROT model
# path="/mnt/hc-storage/users/hprice/Pipeline/"$2"/Qmodel/subs"
# cd $path
# for f in *csv
# do
	# echo '################### Performing QPROT: ' $f
	# Rscript /mnt/hc-storage/users/hprice/Pipeline/subServer.r $f $2 &
# done
# wait

# # # Combines the subset results
# cd /mnt/hc-storage/users/hprice/Pipeline
# echo '################### Combing output'
# Rscript QPROTmodel_dataComb.R $2

# Performs pathway analysis on tTest results
norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
for n in ${norms[@]};
do
	echo '################### tTest PA: ' $n
	Rscript T_PAnoKegg.R $n $2 &
done
wait

# # # FDRs calculated
# norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
	# echo '################### Calculating FDR: ' $n
	# Rscript QPROTmodel_FDR.R $n $2 &
# done
# wait

# # Performs pathway analysis on Qmodel results
norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
for n in ${norms[@]};
do
	echo '################### Qmodel PA: ' $n
	Rscript QPROTmodel_PAnoKegg.R $n $2 &
done
wait

# Evaluates the best  result
echo '################### Analysing results'

Rscript ResultsSummary.r $2

echo '#################'