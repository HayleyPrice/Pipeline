# sets directory
#cd Z:/hprice/Rbatch

# # Performs normalisation of raw data
# echo '################### Normalising data'
cd Normalisation
# Rscript Normalyser_Manual.R
 
# # # runs R scripts to 
# 1. Perform t-test analysis on normalised data 
# 2. Separate each of the normalised data into subsets ready for QPROT modle analysis
for n in *normalized_Log2
do
    echo '################### Processing file: ' $n
    Rscript ServerT.r $n
    Rscript QPROTmodel_dataSplit_Server.R $n PXD004682 25 &
done
wait

# # # Performs DE analysis on each data subset using QPROT model
# cd /mnt/hc-storage/users/hprice/Pipeline/QPROTmodel/test
# for f in *csv
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

# runs R script to separate each of the normalisaed data into subsets
cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T
norms=(None Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
for n in ${norms[@]};
do
	echo '################### Processing normalisation: ' $n
	Rscript T_PAnoKegg.R PXD004682 $n &
done
wait

# # runs R script to separate each of the normalisaed data into subsets
# cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis
# norms=(Loess-G RLR-G VSN-G TI-G MedI-G AI-G Quantile )
# for n in ${norms[@]};
# do
# 	echo '################### Processing normalisation: ' $n
# 	Rscript QPROTmodel_PAnoKegg.R PXD004682 $n &
# done
# wait
# Evaluates the best  result
#echo '################### Analysing results'
cd /mnt/hc-storage/users/hprice/Pipeline/PathwayAnalysis/T/results
Rscript T_results.r PXD004682

# # Evaluates the best  result
# #echo '################### Analysing results'
# cd results
# Rscript QPROTmodel_results.r PXD004682

echo '################### Complete'

#srun -c 40 --pty bash