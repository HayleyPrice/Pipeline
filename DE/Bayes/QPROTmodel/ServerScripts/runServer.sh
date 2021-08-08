# sets directory
#cd Z:/hprice/Rbatch
# runs R script to separate data into subsets
Rscript Ql_dataSplit_Server.R PXD001385_P_PeptideIon_HiN_NG_23.09.16_sum_sumCharge_QL_0 PXD001385 0 2
#
#Rscript subServer.R subset1 PXD001385 0 
#
cd PXD001385\0\subsets
#
for f in *.csv
do
	echo "Processing $f file"
	Rscript subServer.R $f PXD001385 0 
done
#	
#
#setlocal ENABLEDELAYEDEXPANSION
# loops through the  subsets, opens cmd and calls R script to perform sampling
#for /l %%x in (1, 1, 2) do (
#	
#	set ctr=%%x
#	set inFile=subset!ctr!
#	set outFile=out_subset!ctr!
#	start cmd /k Rscript subServer.R !inFile! !outFile! PXD001385 0
#)
#sleep 