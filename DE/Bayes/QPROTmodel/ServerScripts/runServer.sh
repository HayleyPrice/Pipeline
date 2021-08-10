# sets directory
#cd Z:/hprice/Rbatch
# runs R script to separate data into subsets
Rscript QPROTmodel_dataSplit_Server.R PXD004682_QPROTin_test PXD004682 32

cd test
#
for f in *.csv
do
	echo "Processing $f file"
	Rscript subServer.R $f PXD004682
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