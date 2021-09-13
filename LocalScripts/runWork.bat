
:: sets path to R
PATH C:\Program Files\R\R-4.1.0\bin\x64;%path%

:: sets directory
cd E:\OneDrive\PhD\Project\Thesis\4_Pipeline\Pipeline\DE\Bayes\QPROTmodel\LocalScripts

:: runs R script to separate data into subsets
Rscript QPROTmodel_dataSplit_Local.R PXD004682_QPROTin_test PXD004682 32

setlocal ENABLEDELAYEDEXPANSION

:: loops through the 1st 25 subsets, opens cmd and calls R script to perform sampling
for /l %%x in (1, 1, 24) do (
	
	set ctr=%%x
	set inFile=subset!ctr!
	set outFile=out_subset!ctr!
	start cmd /k Rscript subLocal.R !inFile! !outFile! PXD004682
)

:: loops through the 2nd 25 subsets, opens cmd and calls R script to perform sampling
for /l %%x in (25, 1, 32) do (
	
	set ctr=%%x
	set inFile=subset!ctr!
	set outFile=out_subset!ctr!
	start cmd /k Rscript subLocal.R !inFile! !outFile! PXD004682
)


pause 