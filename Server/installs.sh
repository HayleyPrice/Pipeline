#!/usr/bin/env bash
# runs R scripts to intall R packages, Rstan and to compile the Stan model
## please provide path to pipeline as arugument

cd setUp
echo '################### installing dependancies'
Rscript dependancies.r
echo '################### installing stan'
Rscript rStanInstall.r
echo '################### compiling stan model'
Rscript compileStan.r $1

wait

echo '################### Complete'

