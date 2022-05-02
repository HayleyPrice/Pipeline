#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop('At least one argument must be supplied (input file).n', call.=FALSE)
}
norm <-args[1]
dataSet <- args[2]
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/Results', sep = ""))
fileIn <- paste(dataSet, "_", norm, "_DEresults.csv", sep = '')

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server/FDR') 
#inFile <- "PXD004682_QPROTout_ProgNorm_2up_ThreshResults.csv"
outFile <- paste(dataSet, "_", norm, "_FDRresults.csv", sep = '')

example_qprot_fdr <- read.csv(fileIn, header = FALSE)
cols <- example_qprot_fdr[2,]
cols <- cols[,-1]
cols <- cbind('Protein', cols)
example_qprot_fdr <- example_qprot_fdr[c(-1,-2),]

# Column names
colnames(example_qprot_fdr) <- cols
#example_qprot_fdr <- example_qprot_fdr[order(example_qprot_fdr$Zstatistic), ]
#example_qprot_density <- read.table('example_qprot_density', header = TRUE, sep = '\t')

np <- nrow(example_qprot_fdr)
example_qprot_fdr <- example_qprot_fdr[order(-abs(as.numeric(example_qprot_fdr$Zstatistic))), ]

log_gaussian_pdf <- function(x, mu, sigmasq) {
    out <- -0.5*(x-mu)^2/sigmasq - 0.5 * log(2 * pi * sigmasq)
    return(out)
}

gaussian_pdf <- function(x, mu, sigmasq) {
    return(exp(log_gaussian_pdf(x, mu, sigmasq)))
}

## initialising parameters
pi_0 <- 0.9
pi_1 <- 0.1
pi_true <- 0.1
nullmean <- rnorm(1, 0 , 0.1)
nullvar <- exp(rnorm(1, 0 , 0.1))

if(np > 10000) {
    ntick <- 10000
} else {
    ntick <- np
}

zstat <- as.numeric(example_qprot_fdr$Zstatistic)
zstat_sort <- sort(zstat)
xmin <- zstat_sort[1] - 0.1
xmax <- zstat_sort[np] + 0.1

tickMarks <- NULL
f0 <- NULL
f1 <- NULL
f <- NULL

for (i in 1:ntick) {
    tickMarks[i] <- xmin + i * (xmax - xmin) / ntick
}

for (i in 1:ntick) {
    f0[i] <- 0
    f0[i] <- f0[i] + gaussian_pdf(tickMarks[i], nullmean, nullvar)
    f1[i] <- 0.5 * gaussian_pdf(tickMarks[i], -10, 1) +  0.5 * gaussian_pdf(tickMarks[i], 10, 1)
    f[i] <- pi_true * f1[i] + (1.0-pi_true) * f0[i]
}

df <- NULL
df0 <- NULL
df1 <- NULL
for (i in 1:np) {
    df[i] <- 0
    df0[i] <- 0
    df1[i] <- 0
}

## Overall nonparametric density for f

#getBandwidth()
##identifying the range to fit null component to the overall density
mm <- nullmean
vv <- nullvar
temp_sd <- sqrt(vv)
bandwidth <- 1.06 * temp_sd/(np^0.2)

denom <- 0
for (j in 1:np) {
    denom <- denom + 1
}

for (i in 1:ntick) {
    numer <- 0
        for (j in 1:np) {
        numer <- numer + gaussian_pdf((zstat[j] - tickMarks[i])/bandwidth, 0.0, 1.0)
    }
    f[i] <- numer / (denom * bandwidth)
}

## Check mean shift
max_y <- 0
tmp_mean <-NULL
for (i in 1:ntick) {
    if(f[i] > max_y) {
        max_y <- f[i]
        tmp_mean <- tickMarks[i]
    }
}
nullmean <- tmp_mean

## Fit gaussian unif
zz <- matrix(nrow = np, ncol = 2)
z_min <- zstat_sort[1] - 0.1
z_max <- zstat_sort[np] + 0.1
old_mean <- NULL
old_var <- NULL

for (k in 1:np) {
    zz[k, 1] <- 0.1
    zz[k, 2] <- 0.9
    
    zz[k, 1] <- zz[k, 1] * exp(rnorm(1, 0 , 0.1))
    zz[k, 2] <- zz[k, 2] * exp(rnorm(1, 0 , 0.1))
    
    tmp_denom <- 0
    
    tmp_denom <- tmp_denom + zz[k, 1]
    tmp_denom <- tmp_denom + zz[k, 2]
    
    zz[k, 1] <- zz[k, 1] / tmp_denom
    zz[k, 2] <- zz[k, 2] / tmp_denom
}

counter <- 0
tol <- 1000
tmpl <- NULL

while(tol > 0.00001 && counter < 100) {
    old_mean <- nullmean
    old_var <- nullvar
    
    tmp_num <- 0
    tmp_denom <- 0
    for (i in 1:np) {
        tmp_num <- tmp_num + (zstat[i] * zz[i, 1])
    }
    for (i in 1:np) {
        tmp_denom <- tmp_denom + (zz[i, 1])
    }
    # if(tmp_denom > 0) {
    #     nullmean <- tmp_num / tmp_denom
    # } else{
    #     nullmean <- 0
    # }
    
    tmp_num <- 0
    tmp_denom <- 0
    for (i in 1:np) {
        tmp_num <- tmp_num + ((zstat[i] - nullmean) ^ 2) * (zz[i, 1])
    }
    for (i in 1:np) {
        tmp_denom <- tmp_denom + (zz[i, 1])
    }
    if(tmp_denom > 0) {
        nullvar <- tmp_num / tmp_denom
    } else{
        nullvar <- 0
    }
    for (i in 1:np) {
        tmpl[1] <- pi_0 * gaussian_pdf(zstat[i], nullmean, nullvar)
        tmpl[2] <- pi_1 / (z_max - z_min)
        tmp_denom <- 0
        tmp_denom <- tmp_denom + tmpl[1]
        tmp_denom <- tmp_denom + tmpl[2]
        if (tmp_denom > 0) {
            zz[i, 1] <- tmpl[1] / tmp_denom
            zz[i, 2] <- tmpl[2] / tmp_denom
        } else {
            zz[i, 1] <- 0
            zz[i, 2] <- 0
        }
    }
    tmp <- 0
    for (i in 1:np) {
        tmp <- tmp + zz[i, 1]
    }
    tmp <- tmp / np
    pi_0 <- tmp
    
    tmp <- 0
    for (i in 1:np) {
        tmp <- tmp + zz[i, 2]
    }
    tmp <- tmp / np
    pi_1 <- tmp
    
    tol <- 0
    tol <- tol + abs(old_mean - nullmean)
    tol <- tol + abs(old_var - nullvar)
    
    counter <- counter + 1
}

## Evaluate f0 and df0
for (i in 1:ntick) {
    tmp <- 0
    tmp <- tmp + (pi_0 * gaussian_pdf(tickMarks[i], nullmean, nullvar))
    tmp <- tmp / (1.0 - pi_1)
    f0[i] <- tmp
}
for(i in 1:np) {
    tmp <- 0
    tmp <- tmp + (pi_0 * gaussian_pdf(zstat[i], nullmean, nullvar))
    tmp <- tmp / (1.0 - pi_1)
    df0[i] <- tmp
}
pi_true <- pi_1

# Re-scale variances
# set variance rescale gradient: 0.99 -> 0.98 -> ... -> 0.01
# loss function: all f0 values should be under f within reasonable areas

mm <- 0
vv <- 0
mm <- mm + nullmean
vv <- vv + nullvar

z_left = mm - sqrt(vv)
z_right =  mm + sqrt(vv)

f0_temp <- NULL
for (i in 1:ntick) {
    f0_temp[i]<- f0[i]
}

f_ord <- sort(f)
tol_diff <- f_ord[np] * 0.01

ct <- 0
ntry <- 10000
rescale_factor <- 1
is_under <- 1
is_over <- 0

while(rescale_factor > 0 && ct < ntry) {
    
    #if(ct > 0 && ct %% 1000 == 0) {
    #    print('error')
    #}
    for(i in 1:ntick) {
        tmp <- 0.0
        tmp <- tmp + (pi_0 * gaussian_pdf(tickMarks[i], nullmean, (nullvar * rescale_factor)))
        tmp <- tmp / (1.0 - pi_1)
        f0[i] <-  tmp # equalizing the height as it was obtained before rescaling
    }
    
    # f0_factor is updated as we shrink the variances
    i <- 1
    while(tickMarks[i] < nullmean) {
        i <- i + 1
    }
    f0_factor = f[i] / f0[i]
    
    for (i in 1:ntick) {
        f0[i] <- f0[i] * f0_factor
    }
    
    # check if f0_temp falls under the non-parametric estimate of f
    is_under <- 1
    for (i in 1:ntick) {
        if(tickMarks[i] > z_left && tickMarks[i] < z_right) {
            if(f0[i] > (f[i]+tol_diff)) {
                is_under <- 0
            }
        }
    }
    if(is_under==1) {
        nullvar <- nullvar * rescale_factor
        print(paste('Found re-scaling factor: ', rescale_factor, sep = ''))
        break
    }
    rescale_factor <- rescale_factor * 0.99
    ct <- ct + 1
}

if(ct == ntry) {
    print(paste('Unable to rescale variances. Last try: ', rescale_factor, sep = ''))
}
#If it satisfied in the first try, that means we can inflate the variances 
if(ct == 0) {
    print("Inflating variances")
    
    ct <- 0
    rescale_factor <- 1
    while(rescale_factor > 0 && ct < ntry) {
        
        #if(ct > 0 && ct %% 1000 == 0) {
        #    print('error')
        #}
        for(i in 1:ntick) {
            tmp <- 0.0
            tmp <- tmp + (pi_0 * gaussian_pdf(tickMarks[i], nullmean, (nullvar * rescale_factor)))
            tmp <- tmp / (1.0 - pi_1)
            f0[i] <-  tmp # equalizing the height as it was obtained before rescaling
        }
        
        # f0_factor is updated as we shrink the variances
        i <- 1
        while(tickMarks[i] < nullmean) {
            i <- i + 1
        }
        f0_factor = f[i] / f0[i]
        
        for (i in 1:ntick) {
            f0[i] <- f0[i] * f0_factor
        }
        
        # check if f0_temp falls under the non-parametric estimate of f
        is_over <- 0
        for (i in 1:ntick) {
            if(tickMarks[i] > z_left && tickMarks[i] < z_right) {
                if(f0[i] > (f[i]+tol_diff)) {
                    is_over <- 1
                }
            }
        }
        if(is_over==1) {
            nullvar <- nullvar * rescale_factor
            print(paste('Found re-scaling factor: ', rescale_factor, sep = ''))
            break
        }
        rescale_factor <- rescale_factor * 1.05
        ct <- ct + 1
    }
}

for(i in 1:ntick) {
    if(f[i] < f0[i]) {
        f[i] = f0[i]
    }
}

# re-evaluate f0
for(i in 1:ntick) {
    tmp = 0.0
    tmp <- tmp + (pi_0 * gaussian_pdf(tickMarks[i], nullmean, nullvar))
    tmp <- tmp / (1.0 - pi_1)
    f0[i] = tmp
}

#computeProportion()
#get conservative estimates of pi_true: just use the data within a range: (-5,5) for now

#identifying the range to fit null component to the overall density
tmp_r <- NULL

mm <- 0
vv <- 0
mm <- mm + nullmean
vv <- vv + nullvar

tmp <- 1
for(i in 1:ntick) {
    if(abs(tickMarks[i] - mm) < 3.0 * sqrt(vv)) {
        tmp_r = f[i] / f0[i]
        if(tmp_r < tmp) {
            tmp <- tmp_r
        }
    }
}
pi_true = 1 - tmp
print(paste('The estimate of pi(DE) is: ', pi_true, sep = ''))

# Evaluate densities
for(j in 1:np) {
    ## Evaluate f0
    df0[j] <- pi_0 * gaussian_pdf(zstat[j], nullmean, nullvar)
    df0[j] <- df0[j] / (1.0 - pi_1)
    
    ## Evaluate df
    i = 1
    while(tickMarks[i] < zstat[j]  && i <= ntick) {
        i <- i + 1
    }
    df[j] = f[i] + (f[i+1] - f[i]) / (tickMarks[i+1] - tickMarks[i]) * (zstat[j] - tickMarks[i])
}

## Evaluate f1
for(i in 1:ntick) {
    f1[i] <- ( f[i] - (1.0-pi_true) * f0[i] ) / pi_true
    if(f1[i] < 1e-100) {
        f1[i] <- 1e-100
    }
}

#computeFDR
QM_fdr <- NULL
QM_FDRup <- NULL
QM_FDRdown <- NULL

for(j in 1:np) {
    QM_fdr[j] <- (1.0 - pi_true) * df0[j] / df[j]
    if(QM_fdr[j] > 1) {
        QM_fdr[j] <- 1
    } else {
        QM_fdr[j] <- QM_fdr[j]
    }
}
example_qprot_fdr$QM_fdr <- QM_fdr

for(j in 1:np) {
    i <- ntick
    tmpsum <- 0.0
    tmpsum0 <- 0.0
    while(tickMarks[i-1] < zstat[j]) {
        tmpsum <- tmpsum + (0.5 * (f[i] + f[i-1]) * (tickMarks[i] - tickMarks[i-1]))
        tmpsum0 <- tmpsum0 + (0.5 * (f0[i] + f0[i-1]) * (tickMarks[i] - tickMarks[i-1]))
        i <- i - 1
    }
    tmpsum <- tmpsum + (0.5 * (f[i] + df[j]) * (tickMarks[i] - zstat[j]))
    tmpsum0 <- tmpsum0 + (0.5 * (f0[i] + df0[j]) * (tickMarks[i] - zstat[j]))
    
    QM_FDRup[j] = (1.0 - pi_true) * tmpsum0 / tmpsum
    if(QM_FDRup[j] > 1) {
        QM_FDRup[j] <- 1
    } else {
        QM_FDRup[j] <- QM_FDRup[j]
    }
}
example_qprot_fdr$QM_FDRup <- QM_FDRup

for(j in 1:np) {
    i <- 1
    tmpsum <- 0.0
    tmpsum0 <- 0.0
    while(tickMarks[i+1] < zstat[j]) {
        tmpsum <- tmpsum + (0.5 * (f[i+1] + f[i]) * (tickMarks[i+1] - tickMarks[i]))
        tmpsum0 <- tmpsum0 + (0.5 * (f0[i+1] + f0[i]) * (tickMarks[i+1] - tickMarks[i]))
        i <- i + 1
    }
    tmpsum <- tmpsum + (0.5 * (df[j] + f[i]) * (zstat[j]  - tickMarks[i]))
    tmpsum0 <- tmpsum0 + (0.5 * (df0[j] + f0[i]) * (zstat[j]  - tickMarks[i]))
    
    QM_FDRdown[j] = (1.0 - pi_true) * tmpsum0 / tmpsum
    if(QM_FDRdown[j] > 1) {
        QM_FDRdown[j] <- 1
    } else {
        QM_FDRdown[j] <- QM_FDRdown[j]
    }
}
example_qprot_fdr$QM_FDRdown <- QM_FDRdown
setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, '/Qmodel/FDR', sep = ""))
write.csv(example_qprot_fdr, outFile, row.names = FALSE)





