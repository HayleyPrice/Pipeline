# remove.packages("rstan")
# if (file.exists(".RData")) file.remove(".RData")

Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
install.packages("V8", repos = 'http://cloud.r-project.org/')

stopifnot(nchar(Sys.getenv("R_LIBS_USER")) > 0)
dir.create(Sys.getenv('R_LIBS_USER'),
           showWarnings = FALSE,
           recursive = TRUE)
install.packages(c('StanHeaders',
                    'Rcpp',
                    'RcppEigen',
                    'inline',
                    'BH',
                   'rstan'),
                 repos = 'http://cloud.r-project.org/',
                 lib = Sys.getenv("R_LIBS_USER"),
                 dependencies=TRUE)

example(stan_model, package = "rstan", run.dontrun = TRUE)

#install.packages("Rcpp", repos = 'http://cloud.r-project.org/')