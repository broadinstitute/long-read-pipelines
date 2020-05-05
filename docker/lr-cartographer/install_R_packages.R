options(warn = 2)     # treat warnings as errors, otherwise script can fail silently if a package fails to install

InstallPackageFromArchive = function(packageName, packageURL) {
		# make sure to use http not https as this will give an "unsupported URL scheme" error
		if (!(packageName %in% rownames(installed.packages()))) {
				install.packages(packageURL, repos = NULL, type = "source", clean = TRUE)
		}
}

dependencies = c("qqman", 
								 "digest", "gtable", "MASS", "plyr", "reshape2", "scales", "tibble", "lazyeval", # for ggplot2
								 "ggplot2")
#repos <- c("http://cran.mtu.edu")
repos <- c("https://CRAN.R-project.org")
install.packages(dependencies, repos = repos, clean = TRUE)

install.packages("XML", repos = "http://www.omegahat.net/R")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DRIMSeq")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("stageR")


#InstallPackageFromArchive("getopt", "http://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.0.tar.gz")
#InstallPackageFromArchive("optparse", "http://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.3.2.tar.gz")
#InstallPackageFromArchive("data.table", "http://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.10.4-2.tar.gz")
#InstallPackageFromArchive("gsalib", "http://cran.r-project.org/src/contrib/gsalib_2.1.tar.gz")
#InstallPackageFromArchive("ggplot2", "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_2.2.1.tar.gz")
#InstallPackageFromArchive("dplyr", "http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.8.0.1.tar.gz")

q(save = "no")
