# rPeaks

This is a re-write of the R Peaks package maintained by M.Kondrin
based upon a C library written by Miroslav Morhac. rPeaks was
written to use current build techniques from RStudio to build.

The documentation is automatically generated from tags included in the
source files and processed by roxygen2. I did this because of issues
with R-3.2.0 using the Peaks package from CRAN.

Note: the demos from the Peaks package did not work with R-3.2.0. Hadley
Wickham notes that demos are not checked with R CMD check, but
vignettes are. The Peaks demos were converted to vignettes here.
