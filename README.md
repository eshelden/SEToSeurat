# SEToSeurat
R function that accepts a StandardizedExperiment object and returns a Seurat object
use:
msd<-SEToSeurat(combinedResults, sName = "myProject")
Note that Monocle will also import a Seurat object produced by the above using the importCDS() function.
