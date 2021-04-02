# yipCat
### Description
YipCat package can realize some of the more commonly used analysis requirements of transcription, such as 
** trajectory analysis,heatmap, cell interaction,imputeweight,cellphonedb,agescore etc ** , you can view the function description inside the package.
Mainly the implementation style of the diagram is more advanced to look at, in particular, the **mass spectrometry streaming data** is also integrated and 
can achieve the interaction of **cellranger cloupe files**, very convenient.

### Installtation
you can install this package via command:
```r
install.packages("yipCat_1.0.1.tar.gz")
```

### requirement
```r
require(Seurat)
require(dplyr)
require(ggplot2)
require(S4Vectors)
require(nabor)
...
```

### trajectory 
#### The implementation of spline trajectory, the specific method description mainly draws on part of the ArchR method, the diagram is as follows：
you can caculate trajecory by follwing commands,for example
```r
seurat<-addSeuratTrajectory(object=seurat,trajectory=c("memory B","Naive B","Plasma"),groupBy="label_fine",embedding="pca") # caculate trajecory
se<-getSeuratTrajectory(seurat)  # get trajectory object
ArchR::plotTrajectoryHeatmap(se)  # plot trajectory heatmap
```
![trajectory](inst/extdata/testTrajectoyHeatmap_page-0001.jpg)

### cacaulte ImputeWeight
#### computes imputations weights that describe each cell as a linear combination of many cells based on a MAGIC diffusion matrix.
for example 
```r
seurat<-ImputeWeights(seurat,reducedDims="pca",nRep=2,sampleCells=5000)  # caculate impute weight
weight<-getImputeWeight(seurat)  # return impute weight
EmbPlot(seurat,colorBy="matrix",features="S100A9",embedding="umap",imputeWeights=NULL)  # first plot
EmbPlot(seurat,colorBy="matrix",features="S100A9",embedding="umap",imputeWeights=weight)  # second plot
```
first plot  ![notImpute](inst/extdata/notImpute-S100A9_page-0001.jpg)  
second plot ![Impute](inst/extdata/impute-S100A9_page-0001.jpg) 
