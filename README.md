# yipCat
### Description
YipCat package can realize some of the more commonly used analysis requirements of transcription, such as **trajectory analysis,
heatmap, cell interaction, imputeweight etc**. , you can view the function description inside the package.
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
