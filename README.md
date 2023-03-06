# single_cell_Injection

# Presentation

#Containt the package "Injector", which allows to 'inject' or add new single-cells to reference trajectory. It will use a previously defined Seurat object (using Seurat) with the reference cells, and the information of the infered trajectory (using Dynverse), and the new query cells as a Seurat object.

# How to install?

## Dependencies 

First, you need to install "devtools" package from CRAN

```{r}
install.packages("devtools")
```

Normally by installing the package all dependcies are isntalled with it, however some libraries will need to be installed independently.

```{r}
devtools::install_github("dynverse/dyno")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
```

## Installation of Injector 

```{r}
devtools::install_github("jp089/single_cell_Injection/Injector")
library(Injector)
```

# How to use?

A vignette ("single_cell_Injection/Vignette.Rmd") is included within the GitHub repository to present how to use this library, as well as some neccessary files for this test. 
