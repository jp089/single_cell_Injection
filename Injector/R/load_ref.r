#' Loading reference
#'
#' This function loads a reference Seurat object with trajectory
#' inference coordinates as metadata. It assumes that the 'Seurat'
#' object was previously generate by using CreateSeuratObject (Seurat
#' package). Trajectory coordinates (comp_1 and comp_2) were obtained
#' with one dynverse package methods and then added to the 'Seurat' object,
#' by using AddMetaData (Seurat package). If the reference
#' have < 20000 cells, it will use all of them. If it is between 20000
#' and 25000, it will randomly select half of the data. And if it have > 25000
#' cells, it will randomly select one third of the data.
#'
#' @param ref Path to the input file (.h5seurat)
#' @return A Seurat object, with trajectory coordinates as metadata
#' @export
load_ref <- function(ref) {
    reference <- SeuratDisk::LoadH5Seurat(ref)
    if ((length(colnames(reference))) <= 20000) {
        reference <- reference
        } else if ((length(colnames(reference))) > 20000 &&
                   (length(colnames(reference))) <= 25000) {
            random <- base::sample(1:(length(colnames(reference))),
                                   round(((length(colnames(reference))) / 2)))
            reference <- reference[, random]
            } else if ((length(colnames(reference))) > 25000) {
                random <- base::sample(1:(length(colnames(reference))),
                                       round(((length(colnames(reference))) /3))
                                       )
                reference <- reference[, random]
                } else {
                    reference <- reference
                    }
    }