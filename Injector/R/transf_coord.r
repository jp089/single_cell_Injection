#' Transfering coordinates
#'
#' A function that transfer the Trajectory coordinates
#' (comp_1, and comp_2, or and comp_3) within the metadata
#' of 'reference' Seurat object. The function use Seurat package
#' workflow of mapping and annotating query datasets.
#'
#' the 'reference' obtained from load_ref() and
#' the 'query' obtained from load_query(). This function includes
#' preprocessing, normlization, search of variables features and anchors,
#' for data integration of Seurat package workflow.
#'
#' @param reference A Seurat object obtained from load_ref function
#' @param dims number of PCA dimensions used to calculate the reference's UMAP
#' @param n_comps number of components of your trajectory (2 or 3)
#' @param integrated1 A Seurat object obtained from intg1 function
#' @param data_name name of the data (same as used in load_query)
#' @return A Seurat object, with only query data and transfered
#' trajectory coordinates
#'
#' @export
transf_coord <- function(reference, dims, n_comps, integrated1, data_name) {
    if (length(data_name) == 1) {
        reference <- Seurat::RunUMAP(object = reference, dims = 1:dims,
                     verbose = FALSE, return.model = TRUE)
        trf_anchors <- Seurat::FindTransferAnchors(reference = reference,
                       query = integrated1, normalization.method = "LogNormalize",
                       reference.reduction = "pca", dims = 1:50)
        if (n_comps == 2) {
            integrated2 <- Seurat::TransferData(anchorset = trf_anchors,
                           reference = reference, query = integrated1,
                           refdata = list(comp_1 = "comp_1", comp_2 = "comp_2"),
                           store.weights = FALSE)
            } else if (n_comps == 3) {
                integrated2 <- Seurat::TransferData(anchorset = trf_anchors,
                               reference = reference, query = integrated1,
                               refdata = list(comp_1 = "comp_1", comp_2 = "comp_2",
                               comp_3 = "comp_3"), store.weights = FALSE)
                       }
        query_data_name <- paste("q", data_name, sep = "_")
        integrated2 <- integrated2[,grep(query_data_name, colnames(integrated2),
                       value = TRUE)]
    }
}
