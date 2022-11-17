#' Integration - 1
#'
#' A function that integrate both Seurat objects reference and query.
#' The function use the 'reference' obtained from load_ref() and
#' the 'query' obtained from load_query(). This function includes
#' preprocessing, normlization, search of variables features and anchors,
#' for data integration of Seurat package workflow.
#'
#' @param reference A Seurat object obtained from load_ref function
#' @param query A Seurat object obtained from load_query function
#' @return A Seurat object, where the reference and query are integrated
#' @export
intg1 <- function(reference,query) {
    #Listing reference and query
    list1 <- list(reference,query)
    names(list1) <- c("ref", "qry")
    #Normalizing query
    list1[[2]] <- Seurat::NormalizeData(list1[[2]], verbose = FALSE)
    #Variable features for both
    for (i in 1:length( list1)) {
        list1[[i]] <- Seurat::FindVariableFeatures(list1[[i]],
                      selection.method = "vst",
                      nfeatures = 3000, verbose = FALSE)
        }
    #Finding anchors for integration
    intg1_anchors <- Seurat::FindIntegrationAnchors(object.list = list1,
                    anchor.features = 3000)
    #Integration
    integrated1 <- Seurat::IntegrateData(anchorset = intg1_anchors)
    #Scaling
    integrated1 <- Seurat::ScaleData(integrated1, verbose = FALSE)
}