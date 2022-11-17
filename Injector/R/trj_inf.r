#' Trajectory injection
#'
#' trajectory model obtain by infer_trajectory (dynwrap package)
#'
#' @param integrated2 A Seurat object obtained from transf_coord function
#' @param n_comps number of components of your trajectory (2 or 3)
#' @param model_rooted Path to the trajectory model (.Robj)
#' @param root name of the milestone used for rooting the reference trajetory
#'
#' @export
trj_inf <- function(integrated2, n_comps, model_rooted, root) {
    counts_q <- Matrix::t(integrated2@assays$integrated@data)
    gene_expression_q <- Matrix::t(integrated2@assays$integrated@scale.data)
    if ((table(colnames(counts_q) %in% colnames(gene_expression_q)))[[1]] ==
       dim(counts_q)[2]) {
        print("Gene names ok")
        # colnames(counts_q) <- gsub("\\-", "\\.", colnames(counts_q))
        if (n_comps == 2) {
            mst_q <- as.matrix(integrated2@meta.data[, c("predicted.comp1",
                               "predicted.comp_2")])
            } else if (n_comps == 3) {
            mst_q <- as.matrix(integrated2@meta.data[, c("predicted.comp1",
                               "predicted.comp_2", "predicted.comp_3")])
                   }   
        if (dim(counts_q)[2] == dim(gene_expression_q)[2]) {
            dataset_q <- dynwrap::wrap_expression(expression = gene_expression_q,
                                                  counts = counts_q)
            } else {
                counts_q <- counts_q[, (colnames(counts_q) %in%
                                      colnames(gene_expression_q))]
                dataset_q <- dynwrap::wrap_expression(
                             expression = gene_expression_q, counts = counts_q)
            }
        load(model_rooted)
        expression <- dataset_q$expression
        if (n_comps == 2) {
            space <- as.data.frame(mst_q)
            space$predicted.comp1 <- as.numeric(space$predicted.comp1)
            space$predicted.comp_2 <- as.numeric(space$predicted.comp_2)
            } else if (n_comps == 3) {
                space <- as.data.frame(mst_q)
            space$predicted.comp1 <- as.numeric(space$predicted.comp1)
            space$predicted.comp_2 <- as.numeric(space$predicted.comp_2)
            space$predicted.comp_3 <- as.numeric(space$predicted.comp_3)
            }
        centers <- model_rooted$dimred_milestones
        mt_ids <- model_rooted$milestone_ids
        mt_ntw <- model_rooted$milestone_network
        output <- dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
                  dynwrap::add_dimred_projection(
                           milestone_ids = mt_ids,
                           milestone_network = mt_ntw,
                           dimred = as.matrix(space),
                           dimred_milestones = centers)
            
        if (length(root) == 1) {
            library(dplyr)
            output_rooted1 <- output %>% dynwrap::add_root(root_milestone_id = root)
            } else {
                output_rooted1 <- output
                }
        output_rooted1 <- dynwrap::add_pseudotime(output_rooted1, pseudotime = NULL)
        } else {
            print("verify gene names")
        }
    }