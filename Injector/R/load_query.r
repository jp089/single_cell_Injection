#' Loading query
#'
#' This function loads gene expression data (query) by
#' using Read10X (Seurat package), and prepeared to injection
#' on the reference. If ADT data is availible, it could be added.
#' If it exists a list of interest, this could be also take
#' into account to subset both gene expression and ADT data.
#' Finally, a 'data name' is required. Cells in both gene
#' expression and ADT data need to be in the same order.
#'
#' @param query Path to the directorey containing matrix.mtx,
#' genes.tsv, and barcodes.tsv files provided by 10X.
#' @param adt Path to the input file (tab separated format)
#' @param list_of_interest Path to the input file (tab separated format)
#' @param data_name name of the data query
#' @return A Seurat object, with adt data as an 'assay'
#' @export
load_query <- function(query, adt, list_of_interest, data_name) {
    query_data <- Seurat::Read10X(query)
    colnames(query_data) <- paste(data_name, colnames(x = query_data),
                                  sep = "_")
    if (base::file.exists(adt) == TRUE) {
        query_adt <- read.table(adt, header = TRUE, row.names = 1, sep = "\t")
        colnames(query_adt) <- paste(data_name, colnames(x = query_adt),
                                     sep = "_")
        if(base::file.exists(list_of_interest) == TRUE) {
            listInt <- read.table(list_of_interest)$V1
            listInt <- paste(data_name, listInt, sep="_")
            query_data_name <- paste("q", data_name, sep = "_")
            if (length(colnames(query_data) == length(colnames(query_adt)))) {
                query_data <- query_data[, colnames(query_data) %in% listInt]
                colnames(query_data) <- gsub(data_name,query_data_name,
                                             colnames(query_data))
                query_adt <- query_adt[, colnames(query_adt) %in% listInt]
                colnames(query_adt) <- gsub(data_name, query_data_name,
                                            colnames(query_adt))
                } else {
                        print("query and adt don't have 
                               the same number of cells")
                        }
            }
        query <- SeuratObject::CreateSeuratObject(counts = query_data,
                                                      min.cells = 1)
        query[["adt"]] <- SeuratObject::CreateAssayObject(
                          counts = query_adt)
        } else {
                if(base::file.exists(list_of_interest) == TRUE) {
                    listInt <- read.table(list_of_interest)$V1
                    listInt <- paste(data_name, listInt, sep="_")
                    query_data_name <- paste("q", data_name, sep = "_")
                    query_data <- query_data[, colnames(query_data) %in% listInt]
                    colnames(query_data) <- gsub(data_name,query_data_name,
                                                 colnames(query_data))
                } 
        query <- SeuratObject::CreateSeuratObject(counts = query_data,
                                                  min.cells = 1)
                }
    }