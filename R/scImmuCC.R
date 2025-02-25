#' @title scImmuCC
#' @description Creating Hierarchical_annotation for scRNA-Seq data immune cell.
#' @details Input takes a cells-genes matrix with cell unique barcodes as column names and gene names as row names and returns the cells annotation.
#' @param count  a matrix with cell unique barcodes as column names and gene names as row names .
#' @param Non_Immune Whether non-immune cells are included in the matrix.
#' @return Data frames with the barcodes and cell types, and some maps.
#' @import GSVA
#' @importFrom GSVA gsva
#' @import Matrix
#' @export
#' @examples test_data
#' library(GSVA)
#' library(Seurat)

scImmuCC <- function(count, Non_Immune=TRUE, resultsdir="", genelist_dir="data", min_cells=50, ncores=1){

  library(BiocParallel)
  register(MulticoreParam(workers = ncores))

  dir.create(resultsdir, showWarnings = FALSE)
  genelists <- load_genelists(genelist_dir)
  results <- list()

  run_scImmuCC <- function(count, genelist, name, resultsdir="", min_cells=50){
    if (ncol(count) < min_cells) return(NULL)
    if (is.null(genelist)) return(NULL)
    ssGSEA_result <- scImmuCC_main(count, genelist, file.path(resultsdir, name))
    seurat_result <- seurat_Heatmap(count, genelist, ssGSEA_result, file.path(resultsdir, name))
    return(ssGSEA_result)
  }

  # layer 0/1
  if (Non_Immune==TRUE){
    layer0_result <- run_scImmuCC(count, genelists$layer0_genelist, "Layer0", resultsdir=resultsdir, min_cells=min_cells)
    results$Layer0 = layer0_result
    immune <- layer0_result[which(layer0_result[,2]=="Immune"),]
    ssGSEA_result <- run_scImmuCC(count[,immune[,1]], genelists$layer1_genelist, "Layer1", resultsdir=resultsdir, min_cells=min_cells)
    results$Layer1 = ssGSEA_result
  } else {
    ssGSEA_result <- run_scImmuCC(count, genelists$layer1_genelist, "Layer1", resultsdir=resultsdir, min_cells=min_cells)
    results$Layer1 = ssGSEA_result
  }

  cell_type <- unique(ssGSEA_result[,2])

  # layer 2
  layer2_celltypes = intersect(cell_type, c("Tcell", "Bcell", "DC", "NK", "Macrophage", "ILC"))
  results$Layer2 = bpapply(setNames(layer2_celltypes, layer2_celltypes), function(celltype){
    ix = ssGSEA_result[which(ssGSEA_result[,2]==celltype),][,1]
    run_scImmuCC(count[,ix], genelists[[paste0(celltype, "_genelist")]], paste0("Layer2_", celltype), resultsdir=resultsdir, min_cells=min_cells)
  })

  # layer 3
  if ("Tcell" %in% layer2_celltypes){
    genelists$CD4_T_genelist <- genelists$CD4_genelist
    genelists$CD8_T_genelist <- genelists$CD8_genelist
    layer3_celltypes = c("CD4_T", "CD8_T")
    results$Layer3 = bpapply(setNames(layer3_celltypes, layer3_celltypes), function(celltype){
      ix = ssGSEA_result[which(ssGSEA_result[,2]==celltype),][,1]
      run_scImmuCC(count[,ix], genelists[[paste0(celltype, "_genelist")]], paste0("Layer3_", celltype), resultsdir=resultsdir, min_cells=min_cells)
    })

  }

  return(results)
}


#' Load gene lists
#'
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
load_genelists <- function(dir="data"){
  files = list.files(dir, pattern = "_genelist.rda", full.names = TRUE)
  geneenv <- new.env()
  lapply(files, function(file){
    load(file, envir = geneenv)
  })
  varnames <- ls(geneenv)
  lapply(setNames(varnames, varnames), function(var){
    get(var, envir = geneenv)
  })
}










