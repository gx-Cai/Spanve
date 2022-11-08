#!/usr/bin/env Rscript
library(MERINGUE)
library(Matrix)
library(reticulate)
library(argparser)

meringue_wrapper <- function(adata_file_path,workdir) {
    log.file = paste0(workdir,'//','run.log')
    out.file = paste0(workdir,'//','meringue.out.csv')
    
    sc <- import("scanpy")
    adata <- sc$read_h5ad(adata_file_path)
    counts = as.matrix(adata$X)
    coords = as.data.frame(adata$obsm$get('spatial'))

    colnames(counts) = adata$var_names$values
    row.names(counts) = adata$obs_names$values
    colnames(coords) = c('x','y')
    row.names(coords) = as.array(adata$obs_names$values)

    counts = t(counts)
    counts <- Matrix::Matrix(counts, sparse = TRUE)

    # CPM normalize
    mat <- MERINGUE::normalizeCounts(counts = counts, 
        log=TRUE,
        verbose=TRUE
    )
    
    st = Sys.time()
    
    w <- getSpatialNeighbors(coords, filterDist = 150)
    I <- getSpatialPatterns(mat, w)
    results.filter <- filterSpatialPatterns(
        mat = mat,
        I = I,
        w = w,
        adjustPv = TRUE,
        alpha = 1,
        minPercentCells = 0.10,
        verbose = TRUE,
        details = TRUE
    )
    
    et = Sys.time()
    spend = as.numeric(as.POSIXct(Sys.time()) - as.POSIXct(st),units="secs")
    write(x=paste('time cost:',spend),file = log.file,append = T)
    write.csv(x = results.filter,file = out.file)
}


argv = arg_parser("run_meringue")
argv = add_argument(argv,"--data",help = 'adata_file_path')
argv = add_argument(argv,"--workdir",help = 'workdir')
argv = parse_args(argv)

meringue_wrapper(argv$data,argv$workdir)