#!/usr/bin/env Rscript

library(SPARK)
library(Matrix)
library(reticulate)
library(argparser)

sparkx_wrapper <- function(adata_file_path,workdir) {
    log.file = paste0(workdir,'//','run.log')
    out.file = paste0(workdir,'//','sparkx.out.csv')
    
    sc <- import("scanpy")
    adata <- sc$read_h5ad(adata_file_path)
    counts = as.matrix(adata$X)
    coords = as.data.frame(adata$obsm$get('spatial'))
    colnames(counts) = adata$var_names$values
    row.names(counts) = adata$obs_names$values
    colnames(coords) = c('x','y')
    # coords$cell_ID = as.array(adata$obs_names$values)

    sink(file = log.file, append = TRUE)
    st = Sys.time()
    sparkX <- sparkx(t(counts), coords, numCores=56, option="mixture")
    sink()
    et = Sys.time()
    spend = as.numeric(as.POSIXct(Sys.time()) - as.POSIXct(st),units="secs")
    write(x=paste('time cost:',spend),file = log.file,append = T)
    write.csv(x = sparkX$res_mtest,file = out.file)
}


argv = arg_parser("run_sparkx")
argv = add_argument(argv,"--data",help = 'adata_file_path')
argv = add_argument(argv,"--workdir",help = 'workdir')
argv = parse_args(argv)

sparkx_wrapper(argv$data,argv$workdir)