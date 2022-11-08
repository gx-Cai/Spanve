#!/usr/bin/env Rscript

library(Giotto)
library(reticulate)
library(argparser)

giotto_binspect_rank <- function(adata_file_path,workdir,method) {
    py.info = reticulate::py_config()
    sc <- import("scanpy")
    adata <- sc$read_h5ad(adata_file_path)
    counts = as.matrix(adata$X)
    coords = as.data.frame(adata$obsm$get('spatial'))

    colnames(counts) = adata$var_names$values
    row.names(counts) = adata$obs_names$values

    colnames(coords) = c('x','y')
    coords$cell_ID = as.array(adata$obs_names$values)
    
    instructions <- createGiottoInstructions(
        show_plot = F,
        return_plot = F,
        save_plot = F,
        save_dir = workdir,
        python_path = py.info$python
    )

    data <- createGiottoObject(raw_exprs = t(counts), spatial_locs = coords, instructions = instructions )

    ## normalize
    data <- normalizeGiotto(gobject = data, scalefactor = 6000, verbose = T)

    if (method=='rank'){
        st = Sys.time()
        
        data <- createSpatialNetwork(gobject = data, method = 'kNN', k = 5, maximum_distance_knn = 400, name = 'spatial_network')
        spatialgenes <- binSpect(data,bin_method='rank', calc_hub = T, hub_min_int = 5, spatial_network_name = 'spatial_network', do_parallel = TRUE, cores=32)
        
        et = Sys.time()
        spend = as.numeric(as.POSIXct(Sys.time()) - as.POSIXct(st),units="secs")
        write(x=paste('time cost:',spend),file = paste0(workdir,'//','run.log'),append = T)
        write.csv(x=spatialgenes,file=paste0(workdir,'//','gitto.rank.csv'))
    }
    else{
        st = Sys.time()
        spatialgenes = silhouetteRank(gobject = data, expression_values = 'normalized')
        et = Sys.time()
        spend = as.numeric(as.POSIXct(Sys.time()) - as.POSIXct(st),units="secs")
        write(x=paste('time cost:',spend),file = paste0(workdir,'//','run.log'),append = T)
        write.csv(x=spatialgenes,file=paste0(workdir,'//','gitto.sirank.csv'))
    }
}


argv = arg_parser("run_gitto")
argv = add_argument(argv,"--data",help = 'adata_file_path')
argv = add_argument(argv,"--workdir",help = 'workdir')
argv = add_argument(argv,"--method",help = 'rank or sirank')
argv = parse_args(argv)

giotto_binspect_rank(argv$data,argv$workdir,argv$method)
