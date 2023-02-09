
library(BayesSpace)
for (data_id in c('151507','151508','151509','151510','151669','151670','151671','151672','151673','151674','151675','151676')){

    dlpfc <- getRDS("2020_maynard_prefrontal-cortex", data_id)
    saveRDS(dlpfc,paste0('./',data_id,'.rds'))
    sceasy::convertFormat(dlpfc, from="sce", to="anndata",outFile=paste0(data_id,'.h5ad'))
}