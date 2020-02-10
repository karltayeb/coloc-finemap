library(coloc)
library(reticulate)

paths = read.csv('/work-zfs/abattle4/karl/gp_fine_mapping/simulation/single_snp_sim/sim_paths.txt')

run_coloc <- function(path){
    p = as.character(path)
    d1 = py_load_object(paste(p, list.files(p, pattern = '*_data1'), sep = ''))
    d2 = py_load_object(paste(p, list.files(p, pattern = '*_data2'), sep = ''))
    results = coloc.abf(d1, d2)
    write.csv(x = results[2], file = paste(p, 'coloc_abf', sep = ''))
    write.csv(x = results[1], file = paste(p, 'coloc_abf_summary', sep = ''))
}

for (path in paths$X0){
    print(as.character(path))
    run_coloc(path)
}
