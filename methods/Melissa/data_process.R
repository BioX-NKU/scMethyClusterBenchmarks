library(readr)
files <- Sys.glob("../input/cell_files/*")
print(files)
dt_obj <- list()
dt_obj$met <- list()
dt_obj$anno_region <- NULL
dt_obj$opts <- list()
dt_obj$opts$basis_obj <- NULL
dt_obj$opts$N <- NULL
dt_obj$opts$M <- NULL
dt_obj$opts$K <- NULL
dt_obj$opts$pi_k <- NULL
dt_obj$opts$max_cov <- NULL
dt_obj$opts$cluster_var <- NULL
n_cells <- length(files) 
dt_obj$opts$cell_names <- paste("cell_", seq(1, n_cells, by = 1), sep = "")
cnt <- 1
binary_func <- function(x){
    ifelse(x >=  0.5, 1, 0)
}
for (file in files){
    data_file <- read_table(gzfile(file), col_names=c("chr","pos","Met","Total"), col_type = 'cddd')
    data_file[,3] <- data_file[,3] / data_file[,4]
    data_file[3] <- apply(data_file[3], 2, binary_func)
    data_frame <- data.frame(data_file[,1], data_file[,2], data_file[,3])
    data_frame <- data_frame[order(data_frame$chr),]
    dt_obj$met[[cnt]] <- list()
    list_matrix <- split(data_frame, data_frame$chr)
    for(chr_group in list_matrix){
        chr_group <- chr_group[,-1]
        chr_group <- chr_group[order(chr_group$pos),]
        matrix_group <- as.matrix(chr_group)
	dt_obj$met[[cnt]][[length(dt_obj$met[[cnt]]) + 1]] <- matrix_group
    }
    cnt <- cnt + 1
}
names(dt_obj$met) <- dt_obj$opts$cell_names
