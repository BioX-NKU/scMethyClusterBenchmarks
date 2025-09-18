source("Step_1.r")
library(readr)
#library(dplyr)
files <- Sys.glob("../input/cell_files/*")
print(files)
data <- list()
cnt <- 1
binary_func <- function(x){
    ifelse(x >=  0.5, 1, 0)
}
for (file in files){
    data_file <- read_table(gzfile(file), col_names=c("chr","pos", "Met", "Total"), col_type = 'cddd')
    data_file[,3] <- data_file[,3] / data_file[,4]
    data_file[3] <- apply(data_file[3], 2, binary_func)
    data_frame <- data.frame(data_file[,1], data_file[,2], data_file[,3])
    data[[cnt]] <- data_frame
    cnt <- cnt + 1
}
test <- data[[1]]
print(head(test))
dism_pearson <- Output_DISM(k_cpu = 8, data, method = "Pearson")  #k_cpu means the number of processors to be used
print("done_pearson")
dism_cosine  <- Output_DISM(k_cpu = 8, data, method = "Cosine")
print("done_consine")
dism_hamming <- Output_DISM(k_cpu = 8, data, method = "Hamming")
print("done_hamming")
###Save the resulting distance matrices
write.csv(dism_pearson, file = "./distance/dism_pearson.csv")
write.csv(dism_cosine, file = "./distance/dism_cosine.csv")
write.csv(dism_hamming, file = "./distance/dism_hamming.csv")
