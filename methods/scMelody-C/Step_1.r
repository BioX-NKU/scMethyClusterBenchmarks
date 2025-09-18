###Function : Three similarity matrices using genome-wide methylation CpGs were calculated from the input single-cell methylation files###

#'@Aim: Merging chromosome numbers with CpG positions
mat_F <- function(df){
  chr_P <- paste(df[,1], df[,2], sep = "_")
  df <- data.frame(chr_P,df)
  df <- df[,c(1,4)]
  return(df)
  rm(chr_P)
}

#'@Aim: Calculate Pearson similarity matrix
Get_Pearson <- function(j,data){

  cor_cp <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    pv <- cor(dm[,2],dm[,3])
    return(pv)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate Pearson correlation coefficient

  res_cor <- lapply(1:length(data),cor_cp,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}


#'@Aim: Calculate Cosine similarity matrix
Get_Cosine <- function(j,data){

  cor_cos <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    cos <- sum(dm[,2]*dm[,3])/sqrt((sum(dm[,2]^2)*sum(dm[,3]^2)))
    return(cos)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate Cosine correlation coefficient

  res_cor <- lapply(1:length(data),cor_cos,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@Aim: Calculate Hamming similarity matrix
Get_Hamming <- function(j,data){

  cor_ham <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    dm$cot <- ifelse(dm[,2]==dm[,3],1,0)
    dua <- sum(dm$cot)/length(dm$cot)
    return(dua)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate Hamming correlation coefficient

  res_cor <- lapply(1:length(data),cor_ham,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@Aim: Output similarity matrix using parallel computing


library(parallel)
Output_SM <- function(k_cpu,data,method){
  if(method=='Cosine'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Cosine,data)
    stopCluster(cl)
    res_mat <- do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Hamming'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Dual,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Pearson'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Pearson,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  else{
    print('Something wrong with your settings,please check...')
  }
}