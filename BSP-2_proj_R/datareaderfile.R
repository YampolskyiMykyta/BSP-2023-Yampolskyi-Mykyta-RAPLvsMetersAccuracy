# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("gridExtra")
# install.packages("ggpubr")
# install.packages("psych")
library(dplyr)
library(ggplot2)
library(gridExtra)
library("ggpubr")
library(psych)

dataRapl <- read.csv(file = 'dataRAPL.csv')
dataMeter <- read.csv(file = 'dataMETER.csv')


subset_data_rapl <- dataRapl[, c(9:11)]
subset_data_rapl$Total <- rowSums(dataRapl[, c(9:11)])
subset_data_rapl$n.workers = dataRapl$n.workers
subset_data_rapl$n.secondary.workers = dataRapl$n.secondary.workers

subset_data_Meter <- dataMeter[,c(7,8)]


summary(subset_data_rapl["Total"])
summary(subset_data_Meter["Energy.Meter"])

# Get unique number of threads
workers <- unique(dataMeter[0:990,]$n.workers)
size_chunk<- 990

#' get_list_with_avar_workers_Meter
#'
#' Retrieves a list of average workers for the Meter data.
#'
#' @param chunk The chunk number to process.
#'
#' @return A list containing the average workers for each chunk.
#'
get_list_with_avar_workers_Meter <- function(chunk) {
# This function takes a parameter "chunk" and returns a list of average energy consumption values
#   from the Meter data for different numbers of workers. It iterates over each worker and selects
#   the corresponding subset of data based on the chunk parameter. It calculates the mean of the "Energy.Meter"
#   column for each workers subset and stores it in the list_ps variable. Finally, it returns the list_ps containing
#   the average energy consumption values for each worker

  list_ps <- list()
  for (i in 1:length(workers)) {
    
    new_dat <- dataMeter[((size_chunk*(chunk-1))+1):(size_chunk*(chunk)),]
    data <- new_dat[new_dat$n.workers == workers[i],]
    list_ps[i] <- mean(data$Energy.Meter)
  }
  return(list_ps)
}

#' get_list_with_avar_workers_Rapl
#'
#' Retrieves a list of average workers for the Rapl data.
#'
#' @param chunk The chunk number to process.
#'
#' @return A list containing the average workers for each chunk.
#'
get_list_with_avar_workers_Rapl <- function(chunk) {
   #Similar to the previous function, this function also takes a parameter "chunk"
   #and returns a list of average energy consumption values. However, in this case,
   #it works with the Rapl data. It follows the same logic as the previous function
   #but calculates the mean of the "Total" column instead of "Energy.Meter".
  list_ps <- list()
  for (i in 1:length(workers)) {
    
    new_dat <- subset_data_rapl[((size_chunk*(chunk-1))+1):(size_chunk*(chunk)),]
    data <- new_dat[new_dat$n.workers == workers[i],]
    list_ps[i] <- mean(data$Total)
  }
  return(list_ps)
}

#' read_chunks
#'
#' Reads and prints summary statistics for data chunks.
#'
#' @param subset The subset of data to process.
#'
read_chunks <- function(subset) {
   #It defines a list of values "txacts" representing different transactions.
   #The function initializes a counter variable and iterates over the chunks.
   #If the counter is 1, it selects the original chunk from the subset and
   #prints a summary of the chunk. Otherwise, it selects the chunk based on
   #the chunk size and prints a summary with the corresponding transaction value.
  chunk_size <- 990
  txacts <- list(0,1,2,8,64)
  counter <- 1
  for (i in seq(chunk_size, chunk_size*5,by = chunk_size)) {
    if (counter == 1) {
      chunk <- subset[(i-chunk_size):i,]
      print("----- ORIGINAL -----")
      print(summary(chunk))
    } else {
      chunk <- subset[(i-chunk_size):i,]
      print(paste("----- TXACT-",txacts[[counter]]," -----"))
      print(summary(chunk))
    }
    counter <- counter+1
  }
}

read_chunks(subset_data_rapl["Total"])
read_chunks(subset_data_Meter["Energy.Meter"])

#' describe_orig_chunk
#'
#' Describes the original chunk using benchmark data.
#'
#' @param subset1 The first subset of data.
#' @param subset2 The second subset of data.
#' @param chunk The chunk number.
#'
describe_orig_chunk <- function(subset1,subset2,chunk) {
  benchmark_data <- data.frame(threads = numeric(), energy_rapl = numeric(),energy_meter = numeric())
  for (i in 1:30) {
    energy_rapl <- subset1[((33*(i-1))+chunk)] # Example value, replace with actual data
    energy_meter <- subset2[((33*(i-1))+chunk)]# Example value, replace with actual data
    benchmark_data[i,] <- c(workers[chunk],energy_rapl, energy_meter)
  }
  describe(benchmark_data)
}

describe_orig_chunk(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"],1)
describe_orig_chunk(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"],2)
describe_orig_chunk(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"],11)
describe_orig_chunk(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"],12)
describe_orig_chunk(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"],33)


# RAPL summaries
summary(filter(subset_data_rapl, n.secondary.workers == 0)$Total)

summary(filter(subset_data_rapl, n.secondary.workers == 0, n.workers == 1)$Total)
summary(filter(subset_data_rapl, n.secondary.workers == 0, n.workers == 2)$Total)

summary(filter(subset_data_rapl, n.secondary.workers == 1)$Total)
summary(filter(subset_data_rapl, n.secondary.workers == 1, n.workers == 1)$Total)
summary(filter(subset_data_rapl, n.secondary.workers == 1, n.workers == 2)$Total)

summary(subset_data_rapl$Total)

# METER summaries
summary(filter(dataMeter, n.secondary.workers == 0)$Energy.Meter)

summary(filter(dataMeter, n.secondary.workers == 0, n.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 0, n.workers == 2)$Energy.Meter)

summary(filter(dataMeter, n.secondary.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 1, n.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 1, n.workers == 2)$Energy.Meter)

summary(dataMeter$Energy.Meter)


# ------------------------------ METER -----------------------------
# ------------ Normality testing -----------

#' get_some_nsw_data Function
#'
#' This function subsets data based on specified conditions and performs a Shapiro-Wilk test on the subsetted data.
#' It then prints the result of the test and generates a QQ plot if the p-value is equal to 0.05.
#'
#' @param nsw Number of secondary workers
#' @param nw Number of workers
#' @param who Variable to pass data from Meter (Conservative or less) or RAPL
#' @param what Variable(s) to be subsetted from the data frame
#' @return None
#' @export
get_some_nsw_data <- function(nsw,nw,who,what){
    # Subset data for current number of threads
    new_dat <- filter(who, n.secondary.workers == nsw,n.workers==nw)[what]

    shapiro_test <- shapiro.test(unlist(new_dat))
    if(shapiro_test$statistic > 0.1 ){
    if(shapiro_test$p.value > 0.05 ){cat(nw, " p.val = ", shapiro_test$p.value,"       +","\n" )}
    else{cat(nw, " p.val = ", shapiro_test$p.value,"   -","\n" )}
    # for more information:
    # if(shapiro_test$p.value > 0.05 ){cat(nw, " W = ", shapiro_test$statistic," p.val = ", shapiro_test$p.value,"       +","\n" )}
    # else{cat(nw, " W = ",  shapiro_test$statistic," p.val = ", shapiro_test$p.value,"   -","\n" )}
        
    if(shapiro_test$p.value == 0.05){
    qqnorm( unlist(new_dat))
    qqline(unlist(new_dat))
    }}
    else{
      cat(nw,"\n")
    }
}

# WORKING WITH ORIGINAL:

# Perform Shapiro-Wilk normality test for all number of threads
# with q-q plots

for(p in workers){
  get_some_nsw_data(0,p,dataMeter_ver2,"Energy.Meter")
}

# WORKING WITH TXACTs:

# Perform Shapiro-Wilk normality test for all number of threads
# with q-q plots

for(k in list( 1, 2, 8, 64)){
  cat("\n\n---SEC.WORK ",k," ---\n\n")
  for(p in workers){
    get_some_nsw_data(k,p,dataMeter,"Energy.Meter")
    #get_some_nsw_data(k,p,dataMeter_ver2,"Energy.Meter")
  }
}


# ------------------------------ RAPL -----------------------------
# ------------ Normality testing -----------

# WORKING WITH ORIGINAL:

# Perform Shapiro-Wilk normality test for all number of threads equals
# with q-q plots

for(p in workers){
  get_some_nsw_data(0,p,subset_data_rapl,"Total")
}

# WORKING WITH TXACTs:

# Perform Shapiro-Wilk normality test for all number of threads equals
# with q-q plots.

for(k in list(1, 2, 8, 64)){
  cat("\n\n---SEC.WORK ",k," ---\n\n")
  for(p in workers){
    get_some_nsw_data(k,p,subset_data_rapl,"Total")
  }
}


# ------------------------------- METER vs RAPL correlation ----------------------

## NOT AVERAGES
# WORKING WITH ORIGINAL:

# check correlation with Pearson and spearman for threads 1,2,8,32,64

#' get_pear_t_spear_f_plot_new Function
#'
#' This function calculates the Pearson and Spearman correlation coefficients between the Meter and Rapl
#' energy consumption data, and plots a scatter plot with a regression line.
#'
#' @param nsw Numeric value representing the number of secondary workers.
#' @param who Logical value indicating whether to calculate the Spearman correlation coefficient (TRUE) or the Pearson correlation coefficient (FALSE).
#' @param Meterdataset Data frame containing the Meter energy consumption data.
#'
#' @return The function does not return any value. It prints the correlation coefficients and generates a scatter plot.
#'
get_pear_t_spear_f_plot_new <- function(nsw,who,Meterdataset){
  cat("Transactional models = ", nsw, "\n")
  list_dr <- list()
  list_dm <- list()
  
  for(i in workers){
    new_dat <- filter(Meterdataset, n.secondary.workers == nsw,n.workers==i)$Energy.Meter
    list_dm[i] <- median(unlist(new_dat))
    new_dat2 <- filter(subset_data_rapl, n.secondary.workers == nsw,n.workers==i)$Total
    list_dr[i] <- median(unlist(new_dat2))
  }
  
  p <-  cor(unlist(list_dm),unlist(list_dr),method = c("pearson"))
  sp <-  cor(unlist(list_dm),unlist(list_dr),method = c("spearman"))

    cat(p, "   p  + \n")
    cat(sp, "   sp  + \n")

  my_data <- Meterdataset[1:33,c(1:2)]

  my_data$Meter <- unlist(list_dm)
  my_data$Rapl <- unlist(list_dr)
  
  if(who == TRUE){
    ggscatter(my_data, x = "Meter", y = "Rapl",
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "spearman",
              xlab = "Joules M p", ylab = "Joules R p")
  }
   if(who == FALSE){
    ggscatter(my_data, x = "Meter", y = "Rapl",
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Joules M p", ylab = "Joules R p")
  }
}


get_spear_corr_plot(TRUE,dataMeter_ver2)
#get_spear_corr_plot(TRUE,dataMeter)

# threads 1,2,8,32,64
# spearman - TRUE, pearson - FALSE
for(k in list(0, 1, 2, 8, 64)){
  cat("\n\n---SEC.WORK ",k," ---\n\n")
    get_pear_t_spear_f_plot_new(k,TRUE,dataMeter)
    #get_pear_t_spear_f_plot_new(k,TRUE,dataMeter_ver2)
    #get_some_nsw_data(k,p,subset_data_rapl,"Total")
}

# spearman
get_pear_t_spear_f_plot_new(0,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(1,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(2,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(8,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(64,TRUE,dataMeter_ver2)

# pearson
get_pear_t_spear_f_plot_new(0,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(1,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(2,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(8,TRUE,dataMeter_ver2)
get_pear_t_spear_f_plot_new(64,TRUE,dataMeter_ver2)


# If you want to get correlation tests according to the normality:
## NOT AVERAGES
get_pear_t_spear_f_plot_accord_to_norm <- function(nsw,nw){
  list_dr <- list()
  list_dm <- list()
  
  list_dm <- filter(dataMeter, n.secondary.workers == nsw,n.workers==nw)$Energy.Meter
  list_dr <- filter(subset_data_rapl, n.secondary.workers == nsw,n.workers==nw)$Total

  shapiro_test_Meter <- shapiro.test(unlist(list_dm))
  shapiro_test_Rapl <- shapiro.test(unlist(list_dr))
  
  if(shapiro_test_Meter$statistic > 0.7 && shapiro_test_Rapl$statistic >0.7){
    if(shapiro_test_Meter$p.value > 0.05 && shapiro_test_Rapl$p.value > 0.05){

     cat(nw, " W = ", shapiro_test_Meter$statistic," p.val = ", shapiro_test$p.value,"       +","\n" )
      
      cat(nw," ",cor(unlist(list_dm),unlist(list_dr),method = c("pearson")),"\n")
      
      }
    else if(shapiro_test_Meter$p.value < 0.05 && shapiro_test_Rapl$p.value < 0.05){
      
      cat(nw, " W = ",  shapiro_test$statistic," p.val = ", shapiro_test$p.value,"   -","\n" )
      
      cat(nw," ",cor(unlist(list_dm),unlist(list_dr),method = c("spearman")),"\n")
    }
    else{
      cat(nw,"\n")
    }
  }
  
  # my_data <- dataMeter_ver2[1:30,c(1:2)]
  # 
  # my_data$Meter <- unlist(list_dm)
  # my_data$Rapl <- unlist(list_dr)
  # 
  # # plotting 
  # # pearson
  # if(who == TRUE){
  #   ggscatter(my_data, x = "Meter", y = "Rapl", 
  #             add = "reg.line", conf.int = TRUE, 
  #             cor.coef = TRUE, cor.method = "pearson",
  #             xlab = "Joules M p", ylab = "Joules R p")
  # }
  # else{
  #   # spearman
  #   ggscatter(my_data, x = "Meter", y = "Rapl",
  #             add = "reg.line", conf.int = TRUE,
  #             cor.coef = TRUE, cor.method = "spearman",
  #             xlab = "Joules M s", ylab = "Joules R s")
  # }
}
for( i in workers){
get_pear_t_spear_f_plot_accord_to_norm(2,i)
}

## AVERAGES
NEW_get_list_with_avar_workers_Meter <- function(nsw,fromwho) {
  list_ps <- list()
  for (i in workers) {
    new_dat <- filter(fromwho, n.secondary.workers == nsw,n.workers==i)$Energy.Meter
    list_ps[i] <- mean(new_dat)
  }
  return(list_ps)
}


NEW_get_list_with_avar_workers_Rapl <- function(nsw) {
  list_ps <- list()
  for (i in workers) {
    new_dat <- filter(subset_data_rapl, n.secondary.workers == nsw,n.workers==i)$Total
    list_ps[i] <- mean(new_dat)
  }
  return(list_ps)
}

for(k in list(0, 1, 2, 8, 64)){
  cat("ver1 --- THREADS: ",k," ---\n")
  cat("pearson: ",cor(unlist(NEW_get_list_with_avar_workers_Meter(k,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(k)),method = c("pearson")),"\n"  )
  cat("spearman: ",cor(unlist(NEW_get_list_with_avar_workers_Meter(k,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(k)),method = c("spearman")),"\n\n")

}

for(k in list(0, 1, 2, 8, 64)){
  cat("vers2 --- THREADS: ",k," ---\n")
  cat("pearson: ",cor(unlist(NEW_get_list_with_avar_workers_Meter(k,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(k)),method = c("pearson")),"\n"  )
  cat("spearman: ",cor(unlist(NEW_get_list_with_avar_workers_Meter(k,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(k)),method = c("spearman")),"\n\n")

}

# ----------------------------------- SIGNIFICANCE CHECK --------------------------------

# AVERAGES
# works faster
wilcox.test(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)))
wilcox.test(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)))
wilcox.test(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)))
wilcox.test(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)))
wilcox.test(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)))

# same but slower
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(0,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(0)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(1,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(1)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(2,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(2)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(8,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(8)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(64,dataMeter)),unlist(NEW_get_list_with_avar_workers_Rapl(64)))

# version 2
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(0,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(0)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(1,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(1)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(2,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(2)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(8,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(8)))
wilcox.test(unlist(NEW_get_list_with_avar_workers_Meter(64,dataMeter_ver2)),unlist(NEW_get_list_with_avar_workers_Rapl(64)))



# NOT AVERAGES
get_some_nsw_data_Wilcox <- function(nsw,nw,who,what){
  # Subset data for current number of threads
  new_dat <- filter(who, n.secondary.workers == nsw,n.workers==nw)[what]
  new_dat2 <- filter(subset_data_rapl, n.secondary.workers == nsw,n.workers==nw)$Total
  
  wilcox_test <- wilcox.test(unlist(new_dat),unlist(new_dat2), exact=FALSE)
  if(wilcox_test$p.value > 0.05){cat(nw, wilcox_test$p.value,"\n" )}
  else{cat(nw, wilcox_test$p.value,"  -","\n" )}
  
}

for( p in workers){
get_some_nsw_data_Wilcox(0,p,dataMeter,"Energy.Meter")

}

for(k in list(1, 2, 8, 64)){
  for(p in workers){
    get_some_nsw_data_Wilcox(k,p,dataMeter,"Energy.Meter")
  }
  print("-----------")
}

for(p in workers){
  get_some_nsw_data_Wilcox(0,p,dataMeter_ver2,"Energy.Meter")
}

for(k in list(1, 2, 8, 64)){
  for(p in workers){
    get_some_nsw_data_Wilcox(k,p,dataMeter_ver2,"Energy.Meter")
  }
  print("-----------")
}

# wilcox.test(list_M,list_R)
wilcox.test(RandM$Meter,RandM$Total) 


# ----------------------------------- PLOTS SAVING --------------------------------

save_into_file <- function(subfold,namefile,num,func){

  if (!dir.exists("plots_folder")) {
    dir.create("plots_folder")
  }
  if (!dir.exists("plots_folder/Correlation")) {
    dir.create("plots_folder/Correlation")
  }
  if (!dir.exists("plots_folder/QQ-plots")) {
    dir.create("plots_folder/QQ-plots")
  }
  # saving in jpeg format
  # jpeg("plots_folder/saving_jpeg4.jpeg")
  # jpeg(file = paste0("plots_folder/",subfold,"/",namefile, num, ".png"))
  
  ggsave(file = paste0("plots_folder/", subfold, "/", namefile, num, ".jpeg"), plot =  func(num)  ,device = "jpeg",dpi = 100)
  
  # dev.off()
}

save_into_file("Correlation","corr_pear_orig_avar",1,get_pear_t_spear_f_plot_new)
