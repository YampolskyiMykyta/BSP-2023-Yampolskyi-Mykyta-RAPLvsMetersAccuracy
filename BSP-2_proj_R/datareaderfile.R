install.packages("ggplot2")
install.packages("dplyr")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("psych")
library(ggplot2)
library(dplyr)
library(gridExtra)
library("ggpubr")
library(psych)

dataRapl <- read.csv(file = 'dataRAPL.csv')
dataMeter <- read.csv(file = 'dataMETER.csv')


row_sums_rapl <- rowSums(dataRapl[, c(9:11)])
subset_data_rapl <- dataRapl[, c(9:11)]
subset_data_rapl$Total <- row_sums_rapl
subset_data_rapl$n.workers = dataRapl$n.workers
subset_data_rapl$n.secondary.workers = dataRapl$n.secondary.workers

subset_data_Meter <- dataMeter[,c(7,8)]


summary(subset_data_rapl["Total"])
summary(subset_data_Meter["Energy.Meter"])

# Get unique number of threads
workers <- unique(dataMeter[0:990,]$n.workers)
size_chunk<- 990

# Methods to get averages for each chunk in Meter and Rapl
get_list_with_avar_workers_Meter <- function(chunk) {
  list_ps <- list()
  for (i in 1:length(workers)) {
    
    new_dat <- dataMeter[((size_chunk*(chunk-1))+1):(size_chunk*(chunk)),]
    data <- new_dat[new_dat$n.workers == workers[i],]
    list_ps[i] <- mean(data$Energy.Meter)
  }
  return(list_ps)
}

get_list_with_avar_workers_Rapl <- function(chunk) {
  list_ps <- list()
  for (i in 1:length(workers)) {
    
    new_dat <- subset_data_rapl[((size_chunk*(chunk-1))+1):(size_chunk*(chunk)),]
    data <- new_dat[new_dat$n.workers == workers[i],]
    list_ps[i] <- mean(data$Total)
  }
  return(list_ps)
}


read_chunks <- function(subset) {
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
    #chunks[[counter]] <- chunk
    counter <- counter+1
  }
}

read_chunks(subset_data_rapl["Total"])
read_chunks(subset_data_Meter["Energy.Meter"])


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



# NOT AVERAGES
# Perform Shapiro-Wilk normality test for each number of threads
for( j in 1:5){

list_ps <- list()
  # cat("-----", j)
  for (i in 1:length(workers)) {
  # Subset data for current number of threads
  # print(workers[i])
  new_dat <- dataMeter[((size_chunk*(j-1))+1):(size_chunk*(j)),]
  data <- new_dat[new_dat$n.workers == workers[i],]
   # print(dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),][dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),]$n.workers == workers[i],]$n.workers)
  # Perform Shapiro-Wilk normality test on energy consumption (RAPL)
  shapiro_test <- shapiro.test(data$Energy.Meter)
  list_ps[i]<- shapiro_test$p.value
  # cat("--", i)
  # print(list_ps[i])
  # Print results
  # cat("Number of threads:", workers[i], "\n")
  # cat("W =", shapiro_test$statistic, ", p-value =", shapiro_test$p.value, "\n\n")
}
print(mean(unlist(list_ps)))

}

# Perform Shapiro-Wilk normality test for each number of threads
for( j in 1:5){
  
  list_ps <- list()
  
  for (i in 1:length(workers)) {
    # Subset data for current number of threads
    # print(workers[i])
    new_dat <- subset_data_rapl[((size_chunk*(j-1))+1):(size_chunk*(j)),]
    data <- new_dat[new_dat$n.workers == workers[i],]
    # print(dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),][dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),]$n.workers == workers[i],]$n.workers)
    # Perform Shapiro-Wilk normality test on energy consumption (RAPL)
    shapiro_test <- shapiro.test(data$Total)
    list_ps[i]<- shapiro_test$p.value
    # Print results
    # cat("Number of threads:", workers[i], "\n")
    # cat("W =", shapiro_test$statistic, ", p-value =", shapiro_test$p.value, "\n\n")
  }
  print(mean(unlist(list_ps)))
}

# AVERAGES
# METER
for(i in 1:5){
  shapiro_test_iso <- shapiro.test(unlist(get_list_with_avar_workers_Meter(i)))
  print(shapiro_test_iso$p.value)
}
# RAPL
for(i in 1:5){
  shapiro_test_iso <- shapiro.test(unlist(get_list_with_avar_workers_Rapl(i)))
  print(shapiro_test_iso$p.value)
}

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


# Correlation tests

print(get_list_with_avar_workers_Rapl(3))


# CORRELATION for different chunks, method: PEARSON

cor(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)),method = c("pearson"))

cor(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)),method = c("pearson"))
cor(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)),method = c("pearson"))
cor(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)),method = c("pearson"))
cor(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)),method = c("pearson"))

# General
cor(subset_data_Meter$Energy.Meter,subset_data_rapl$Total,method = c("pearson")) #same

cor.test(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)),method = c("pearson"))

cor.test(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)),method = c("pearson"))
cor.test(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)),method = c("pearson"))
cor.test(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)),method = c("pearson"))
cor.test(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)),method = c("pearson"))

cor.test(subset_data_Meter$Energy.Meter,subset_data_rapl$Total) #same


# CORRELATION for different chunks, method: SPEARMAN

cor(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)),method = c("spearman"))

cor(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)),method = c("spearman"))
cor(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)),method = c("spearman"))
cor(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)),method = c("spearman"))
cor(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)),method = c("spearman"))

# General
cor(subset_data_Meter$Energy.Meter,subset_data_rapl$Total) #same

cor.test(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)),method = c("spearman"))

cor.test(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)),method = c("spearman"))
cor.test(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)),method = c("spearman"))
cor.test(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)),method = c("spearman"))
cor.test(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)),method = c("spearman"))

cor.test(subset_data_Meter$Energy.Meter,subset_data_rapl$Total) #same


# PEARSON plots

give_me_plot_pear <- function(chunk){

my_data <- dataMeter[(1+(990*(chunk-1))):(33+(990*(chunk-1))),c(1:2)]
print(my_data$version[1])

my_data$Meter <- unlist(get_list_with_avar_workers_Meter(chunk))
my_data$Rapl <- unlist(get_list_with_avar_workers_Rapl(chunk))
# plotting 

ggscatter(my_data, x = "Meter", y = "Rapl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Joules M", ylab = "Joules R")
}
give_me_plot_pear(1)
give_me_plot_pear(2)
give_me_plot_pear(3)
give_me_plot_pear(4)
give_me_plot_pear(5)


# SPEARMAN plots

give_me_plot_spear <- function(chunk){
  
  my_data <- dataMeter[(1+(990*(chunk-1))):(33+(990*(chunk-1))),c(1:2)]
  print(my_data$version[1])
  
  my_data$Meter <- unlist(get_list_with_avar_workers_Meter(chunk))
  my_data$Rapl <- unlist(get_list_with_avar_workers_Rapl(chunk))
  # plotting 
  
  ggscatter(my_data, x = "Meter", y = "Rapl", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Joules M", ylab = "Joules R")
}
give_me_plot_spear(1)
give_me_plot_spear(2)
give_me_plot_spear(3)
give_me_plot_spear(4)
give_me_plot_spear(5)


# General
RandM <- subset_data_rapl[,4:6]
RandM$Meter <- subset_data_Meter$Energy.Meter

ggscatter(RandM, x = "Meter", y = "Total", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Joules M", ylab = "Joules R")



get_qq_plot_MR <- function(chunk){

qqnorm( unlist(get_list_with_avar_workers_Meter(chunk)))
qqline(unlist(get_list_with_avar_workers_Meter(chunk)))

qqnorm( unlist(get_list_with_avar_workers_Rapl(chunk)))
qqline( unlist(get_list_with_avar_workers_Rapl(chunk)))

}
get_qq_plot_MR(1)
get_qq_plot_MR(2)
get_qq_plot_MR(3)
get_qq_plot_MR(4)
get_qq_plot_MR(5)

# General
qqnorm(RandM$Meter)
qqline(RandM$Meter)

qqnorm(RandM$Total)
qqline(RandM$Total)


# Hyp_testing
test <- function(data1, data2) {
  wilcox.test(data1, data2)
}

test(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)))

test(RandM$Meter,RandM$Total) # --


# ------------------------------ METER -----------------------------

# WORKING WITH ORIGINAL:

# Perform Shapiro-Wilk normality test for number of threads equals to 1, 2, 8, 32, 64
# with q-q plots
for(j in list(1,2,5,17,33)){
  list_ps <- list()
  for (i in 1:30) {
    # Subset data for current number of threads
    new_dat <- dataMeter[j+ 33*(i-1),]
    list_ps[i] <- new_dat$Energy.Meter
  }
  shapiro_test <- shapiro.test(unlist(list_ps))
  print(shapiro_test$p.value )
  
  qqnorm( unlist(list_ps))
  qqline(unlist(list_ps))
}

# WORKING WITH TXACTs:

# Perform Shapiro-Wilk normality test for number of threads equals to 1, 2, 8, 32, 64
# with q-q plots
for(p in 1:4){
for(j in list(1,2,5,17,33)){
  list_ps <- list()
  for (i in 1:30) {
    # Subset data for current number of threads
    new_dat <- dataMeter[size_chunk*p+ j+ 33*(i-1),]
    list_ps[i] <- new_dat$Energy.Meter
  }
  shapiro_test <- shapiro.test(unlist(list_ps))
  print(shapiro_test$p.value )
  
  qqnorm( unlist(list_ps))
  qqline(unlist(list_ps))
}
  print("-----------")
}


# ------------------------------ RAPL -----------------------------

# WORKING WITH ORIGINAL:

# Perform Shapiro-Wilk normality test for number of threads equals to 1, 2, 8, 32, 64
# with q-q plots
for(j in list(1,2,5,17,33)){
  list_ps <- list()
  for (i in 1:30) {
    # Subset data for current number of threads
    new_dat <- subset_data_rapl[j+ 33*(i-1),]
    list_ps[i] <- new_dat$Total
  }
  shapiro_test <- shapiro.test(unlist(list_ps))
  print(shapiro_test$p.value )
}

# WORKING WITH TXACT8:

# Perform Shapiro-Wilk normality test for number of threads equals to 1, 2, 8, 32, 64
# with q-q plots
for(p in 1:4){
  for(j in list(1,2,5,17,33)){
    list_ps <- list()
    for (i in 1:30) {
      # Subset data for current number of threads
      new_dat <- subset_data_rapl[size_chunk*p+ j+ 33*(i-1),]
      list_ps[i] <- new_dat$Total
    }
    shapiro_test <- shapiro.test(unlist(list_ps))
    print(shapiro_test$p.value )
  }
  print("-----------")
}


# ------------------------------- METER vs RAPL correlation ----------------------

# WORKING WITH ORIGINAL:
# check correlation with Pearson and spearman for threads 1,2,8,32,64

get_pear_t_spear_f_plot <- function(j,who){
  list_dm <- list()
  list_dr <- list()
  for (i in 1:30) {
    list_dm[i] <- dataMeter[j+ 33*(i-1),"Energy.Meter"]
    list_dr[i] <- subset_data_rapl[j+ 33*(i-1),"Total"]
  }
  print(cor(unlist(list_dm),unlist(list_dr),method = c("pearson")))
  print(cor(unlist(list_dm),unlist(list_dr),method = c("spearman")))
  print("--")
  print(test(unlist(list_dm),unlist(list_dr)))
  
  my_data <- dataMeter[1:30,c(1:2)]
  
  my_data$Meter <- unlist(list_dm)
  my_data$Rapl <- unlist(list_dr)
  
 
  # plotting 
  # pearson
  if(who == TRUE){
  ggscatter(my_data, x = "Meter", y = "Rapl", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Joules M p", ylab = "Joules R p")
  }
  else{
  # spearman
  ggscatter(my_data, x = "Meter", y = "Rapl",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Joules M s", ylab = "Joules R s")
  }
}

# threads 1,2,8,32,64 => index: 1,2,5,17,33
# pearson
get_pear_t_spear_f_plot(1,TRUE)
get_pear_t_spear_f_plot(2,TRUE)
get_pear_t_spear_f_plot(5,TRUE)
get_pear_t_spear_f_plot(17,TRUE)
get_pear_t_spear_f_plot(33,TRUE)

# spearman
get_pear_t_spear_f_plot(1,FALSE)
get_pear_t_spear_f_plot(2,FALSE)
get_pear_t_spear_f_plot(5,FALSE)
get_pear_t_spear_f_plot(17,FALSE)
get_pear_t_spear_f_plot(33,FALSE)
