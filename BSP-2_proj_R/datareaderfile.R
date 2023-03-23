# rm(list = ls())

library(psych)

dataRapl <- read.csv(file = 'dataRAPL.csv')
dataMeter <- read.csv(file = 'dataMETER.csv')


row_sums_rapl <- rowSums(dataRapl[, c(9:11)])
subset_data_rapl <- dataRapl[, c(9:11)]
subset_data_rapl$Total <- row_sums_rapl
subset_data_Meter <- dataMeter[,c(7,8)]


summary(subset_data_rapl["Total"])
summary(subset_data_Meter["Energy.Meter"])

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

describe_orig <- function(subset1,subset2) {
  
benchmark_data <- data.frame(energy_rapl = numeric(),energy_meter = numeric())
  
  # Run the benchmark 30 times and collect data
  for (i in 1:30) {
    # print(i)
    # energy_rapl <- sum(subset1[((33*(i-1))+1):(33*(i))]) # Example value, replace with actual data
    # print(subset1[((33*(i-1))+1):(33*(i))])
    # energy_meter <- sum(subset2[((33*(i-1))+1):(33*(i))]) # Example value, replace with actual data
    # 
    energy_rapl <- subset1[((33*(i-1))+1)] # Example value, replace with actual data
    # print(subset1[((33*(i-1))+1)])
    energy_meter <- subset2[((33*(i-1))+1)]# Example value, replace with actual data
    # 
    # Add data to the benchmark data frame
    benchmark_data[i,] <- c(energy_rapl, energy_meter)
  }
  
  describe(benchmark_data)
  
}

print(workers[30])

describe_orig2 <- function(subset1,subset2) {
  benchmark_data <- data.frame(threads = numeric(), energy_rapl = numeric(),energy_meter = numeric())
  for (i in 1:30) {
    energy_rapl <- subset1[((33*(i-1))+1)] # Example value, replace with actual data
    energy_meter <- subset2[((33*(i-1))+1)]# Example value, replace with actual data
    benchmark_data[i,] <- c(workers[1],energy_rapl, energy_meter)
  }
  describe(benchmark_data)
}

describe_orig3 <- function(subset1,subset2) {
  benchmark_data <- data.frame(threads = numeric(), energy_rapl = numeric(),energy_meter = numeric())
  for (i in 1:30) {
    energy_rapl <- subset1[((33*(i-1))+2)] # Example value, replace with actual data
    energy_meter <- subset2[((33*(i-1))+2)]# Example value, replace with actual data
    benchmark_data[i,] <- c(workers[2],energy_rapl, energy_meter)
  }
  describe(benchmark_data)
}

# ...
describe_orig2(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"])
describe_orig3(subset_data_rapl[0:990,"Total"],subset_data_Meter[0:990,"Energy.Meter"])


# Get unique number of threads
workers <- unique(dataMeter[0:990,]$n.workers)
size_chunk<- 990
  
  
# Perform Shapiro-Wilk normality test for each number of threads
for( j in 1:5){

list_ps <- list()

  for (i in 1:length(workers)) {
  # Subset data for current number of threads
  # print(workers[i])
  new_dat <- dataMeter[((size_chunk*(j-1))+1):(size_chunk*(j)),]
  data <- new_dat[new_dat$n.workers == workers[i],]
   # print(dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),][dataMeter[(size_chunk*(j-1)):(size_chunk*(j)),]$n.workers == workers[i],]$n.workers)
  # Perform Shapiro-Wilk normality test on energy consumption (RAPL)
  shapiro_test <- shapiro.test(data$Energy.Meter)
  list_ps[i]<- shapiro_test$p.value
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

# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("gridExtra")
library(ggplot2)
library(dplyr)
library(gridExtra)


# RAPL summaries
summary(filter(subset_data_rapl, n.s.workers == 0)$Total)

summary(filter(subset_data_rapl, n.s.workers == 0, n.workers == 1)$Total)
summary(filter(subset_data_rapl, n.s.workers == 0, n.workers == 2)$Total)

summary(filter(subset_data_rapl, n.s.workers == 1)$Total)
summary(filter(subset_data_rapl, n.s.workers == 1, n.workers == 1)$Total)
summary(filter(subset_data_rapl, n.s.workers == 1, n.workers == 2)$Total)


summary(subset_data_rapl$Total)

# METER summaries
summary(filter(dataMeter, n.secondary.workers == 0)$Energy.Meter)

summary(filter(dataMeter, n.secondary.workers == 0, n.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 0, n.workers == 2)$Energy.Meter)

summary(filter(dataMeter, n.secondary.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 1, n.workers == 1)$Energy.Meter)
summary(filter(dataMeter, n.secondary.workers == 1, n.workers == 2)$Energy.Meter)


summary(dataMeter$Energy.Meter)


# Correlation check
  
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

print(get_list_with_avar_workers_Rapl(3))


cor(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)),method = c("pearson", "kendall", "spearman"))

cor(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)))
cor(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)))
cor(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)))
cor(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)))

# General
cor(subset_data_Meter$Energy.Meter,subset_data_rapl$Total) #same

cor.test(unlist(get_list_with_avar_workers_Meter(1)),unlist(get_list_with_avar_workers_Rapl(1)))

cor.test(unlist(get_list_with_avar_workers_Meter(2)),unlist(get_list_with_avar_workers_Rapl(2)))
cor.test(unlist(get_list_with_avar_workers_Meter(3)),unlist(get_list_with_avar_workers_Rapl(3)))
cor.test(unlist(get_list_with_avar_workers_Meter(4)),unlist(get_list_with_avar_workers_Rapl(4)))
cor.test(unlist(get_list_with_avar_workers_Meter(5)),unlist(get_list_with_avar_workers_Rapl(5)))

cor.test(subset_data_Meter$Energy.Meter,subset_data_rapl$Total) #same

library("ggpubr")

give_me_plot <- function(chunk){

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

give_me_plot(5)

#RandM$Meter <- subset_data_Meter[,2]


# General
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

get_qq_plot_MR(3)


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
