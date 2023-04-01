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