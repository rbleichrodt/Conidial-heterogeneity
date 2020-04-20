    library(plyr)
    setwd("D:/Data/FACS/ISE single reps x3 cbind/Lectins in time/")
  
  folder = getwd()
  dirs = list.dirs(path = folder, full.names = TRUE, recursive = TRUE)
  len1=0
  len2=0
  len3=0
  dat1=0
  dat2=0
  dat3=0
  
  #define no of replicates
  reps = 2
  
  for (j in 1:1){
    setwd(dirs[j])
    folder = getwd()
  
  
  file_list <- list.files(path = folder, pattern = "txt", all.files = FALSE,
                          full.names = FALSE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  print(file_list)
  
  for(i in 1:(length(file_list)/reps)){
    for(j in 1:reps){
      if(j == 1){
      dat1 <- read.table(file_list[(i*reps)+j-reps], header = FALSE)
      }
      if(j == 2){
        dat2 <- read.table(file_list[(i*reps)+j-reps], header = FALSE)
      }
      if(j == 3){
        dat3 <- read.table(file_list[(i*reps)+j-reps], header = FALSE)
      }
    }
    dat1_v <- as.vector(dat1[['V1']])
    dat2_v <- as.vector(dat2[['V1']])
    dat3_v <- as.vector(dat3[['V1']])
    
    n <- max(length(dat1_v), length(dat2_v), length(dat3_v))
    length(dat1_v) <- n                      
    length(dat2_v) <- n
    length(dat3_v) <- n
    data <- cbind(dat1_v,dat2_v,dat3_v)
      
    name_end <- nchar(file_list[i*reps])
      name <- paste0(substr(file_list[i*reps], 1, name_end-reps),"_cbind.txt")
      print(name)
      write.table(data, file = name, row.names = FALSE, col.names = FALSE)
  }
  }