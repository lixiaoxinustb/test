test<-read.table( file="test.txt",header=T)
row1 <- c()
row2 <- "normalized_count"
split<-function(sample){
    for(i in 1:nrow(sample)){ 
      m <- as.character(sample[i, 1])
      l <- strsplit(m, '')
      j = length(l[[1]])
      n <- l[[1]][1]
      for(k in 2:j){
        
        if(l[[1]][k]=="|"){
          break
        }
        n <- c(n,l[[1]][k])
      }
      h <- paste(n[1:k-1],collapse="")
      
      row1 <-c(row1,h)
      
      
    }
  print(row1)
  }
  
  split(test)
  

