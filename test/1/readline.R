test<-read.table( file="unc.edu.a14d258a-87c9-4c09-81ae-0e3d7d5afca0.1596972.rsem.genes.normalized_results")

con <- file("unc.edu.a14d258a-87c9-4c09-81ae-0e3d7d5afca0.1596972.rsem.genes.normalized_results", "r")
lineCnt = 0
while(1){
  oneline = readLines(con, n = 1)
  if(length(oneline) == 0){
    break
  }
  lineCnt = lineCnt+1
}
close(con)
con <- file("unc.edu.a14d258a-87c9-4c09-81ae-0e3d7d5afca0.1596972.rsem.genes.normalized_results", "r")
line=readLines(con,n=1)
while( length(line) != 0 ) {
  print(line)
  line=readLines(con,n=1)
}
close(con)