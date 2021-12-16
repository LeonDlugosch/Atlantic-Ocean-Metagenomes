ExtractField = function(x, n, sep = "_"){
  return(unlist(lapply(strsplit(x, sep), "[", n)))
}

library(stringr)
path = Sys.getenv("outDir")
path = paste(path, "/06_TabularFiles/faa/", sep = "")
print(path)
setwd(path)
file.names = dir(pattern = ".faa")
minCov = Sys.getenv("minCov")
print(paste("Minimal coverage:", minCov), quote = F)

for (i in 1:length(file.names)) {
  Sample = substr(file.names[i], start = 1, stop = nchar(file.names[i])-4)
  print(paste("Reading ", Sample, ".faa ...", sep = ""), quote = F)
  
  fasta.tab = read.delim(file.names[i],
                         header = F,
                         stringsAsFactors = F)
  print("Renaming Sequences...", quote = F)
  names(fasta.tab) = c("ID", "Sequence")
  fasta.tab$ID = str_replace(fasta.tab$ID, " ", "_")
  split = str_split(fasta.tab[, 1], "_")
  ContigNr = unlist(lapply(split, "[[", 2))
  GeneNR = unlist(lapply(split, "[[", 7))
  
  fasta.tab$Coverage = as.numeric(unlist(lapply(split, "[[", 6)))
  fasta.tab$ID_short = paste(Sample, ContigNr, GeneNR, sep = "_")
  fasta.tab$complete = grepl(fasta.tab$ID, pattern = "partial=00")
  print("Running coverage filter...", quote = F)
  fasta.tab.filtered = fasta.tab[which(fasta.tab$Coverage >= minCov), ]
  fasta.tab.filtered = fasta.tab.filtered[, -3]
  fasta.tab.filtered$Seq = fasta.tab.filtered$Sequence
  fasta.tab.filtered = fasta.tab.filtered[, -c(1:2)]
  fasta.tab.filtered.c = fasta.tab.filtered[which(fasta.tab.filtered$complete == T),]
  fasta.tab.filtered.c = fasta.tab.filtered.c[,-2] 
  fasta.tab.filtered = fasta.tab.filtered[,-2]
  file = paste(Sample, "_FR.faa", sep = "")
  file.c = paste(Sample, "_FR_complete.faa", sep = "")
  
  print(paste("Writing", file, "...", sep = " "), quote = F)
  
  write.table(
    file = file,
    fasta.tab.filtered,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  write.table(
    file = file.c,
    fasta.tab.filtered.c,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  ) 
  print("============================================================",quote = F)
}