counts <- read.table("urothelium_counts.txt",header=T)
barplot(counts, main="Fig 1D 1E",
  xlab="Number of SNV", col=c("darkblue","blue"),
  legend = rownames(counts))