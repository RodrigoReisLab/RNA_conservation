if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2",dependencies = TRUE)
}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {
  install.packages("gridExtra",dependencies = TRUE)
}
if("qpcR" %in% rownames(installed.packages()) == FALSE) {
  install.packages("qpcR",dependencies = TRUE)
}

library(ggplot2)
library(gridExtra)
library(qpcR)

dirPath <- getwd()
setwd(dirPath)

files <- Sys.glob("*_oneLine_nuclInfo.csv")
#files <- head(all_files, n=2)

first_file <- as.data.frame(read.csv(file = files[1], header = FALSE, sep = "\t"))
baseN <- paste(tools::file_path_sans_ext(files[1]))
baseName <- sub(baseN, pattern = "_oneLine_nuclInfo$", replacement = "")

colnames(first_file) <- c("chrID", "genomeSize", "Ncount", "percentNs")
first_file$genome_name <- baseName
df <- cbind(first_file[,5], stack(first_file[2:3]))
df[,2] <- df[,2]/10^6

for (f in files) {
  if(f != files[1])
  {
    baseN <- paste(tools::file_path_sans_ext(f))
    baseName <- sub(baseN, pattern = "_oneLine_nuclInfo$", replacement = "")

    rf <- as.data.frame(read.csv(file = f, header = FALSE, sep = "\t"))
    colnames(rf) <- c("chrID", "genomeSize", "Ncount", "percentNs")
    rf$genome_name <- baseName
    stacks <- cbind(rf[,5], stack(rf[2:3]))
    stacks[,2] <- stacks[,2]/10^6 
    df <- qpcR:::rbind.na(df, stacks)
  }
  
  rf = data.frame()
  stacks = data.frame()
}
colnames(df) <- c("genomeName","count_size","index")
df[df == 0] <- NA

#Plot data
p <- ggplot(na.omit(df), aes(x = genomeName, y=count_size, fill=index, color=index)) + theme_classic() + 
  geom_boxplot() + ylab("Size (in million)") + 
  scale_y_continuous(limits=c(0,100), breaks = seq(0, 100, 10)) +	
  theme(axis.text.x = element_text(size = 12, angle=45, colour="black", hjust=1), 
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 12, colour="black"),
        legend.position="top")

###############################################################################
#Plot percentage of Ns in the genomes
df_per <- data.frame()
df_per <- as.data.frame(first_file[,4:5])

for (f in files) {
  if(f != files[1])
  {
    baseN <- paste(tools::file_path_sans_ext(f))
    baseName <- sub(baseN, pattern = "_oneLine_nuclInfo$", replacement = "")
    
    rf <- as.data.frame(read.csv(file = f, header = FALSE, sep = "\t"))
    colnames(rf) <- c("chrID", "genomeSize", "Ncount", "percentNs")
    rf$genome_name <- baseName
    stacks_per <- as.data.frame(rf[,4:5])
    df_per <- as.data.frame(qpcR:::rbind.na(df_per, stacks_per))
  }
  
  rf = data.frame()
  stacks_per = data.frame()
}
colnames(df_per) <- c("percentNs","genomeName")
df_per[df_per == 0] <- NA

s <- ggplot(na.omit(df_per), aes(x = genomeName, y=percentNs)) + 
  theme_classic() + geom_boxplot(fill="#f69806", color="#f69806") + 
  scale_y_continuous(limits=c(0,100), breaks = seq(0, 100, 10)) + 
  ylab("Number of positions as 'Ns' (%)") +
  xlab("Plant Genomes") +
  theme(axis.text.x = element_text(size = 12, angle=45, colour="black", hjust=1), 
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 12, colour="black"))

#Plot data
pdf(file = "genomeSize_Ncontent_distribution.pdf", paper="a4r", height=20, width = 100)
grid.arrange(p, s, ncol=1)
dev.off()
