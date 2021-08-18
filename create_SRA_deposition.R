library(dplyr)
library(openxlsx)
library(parallel)
d <- read.table('AAVengeR/configs/Sabatino.samples.config', sep = ',', header =  TRUE)
d <- select(d, subject, sample)
a <- read.table('data/sampleDetails.tsv', sep = '\t', header = TRUE)

d <- left_join(d, a, by = 'sample')



cluster <- makeCluster(30)
d$n <- 1:nrow(d)

# Create mock quality scores because AAVenger does not save this information.
# Raw reads and AAVenger software archived at Zenodo.
invisible(parLapply(cluster, split(d, d$n), function(x){
  library(ShortRead)
  setwd('/home/everett/canine_hemophilia_AAV')
  sampleReads <- list.files('/home/everett/canine_hemophilia_AAV/AAVengeR/outputs/canFam3/sampleReads')
  system(paste0('cp AAVengeR/outputs/canFam3/sampleReads/', 
                 sampleReads[grepl(x$sample, sampleReads) & grepl('\\.breakReads\\.', sampleReads)], ' SRA/',
                 x$sample, '.R1.fasta'))
  
  o <- readFasta(paste0('SRA/', x$sample, '.R1.fasta'))
  write(paste0('@', as.character(o@id), '\n', as.character(o@sread), '\n+\n',  unlist(lapply(width(o), function(x) paste0(rep('I', x), collapse = '')))),
        file = paste0('SRA/', x$sample, '.R1.fastq'))

  system(paste0('cp AAVengeR/outputs/canFam3/sampleReads/', 
                sampleReads[grepl(x$sample, sampleReads) & grepl('\\.virusReads\\.', sampleReads)], ' SRA/',
                x$sample, '.R2.fasta'))
  
  
  o <- readFasta(paste0('SRA/', x$sample, '.R2.fasta'))
  write(paste0('@', as.character(o@id), '\n', as.character(o@sread), '\n+\n',  unlist(lapply(width(o), function(x) paste0(rep('I', x), collapse = '')))),
        file = paste0('SRA/', x$sample, '.R2.fastq'))
  
}))

system('gzip SRA/*.fastq')
system('rm SRA/*.fasta')
write.xlsx(d, file = 'SRA/sampleData.xlsx')
