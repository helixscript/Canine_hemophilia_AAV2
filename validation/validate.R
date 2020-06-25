library(ShortRead)
library(dplyr)
library(readr)
library(parallel)
library(yaml)
library(ggplot2)
options(stringsAsFactors = FALSE)

invisible(lapply(list.dirs()[!grepl('internal', list.dirs())], function(workDir){
if(workDir == '.') return()
setwd(workDir)
message('\n\nStarting ', workDir, '\n\n')
### if(file.exists('done')) return()
#--------------------------------------------------------------------------------------------------

config <- read_yaml('config.yml')
config$blat.path <- '/home/everett/ext/blat'
config$megahit.path  <- '/home/everett/ext/megahit/bin/megahit'
config$blastBin.path <- '/home/everett/ext/blast+/bin'
config$vectorDB.blat.path <- '/home/everett/projects/canine_hemophilia_AAV/data/sequenceData/vectorSeqs_ITR_to_ITR.2bit'



shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}

syncReads <-function(...){
  arguments <- list(...)
  
  # Create a list of read IDs common to all read arguments.
  n <- Reduce(base::intersect, lapply(arguments, names))
  
  lapply(arguments, function(x){ 
    x <- x[names(x) %in% n]; 
    x[order(names(x))]
  })
}

parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble())
  b <- read_delim(f, delim = '\t', col_names = FALSE, col_types = cols())
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}

R1 <- readFastq(config$R1.path)
R2 <- readFastq(config$R2.path)

R1 <- trimTailw(R1, 2, '5', 5)
R2 <- trimTailw(R2, 2, '5', 5)

R1 <- shortRead2DNAstringSet(R1)
R2 <- shortRead2DNAstringSet(R2)

r <- syncReads(R1, R2)
R1 <- r[[1]]; R2 <- r[[2]]

cluster <- makeCluster(config$CPUs)
clusterExport(cluster, 'config', envir = environment())

if(! dir.exists('tmp')) dir.create('tmp')
invisible(file.remove(list.files('tmp', full.names = TRUE))) # Clear out previous runs.

parallelBlat <- function(x, CPUs, dbConfigKey, label = 'x'){
  invisible(parLapply(cluster, split(x, ntile(1:length(x), CPUs)), function(x){
    library(ShortRead)
    f <- paste0(tempfile(pattern = 'tmp', tmpdir = './tmp'), '.', label)
    writeFasta(x, file = paste0(f, '.fasta'))
    system(paste0(config$blat.path, ' ', config[[dbConfigKey]], ' ', paste0(f, '.fasta'), ' ', paste0(f, '.psl'), ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead'))
    file.remove(paste0(f, '.fasta'))
   }))
}

# Blast all unique reads against the vector sequences.
parallelBlat(R1, config$CPUs, 'vectorDB.blat.path', label = 'R1')
parallelBlat(R2, config$CPUs, 'vectorDB.blat.path', label = 'R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R1.psl', full.names = TRUE), collapse = ' '), ' > R1.vector.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R2.psl', full.names = TRUE), collapse = ' '), ' > R2.vector.psl'))

R1.vector.blat <- parseBLAToutput('R1.vector.psl')
R2.vector.blat <- parseBLAToutput('R2.vector.psl')

# Reduce the reads to only reads that partially align to the vectors.
R1 <- R1[names(R1) %in% R1.vector.blat$qName]
R2 <- R2[names(R2) %in% R2.vector.blat$qName]
writeFasta(R1, file = 'R1.vectorAligned.fasta')
writeFasta(R2, file = 'R2.vectorAligned.fasta')

r <- syncReads(R1, R2)
writeFasta(r[[1]], file = 'R1.vectorAligned.paired.fasta')
writeFasta(r[[2]], file = 'R2.vectorAligned.paired.fasta')

# Build contigs with different min-count values using both paired-end and single read approaches.
contigs <- lapply(c(2, 3, 4), function(n){
  dir <- paste0('MEGAHIT_pairedEnd_output', n)
  
  system(paste0(config$megahit.path, ' -1 R1.vectorAligned.paired.fasta -2 R2.vectorAligned.paired.fasta -o ', dir, 
                ' --min-count ', n, ' --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))

  c1 <- readFasta(file.path(dir, 'final.contigs.fa'))
  c1@id <- BStringSet(paste0(as.character(c1@id), '_pairedEnd_minCount_', n))
  unlink(dir, recursive = TRUE)

   dir <- paste0('MEGAHIT_singleEnd_output', n)
   system(paste0(config$megahit.path, ' -r R1.vectorAligned.fasta,R2.vectorAligned.fasta -o ', dir,
                 ' --min-count ', n, ' --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))
  
   c2 <- readFasta(file.path(dir, 'final.contigs.fa'))
   c2@id <- BStringSet(paste0(as.character(c2@id), '_singleEnd_minCount_', n))
   unlink(dir, recursive = TRUE)

  list('pairedEnd' = c1, 'singleEnd' = c2)
})


se <- lapply(contigs, '[[', 'singleEnd')
se <- se[which(unlist(lapply(se, length)) > 0)]

pe <- lapply(contigs, '[[', 'pairedEnd')
pe <- pe[which(unlist(lapply(pe, length)) > 0)]

se <- Reduce('append', se)
pe <- Reduce('append', pe)
r  <- Reduce('append', (list(se, pe)))
r  <- r[width(r) >= 50]
ids <- gsub('\\s', '_', as.character(r@id))
r   <- sread(r)
names(r) <- ids
r  <- unique(r)
r  <- r[order(width(r), decreasing = TRUE)]
writeFasta(r, file = 'contigs.fasta')

system(paste0(config$blastBin.path, '/blastn -db  ../../data/sequenceData/plasmids -query contigs.fasta -word_size 5 -evalue 100 -out contigs.blast -outfmt 6'))

tbl <- read_delim('contigs.blast', '\t', col_names = FALSE)
names(tbl) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


vectorFeatures <- read.csv('../../data/sequenceData/plasmidsFeatures.csv', header = TRUE)


d <- unlist(strsplit(getwd(), '/'))
d <- unlist(strsplit(d[length(d)], '_'))
if(d[1] %in% c('M50', 'M06', 'M66')){
  tbl <- subset(tbl, sseqid == 'SC')
  vectorFeatures <- subset(vectorFeatures, vector == 'SC')
} else {
  tbl <- subset(tbl, sseqid != 'SC')
  vectorFeatures <- subset(vectorFeatures, vector != 'SC')
  vectorFeatures$feature <- paste0(vectorFeatures$vector, '/', vectorFeatures$feature)
}


tbl$contigLength = as.integer(unlist(lapply(stringr::str_match_all(tbl$qname, '_len=(\\d+)'), '[', 2)))
tbl <- subset(tbl, pident >= 95 & gapopen == 0)
tbl$strand <- ifelse(tbl$sstart > tbl$send, '-', '+')
tbl$contigStart <- ifelse(tbl$qstart > tbl$qend, tbl$qend, tbl$qstart)
tbl$contigEnd   <- ifelse(tbl$qstart > tbl$qend, tbl$qstart, tbl$qend)
tbl$vectorStart <- ifelse(tbl$sstart > tbl$send, tbl$send, tbl$sstart)
tbl$vectorEnd   <- ifelse(tbl$sstart > tbl$send, tbl$sstart, tbl$send)
tbl$width <- tbl$contigEnd - tbl$contigStart
tbl <- select(tbl, qname, contigStart, contigEnd, vectorStart, vectorEnd, strand, sseqid, width, contigLength)
tbl$n <- 1:nrow(tbl)
tbl$contigID <- paste0('Contig-', group_indices(tbl, qname), '-', 
                       stringr::str_extract(tbl$qname, '(single|paired)End'), '-',
                       unlist(lapply(stringr::str_match_all(tbl$qname, 'len=(\\d+)'), '[[', 2)), 'NT')


tbl$mappID <- paste(tbl$contigID, tbl$sseqid, tbl$vectorStart, tbl$vectorEnd)

save(list = ls(all.names = TRUE), file = "dev.RData", envir = environment())

browser()

# 	Contig 3 HC 3874 3952
r <- bind_rows(lapply(split(tbl, tbl$contigID), function(contig){
       q <- IRanges(start = contig$vectorStart, end=contig$vectorEnd, names=contig$sseqid)

       r <- arrange(bind_rows(lapply(split(q, 1:length(q)), function(x){
              vectorFeatures <- subset(vectorFeatures, vector == names(x))
              s <- IRanges(start = vectorFeatures$start, end=vectorFeatures$end, names = vectorFeatures$feature)
              o <- data.frame(IRanges::findOverlaps(x, s))

              bind_rows(lapply(o$subjectHits, function(y){
                i <- data.frame(IRanges::intersect(x, s[y]))
                i$feature <- names(s[y])
                i$contigID <- contig$contigID[1]
                i$mappID <- paste(contig$contigID[1], names(x), start(x), end(x))
                i
              }))
            })), start)
  
      r[!duplicated(paste(r$start, r$end, r$feature)),]
    }))

o <- left_join(tbl, r, by = 'mappID')
o$start2 <- o$start - (o$vectorStart - o$contigStart)
o$end2   <- o$end - (o$vectorStart - o$contigStart)
o <- o[!is.na(o$start2),]

invisible(lapply(split(o, o$contigID.x), function(x){
  x <- arrange(x,desc(width.x), start2) %>% mutate(y = 1:n())
  x$feature <- factor(x$feature, levels = vectorFeatures[vectorFeatures$feature %in% x$feature,]$feature)
  colors <- paste0('#', vectorFeatures[vectorFeatures$feature %in% x$feature,]$color)
  title <- paste0(paste0(d, collapse = '_'), '_', x$contigID.x[1])
  p <- ggplot() + 
  ggtitle(title) +
  labs(x = 'Contig position') +
  theme_bw() +
  geom_text(data = x, aes(end2+5, y, label = feature), hjust = 'left', size = 2)+
  scale_color_manual(values = colors) +
  xlim(c(0, max(x$contigLength+100))) +
  geom_vline(xintercept = x$contigLength[1], color = 'gray50', linetype = 'dotted') +
  geom_vline(xintercept = 0, color = 'gray50', linetype = 'dotted') +
  scale_linetype_manual(values = c( 'dotted', 'solid')) +
  geom_segment(aes(x = start2, y = y, xend = end2, yend = y, color = feature, linetype = strand), data = x, size = 1.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  ggsave(p, file = paste0(title, '.pdf'))
}))



# Blast all unique reads against the genome.
invisible(file.remove(list.files('tmp', full.names = TRUE))) 
parallelBlat(R1, config$CPUs, 'genomeDB.blat.path', label = 'R1')
parallelBlat(R2, config$CPUs, 'genomeDB.blat.path', label = 'R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R1.psl', full.names = TRUE), collapse = ' '), ' > R1.genome.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R2.psl', full.names = TRUE), collapse = ' '), ' > R2.genome.psl'))

stopCluster(cluster)

R1.blat <- parseBLAToutput('R1.genome.psl')
R2.blat <- parseBLAToutput('R2.genome.psl')

R1.blat.flankingDownStream <- subset(R1.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingDownStream <- subset(R2.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingUpStream   <- subset(R2.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))
R1.blat.flankingUpStream   <- subset(R1.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))


downStream.bridgingReads <- c(R1[names(R1) %in% R1.blat.flankingDownStream$qName],
                              R2[names(R2) %in% R2.blat.flankingDownStream$qName])

upStream.bridgingReads <- c(R1[names(R1) %in% R1.blat.flankingUpStream$qName],
                               R2[names(R2) %in% R2.blat.flankingUpStream$qName])

invisible(file.remove(list.files('tmp', full.names = TRUE))) 

writeFasta(downStream.bridgingReads, file = 'downStream.bridgingReads.fasta')
writeFasta(upStream.bridgingReads,   file = 'upStream.bridgingReads.fasta')


system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' ', config$R1.path, ' ', config$R2.path, ' > R1R2.vector.sam'))
system('samtools view -S -b R1R2.vector.sam > R1R2.vector.bam')     
system('samtools sort -f R1R2.vector.bam R1R2.vector.sorted.bam')                          
system('samtools index R1R2.vector.sorted.bam') 
file.remove(c('R1R2.vector.sam', 'R1R2.vector.bam'))


if(length(downStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' downStream.bridgingReads.fasta > downStream.bridgingReads.sam'))
  system('samtools view -S -b downStream.bridgingReads.sam > downStream.bridgingReads.bam')     
  system('samtools sort -f downStream.bridgingReads.bam downStream.bridgingReads.sorted.bam')                          
  system('samtools index downStream.bridgingReads.sorted.bam')                    
  file.remove(c('downStream.bridgingReads.sam', 'downStream.bridgingReads.bam'))

  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' downStream.bridgingReads.fasta > downStream.bridgingReads.vector.sam'))
  system('samtools view -S -b downStream.bridgingReads.vector.sam > downStream.bridgingReads.vector.bam')     
  system('samtools sort -f downStream.bridgingReads.vector.bam downStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index downStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('downStream.bridgingReads.vector.sam', 'downStream.bridgingReads.vector.bam'))
}


if(length(upStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' upStream.bridgingReads.fasta > upStream.bridgingReads.sam'))
  system('samtools view -S -b upStream.bridgingReads.sam > upStream.bridgingReads.bam')     
  system('samtools sort -f upStream.bridgingReads.bam upStream.bridgingReads.sorted.bam')                          
  system('samtools index upStream.bridgingReads.sorted.bam')                    
  file.remove(c('upStream.bridgingReads.sam', 'upStream.bridgingReads.bam'))
  
  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' upStream.bridgingReads.fasta > upStream.bridgingReads.vector.sam'))
  system('samtools view -S -b upStream.bridgingReads.vector.sam > upStream.bridgingReads.vector.bam')     
  system('samtools sort -f upStream.bridgingReads.vector.bam upStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index upStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('upStream.bridgingReads.vector.sam', 'upStream.bridgingReads.vector.bam'))
}

unlink('tmp', recursive = TRUE)
dir.create('alignments')
system('mv *.bam alignments/')
system('mv *.bai alignments/')
invisible(file.remove(list.files(pattern = '*.psl')))

write(c('new',
      'genome canFam3',
      paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.sorted.bam'),
      paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.sorted.bam'),
      paste0('snapshotDirectory ', getwd()),
      paste0('goto ', config$siteChromosome, ':', config$sitePosition-1000, '-',  config$sitePosition+1000),         
      'sort position',
      'collapse',
      'snapshot genome.bridgingReads.png',
      'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto HeavyChainVector:1-4042',         
        'sort position',
        'collapse',
        'snapshot vector.heavy.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto LightChainVector:1-3904',         
        'sort position',
        'collapse',
        'snapshot vector.light.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto SingleChainVector:1-5411',         
        'sort position',
        'collapse',
        'snapshot vector.single.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')




write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto HeavyChainVector:1-4042',         
        'sort position',
        'collapse',
        'snapshot vector.heavy.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto LightChainVector:1-3904',         
        'sort position',
        'collapse',
        'snapshot vector.light.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')



write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto SingleChainVector:1-5411',         
        'sort position',
        'collapse',
        'snapshot vector.single.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


invisible(file.remove(list.files('tmp', full.names = TRUE))) 
dir.create('workFiles')
system('mv batchFile alignments *.blast *.fasta workFiles')

write('done', file = 'done')
#--------------------------------------------------------------------------------------------------
setwd('..')
}))





