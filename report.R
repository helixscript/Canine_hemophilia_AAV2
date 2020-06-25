library(gt23)
library(parallel)
library(BSgenome.Cfamiliaris.UCSC.canFam3)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ShortRead)
library(BSgenome)
library(BSgenome.Cfamiliaris.UCSC.canFam3)
library(scales)

Rscript_path <- '/home/opt/R-3.4.0/bin/Rscript'
createColorPalette <- function(n) grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }
blastBin <- '/home/everett/ext/blast+/bin' 

# Create experimental data save file from the AAVengerR output file.
if(! file.exists('data/data.RData')){
  load('data/sites.RData')
  intSites           <- makeGRangesFromDataFrame(sites, start.field = 'position', end.field = 'position', keep.extra.columns = TRUE)
  intSites$refGenome <- 'canFam3'
  intSites$patient   <- intSites$subject
  intSites           <- annotateIntSites(intSites)
  intSites$subject   <- sub('^p', '', intSites$subject)
  save(list = ls(all.names = TRUE), file = 'data/data.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
} else {
  load('data/data.RData')
}


# Define convenient data subsets and values.
intSites <- data.frame(intSites)

# Remove multihits.
intSites <- subset(intSites, seqnames != 'chrUn')


# Manually review cases where the same intSite position was found more than once.
# and decide if any sites are likely instances of cross contamination based on read counts. 
#
# o <- table(intSites$posid)
# o <- o[o > 1]
# posids <- names(o)[-grep('chrX[+-]12[23]',names(o))]
# o <- subset(intSites, posid %in% posids)
# View(dplyr::select(o, posid, subject, reads, sample, estAbund))

# Removed apparent sample contaminations
i <- c(which(intSites$posid == 'chr25-34548677' & intSites$sample == 'GTSP2178'), 
       which(intSites$posid == 'chr25+34592325' & intSites$sample == 'GTSP2175'))

intSites <- intSites[-i,]



# Create UCSC browser session file and cp it to the Bushman lab web server.
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=canFam3&hgt.customText=http://microb120.med.upenn.edu/UCSC/sabatino/AAVengeR.ucsc
system(paste0('scp AAVengeR/outputs/canFam3_final/AAVengeR.ucsc  microb120:/usr/share/nginx/html/UCSC/sabatino/'))



# Read in information about the sites that were tested in the laboratory.
d <- read.csv('data/siteValidations.csv', header = TRUE)
intSites$validation <- NA
intSites[apply(d, 1, function(x){ which(intSites$sample == x[1] & intSites$posid == x[2]) }),]$validation <- as.character(d$highestMethods)



# Dual-dections. Use the output below to create data/duelDetections.csv file to merge together duel detections.
#--------------------------------------------------------------------------------------------------
# d <- select(intSites, seqnames, strand, start, sample, subject, estAbund, posid, nearestFeature, nearestFeatureDist) 
# distFromSite <- 100
# d <- subset(d, nearestFeature != 'F8')
# dualDetection <- lapply(split(d, paste(d$sample, d$posid)), function(x){
#   o <- subset(d, sample == x$sample & start >= (x$start - distFromSite) & start <= (x$start + distFromSite))
#   if(n_distinct(o$strand) == 2){
#     o$distBetween <- abs(x$start - o$start)
#     return(select(o, posid, distBetween, sample, estAbund, nearestFeature, nearestFeatureDist))
#   } else {
#     return()
#   }
# })
# dualDetection[sapply(dualDetection, is.null)] <- NULL


# Read in the duel detection table created with the commented code block above 
# and merge sites into the mean position.
dualDetection <- read.csv('data/duelDetections.csv', header = TRUE, comment.char = "#")
for(x in 1:nrow(dualDetection)){
  x <- dualDetection[x,]
  a <- as.integer(unlist(strsplit(as.character(x$pos1), '[\\+\\-]'))[2])
  b <- as.integer(unlist(strsplit(as.character(x$pos2), '[\\+\\-]'))[2])
  newPos <-  floor(mean(c(a,b)))
  
  i1 <- which(intSites$posid == x$pos1 & intSites$sample == x$sample)
  i2 <- which(intSites$posid == x$pos2 & intSites$sample == x$sample)
  message(i1, ' ', x$pos1, ' ', i2, ' ', x$pos2, ' New position: ', newPos, ' Diff: ', abs(a-b))
  if(length(c(i1, i2)) != 2){
    browser()
    stop('Duel detection error.')
  }
  
  intSites[i1,]$start    <- newPos
  intSites[i1,]$end      <- newPos
  intSites[i1,]$reads    <- sum(intSites[c(i1, i2),]$reads)
  intSites[i1,]$strand   <- '*'
  intSites[i1,]$posid    <- paste0(intSites[i1,]$seqnames, intSites[i1,]$strand, newPos)
  intSites[i1,]$estAbund <- max(c(x$estAbund1, x$estAbund2))
  v <- intSites[c(i1, i2),]$validation
  v <- v[!is.na(v)]
  if(length(v) == 0) v <- NA
  intSites[i1,]$validation <- paste0(v, collapse = ',')
  intSites <- intSites[-i2,]; message('Removed ', i2)
}


# Add sample details to intSite data.
d <- read.table('data/sampleDetails.tsv', sep = '\t', header = TRUE)
intSites <- left_join(intSites, d, by = 'sample')
intSites$posid2 <- paste(intSites$subject, intSites$posid)


# Create a series of subsets of the data for analysis.
expIntSites <- filter(intSites, ! subject %in% c('HO2', 'M12'))
expIntSites.noF8  <- filter(expIntSites, nearestFeature != 'F8')
expIntSites.abund.gte2  <- filter(expIntSites, estAbund >= 2)
expIntSites.abund.gte2.noF8 <- filter(expIntSites.abund.gte2, nearestFeature != 'F8')
expIntSites.noF8.inTU       <- subset(expIntSites.noF8, nearestFeatureDist == 0)
expIntSites.noF8.inTU.inOncogene <- subset(expIntSites.noF8, nearestOncoFeatureDist == 0)
expIntSites.noF8.nearOncogene    <- subset(expIntSites.noF8, abs(nearestOncoFeatureDist) <= 25000)
expIntSites.noF8.percentInTU     <- n_distinct(expIntSites.noF8.inTU$posid) / n_distinct(expIntSites.noF8$posid)
expIntSites.noF8.inTU.percentInOnco <- n_distinct(expIntSites.noF8.inTU.inOncogene$posid) / n_distinct(expIntSites.noF8.inTU$posid)
expIntSites.noF8.percentNearOnco    <- n_distinct(expIntSites.noF8.nearOncogene$posid) / n_distinct(expIntSites.noF8$posid)
expSitesInTU.inOnco  <- n_distinct(expIntSites.noF8.inTU.inOncogene$posid) / n_distinct(expIntSites.noF8.inTU$posid)


# Same clones across samples
#--------------------------------------------------------------------------------------------------
clonesAcrossSamples <- 
  group_by(expIntSites.noF8, posid) %>% 
  mutate(n = n_distinct(sample)) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>%
  dplyr::select(posid, subject, sample, sampleName, estAbund) %>%
  arrange(posid)




# Write out vector sequences for supp. materials.
a <- readFasta('data/sequenceData/singleChainVector.fasta')
b <- readFasta('data/sequenceData/lightChainVector.fasta')
c <- readFasta('data/sequenceData/heavyChainVector.fasta')
d <- Reduce('append', list(a, b, c))
writeFasta(d, 'tables_and_figures/vectorSequences.fasta')



# Tables
#--------------------------------------------------------------------------------------------------
sitesTable <- mutate(expIntSites.abund.gte2.noF8, ITR_reminant = nchar(ltrRepSeq)) %>%
  select(subject, timePoint, posid, estAbund, ITR_reminant, nearestFeature, 
         nearestFeatureDist, nearestOncoFeature, nearestOncoFeatureDist, validation) %>%
  arrange(desc(estAbund))

names(sitesTable) <- c('Dog',  'Timepoint', 'Position', 'Abundance', 'ITR Remnant length', 
                       'Nearest Gene', 'Distance', 'Nearest oncogene', 'Distance', 'Validated')

openxlsx::write.xlsx(sitesTable, file = 'tables_and_figures/sites_abund_gte2.xlsx')





# Clustering schematics
#--------------------------------------------------------------------------------------------------

createSiteClusterDiagram <- function(d, genes){
  d <- dplyr::select(d, subject, start, estAbund) %>%
       dplyr::arrange(start, estAbund) %>%
       dplyr::mutate(subject = factor(subject, levels = c('H19', 'J60', 'Linus', 'M06', 'M50', 'M66')))
  d$abundBin <- '1 - 4 cells'
  d$abundBin <- ifelse(d$estAbund >= 5 & d$estAbund <= 15, '5 - 15 cells', d$abundBin)
  d$abundBin <- ifelse(d$estAbund > 15, '> 15 cells', d$abundBin)
  d$abundBin <- factor(d$abundBin, levels = c('1 - 4 cells', '5 - 15 cells', '> 15 cells'))
  
  colors <- colorRampPalette(brewer.pal(12, "Paired"))(6)
  names(colors) <- c('H19', 'J60', 'Linus', 'M06', 'M50', 'M66')
  
  
  ggplot(d, aes(start, 0, fill = subject, shape = abundBin)) +
  theme_bw() +
  geom_hline(yintercept = 0, size = 0.5) +
  scale_fill_manual(values = colors) + 
  scale_shape_manual(name = 'Abundance', values = c(21, 22, 24)) +
  geom_rect(data = genes, inherit.aes = FALSE, mapping=aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5), 
            color="black", alpha=0.5, size = 0.5) + 
  geom_point(size = 2) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
}


EGR2cluster <- createSiteClusterDiagram(subset(expIntSites.noF8, seqnames == 'chr4' & start >= 14512504 & end <= 14838823),
                                        data.frame(start = c(14705384, 14698266), end = c(14709742, 14701927), 
                                                   names = c('EGR2', 'ADO'), fill = c('gray40', 'gray60')))
ggsave(EGR2cluster, file = 'tables_and_figures/EGR2cluster.pdf', units = 'in', width = 10)

EGR3cluster <- createSiteClusterDiagram(subset(expIntSites.noF8, seqnames == 'chr25' & start >= 34587341 & end <= 34597415),
                                        data.frame(start = 34591290, end = 34597144, 
                                                   names = 'EGR3', fill = 'gray40'))
ggsave(EGR3cluster, file = 'tables_and_figures/EGR3cluster.pdf', units = 'in', width = 10)

CCND1_cluster <- createSiteClusterDiagram(subset(expIntSites.noF8, seqnames == 'chr18' & start >= 48461911 & end <= 48829811),
                                          data.frame(start = 48501886, end = 48509663, 
                                                     names = 'CCND1', fill = 'gray40'))
ggsave(CCND1_cluster, file = 'tables_and_figures/CCND1_cluster.pdf', units = 'in', width = 10)





# Identify read fragments which cross exon boundaries when the experimental data
# is aligned to the single chain vector sequence. 
#--------------------------------------------------------------------------------------------------

singleChainVectorSeq <- as.character(readFasta('data/sequenceData/singleChainVector.fasta')@sread)

e <- dplyr::select(data.frame(subset(gt23::canFam3.refSeqGenesGRanges.exons, name2 == 'F8')), start, end) %>%
     dplyr::mutate(exonStarts = as.integer(start), exonEnds = as.integer(end))

e <- bind_rows(lapply(split(e, nrow(e):1), function(x){
  x$seq <- getSeq(BSgenome.Cfamiliaris.UCSC.canFam3, 'chrX', x$exonStarts, x$exonEnds, as.character = TRUE)
  x$revCompSeq <- as.character(reverseComplement(DNAString(x$seq)))
  o <- matchPattern(x$revCompSeq, singleChainVectorSeq, max.mismatch = 3)
  if(length(start(o))){
    x$start <- start(o)
    x$end <- end(o)
  }
  x
}))

# Manually add exons 14 and 26 because its center has been deleted in the vector.
e[14,]$start <- 2835
e[14,]$end   <- 3257
e[26,]$start <- 4940
e[26,]$end   <- 5094

f <- readRDS('AAVengeR/outputs/singleChainVector/readFrags.rds')
f <- f[! duplicated(paste(f$start, f$strand, f$end)),]

# Identify reads which cross exon boundaries and exceed the boundaries by >= 10 NT.
fe <- bind_rows(lapply(split(e, 1:nrow(e)), function(x){
  bind_rows(f[f$start < x$start & f$end > x$start & x$start - f$start >= 10 & f$end - x$start >= 10,],
            f[f$start < x$end & f$end > x$end & x$end - f$start >= 10 & f$end - x$end >= 10,])
})) 

# Remove duplicate fragments.
fe <- fe[! duplicated(paste(fe$start, fe$strand, fe$end)),] %>% arrange(start, reads) %>% mutate(y = 1:n())

# 71.9% of fragments cross an exon boundary
message(sprintf("%.1f%%", (nrow(fe) / nrow(f))*100), ' of fragments cross an exon boundary')

fe$reads <- ifelse(fe$reads >= 50, 50, fe$reads)

vectorExonCrossings <-                       
  ggplot() +
  theme_bw() +
  labs(x = 'Vector position', y = 'Fragments') +
  geom_segment(aes(x = start, y = 1, xend = start, yend = max(fe$y)), data = e, linetype = 'dashed', color = 'gray75') +
  geom_segment(aes(x = end, y = 1, xend = end, yend = max(fe$y)), data = e, linetype = 'dashed', color = 'gray75') +
  geom_segment(aes(x = start, y = y, xend = end, yend = y, colour = reads), data = fe) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_colour_gradient2(low  = "dodgerblue", mid  = "gold", high = "red", midpoint = 25) +
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(vectorExonCrossings, file = 'tables_and_figures/vectorExonCrossings.pdf', units = 'in', width = 8)




# Create site counts considering sites within vectors.
#--------------------------------------------------------------------------------------------------

# Load AAVenger sites from runs where the vector was used as the reference geneome.
vectorSitesSingle <- new.env()
vectorSitesLight  <- new.env()
vectorSitesHeavy  <- new.env()

load('AAVengeR/outputs/singleChainVector/sites.RData', envir = vectorSitesSingle)
load('AAVengeR/outputs/lightChainVector/sites.RData',  envir = vectorSitesLight)
load('AAVengeR/outputs/heavyChainVector/sites.RData',  envir = vectorSitesHeavy)


# Here we loop through each sample and compare the number of sites found outside of F8 to the
# number of sites found within the vector. We are assuming all sites within F8 are vector sites 
# based on previous reusults. The split chain vectors needs to be handled differently
# than then single chain vector because the split chain vectors share 5' ITR, promoter, polyA and 3' ITR sequences.
# For the split vectors, for each sample, we use the average number of sites from the shared regions. 
# The number of within vector sites are likely a lower bound estimate because the position standardization 
# algorithm is likely combining unique sites because some many sites are located close to one another.

siteCounts <- bind_rows(lapply(split(expIntSites.noF8, expIntSites.noF8$sample), function(x){
  #vectorSites <- tibble()
  expType <- as.character(unique(x$experimentType))
  
  if(expType == 'SingleChain'){
    vectorSites <- n_distinct(subset(vectorSitesSingle$sites, sample == x$sample[1])$posid)
    vectorSiteList <- unique(subset(vectorSitesSingle$sites, sample == x$sample[1])$posid)
    
  } else if(expType == 'SplitChain'){
    
    # lightBlat 1115  
    # heavyBlat 1100
    lightTransgeneSites <- subset(vectorSitesLight$sites, sample == x$sample[1] & position >= 1122 & position <= (1122 + 2385))
    lightBackboneSites  <- subset(vectorSitesLight$sites, sample == x$sample[1] & ! posid %in% lightTransgeneSites$posid)
    heavyTransgeneSites <- subset(vectorSitesHeavy$sites, sample == x$sample[1] & position >= 1106 & position <= (1106 + 2546))
    heavyBackboneSites  <- subset(vectorSitesHeavy$sites, sample == x$sample[1] & ! posid %in% heavyTransgeneSites$posid)
    
    # Visual check that we are parsing the vector correctly.
    # browser()
    # ggplot(bind_rows(tibble(pos = heavyTransgeneSites$position, source = 'transgene'), 
    #                  tibble(pos = heavyBackboneSites$position, source = 'backbone')), aes(pos, 1, color = source)) + geom_point()
    
    vectorSites <- floor(mean(c(n_distinct(lightBackboneSites$posid), n_distinct(heavyBackboneSites$posid)))) + 
      n_distinct(lightTransgeneSites$posid) + n_distinct(heavyTransgeneSites$posid)
    
    vectorSiteList <- unique(c(lightBackboneSites$posid, heavyBackboneSites$posid, lightTransgeneSites$posid, heavyTransgeneSites$posid))
    
  } else {
    stop('Experiment selection error.')
  }
  
  tibble(dog = x$subject[1], sample = x$sample[1], liverLobe = stringr::str_extract(x$sampleName[1], '\\d+$'), 
         timePoint = x$timePoint[1],expType = expType, VCN = x$VCN[1], sampleMass = x$sampleMass[1],
         nSites = n_distinct(x$posid), nVectorSites = vectorSites, 
         nSitesMassCorrected = n_distinct(x$posid)/x$sampleMass[1],
         nSitesVCNCorrected = n_distinct(x$posid)/x$VCN[1],
         percentSitesInVector =  nVectorSites / (n_distinct(x$posid) + nVectorSites)*100)
}))

siteTable <- arrange(siteCounts, expType)
openxlsx::write.xlsx(siteTable, file = 'tables_and_figures/siteCountTable.xlsx')


# Manuscript stats
#--------------------------------------------------------------------------------------------------

# Unique sites in genome.
n_distinct(expIntSites.noF8$posid)

(n_distinct(subset(expIntSites.noF8, estAbund >= 2)$posid) / n_distinct(expIntSites.noF8$posid))*100
(n_distinct(subset(expIntSites.noF8, estAbund >= 3)$posid) / n_distinct(expIntSites.noF8$posid))*100
(n_distinct(subset(expIntSites.noF8, estAbund >= 4)$posid) / n_distinct(expIntSites.noF8$posid))*100
(n_distinct(subset(expIntSites.noF8, estAbund >= 5)$posid) / n_distinct(expIntSites.noF8$posid))*100


# Test for differences in the number of sites recovered between single chain and split chain dogs.
#--------------------------------------------------------------------------------------------------

# Difference between number of sites recovered from single chain vs split chain (sample level p-val)
wilcox.test(subset(siteCounts, expType == 'SingleChain')$nSites, subset(siteCounts, expType == 'SplitChain')$nSites)$p.val

# Difference between number of sites recovered from single chain vs split chain (sample level, mass-corrected p-val)
wilcox.test(subset(siteCounts, expType == 'SingleChain')$nSitesMassCorrected, subset(siteCounts, expType == 'SplitChain')$nSitesMassCorrected)$p.val

# Difference between number of sites recovered from single chain vs split chain (sample level, VCN-corrected p-val)
wilcox.test(subset(siteCounts, expType == 'SingleChain')$nSitesVCNCorrected, subset(siteCounts, expType == 'SplitChain')$nSitesVCNCorrected)$p.val

# Difference between percent sites in vector from single chain vs split chain (sample level p-val)
wilcox.test(subset(siteCounts, expType == 'SingleChain')$percentSitesInVector, subset(siteCounts, expType == 'SplitChain')$percentSitesInVector)$p.val







# Gene enrichment from abundant clones.
#--------------------------------------------------------------------------------------------------

if(! file.exists('tables_and_figures/KEGG_enrichment.xlsx')){
  library(STRINGdb)
  string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0)
  o <- subset(expIntSites.noF8, abs(nearestFeatureDist) <= 50000 & estAbund >= 2)
  d <- data.frame(gene = unlist(strsplit(o$nearestFeature, ',')))

  # Change gene names to older symbols so that stringDB will recognize them.
  d$gene <- gsub('LTO1', 'ORAOV1', d$gene)
  d$gene <- gsub('CFAP299', 'C4orf22', d$gene)
  d$gene <- gsub('MINDY2', 'FAM63B', d$gene)
  d$gene <- gsub('SUGCT', 'C7orf10', d$gene)
  d$gene <- gsub('UBE2E4P', 'UbcM2', d$gene)
  d$gene <- gsub('ADGRL2', 'LPHH1', d$gene)
  d <- string_db$map(d, "gene", removeUnmappedRows = FALSE)
  e <- string_db$get_enrichment(d$STRING_id, category = 'KEGG', methodMT = "fdr", iea = FALSE)
  e <- subset(e, pvalue_fdr <= 0.05)
  openxlsx::write.xlsx(e, file = 'tables_and_figures/KEGG_enrichment.xlsx')
}




# Sample tables for manuscript supp.
#--------------------------------------------------------------------------------------------------

sampleTable <- group_by(intSites, sampleName) %>%
  summarise(dog = subject[1], sampleID = sample[1], nSites = n_distinct(posid), 
            timePoint = timePoint[1], liverLobe = stringr::str_extract(sampleName[1], '\\d+$'), VCN = VCN[1], inputMass = sampleMass[1]) %>%
  ungroup() %>%
  arrange(desc(dog)) %>%
  dplyr::select(-sampleName) 

openxlsx::write.xlsx(sampleTable, file = 'tables_and_figures/sampleTable.xlsx')



# Abundant clones
#--------------------------------------------------------------------------------------------------
expIntSites.noF8$lobe <- unlist(lapply(str_match_all(expIntSites.noF8$sampleName, '_(\\d+)'), '[', 2))

abundantClones_gte5 <- 
  filter(expIntSites.noF8, estAbund >= 5) %>%
  select(subject, posid, nearestFeature, estAbund, lobe) %>%
  mutate(posid = sub('[^+]+_', '', posid)) %>%
  group_by(subject) %>%
  top_n(n = 15, wt = estAbund) %>%
  ungroup() %>%
  arrange(desc(estAbund)) %>%
  mutate(posid = factor(posid, levels = unique(posid)),
         subject = sub('^p', '', subject),
         subject = factor(subject, levels = rev(c('M06', 'M66', 'M50', 'J60', 'Linus')))) %>%
         {
           ggplot(., aes(posid, estAbund, fill = nearestFeature, label = lobe)) +
             geom_bar(color = 'black', stat='identity', position="dodge", width = 1) + 
             scale_fill_manual(name = 'Nearest Feature', values =   createColorPalette(n_distinct(.$nearestFeature))) +
             facet_grid(subject~., scales = 'free', space = "free", switch = 'x') +
             theme_bw() +
             labs(x = '', y = 'Clonal abundance') +
             coord_flip() +
             geom_text(nudge_y = 5) +
             theme(strip.text.y = element_text(angle = 0, size = 12),
                   strip.text.x = element_text(angle = 0, size = 10),
                   axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
                   legend.title=element_text(size = 14), legend.text=element_text(size = 12),
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
         }

ggsave(abundantClones_gte5, file = 'tables_and_figures/abundantClones_gte5.pdf', height = 10, width = 8, units = 'in')




# Random site tests
#--------------------------------------------------------------------------------------------------

# Create a collection of random sites in the canine genome.
# More sites are requested than sampled because less than the number requested will be returned 
# because sites annoated as N are removed from the returned vector.
# Here we break the random sites into three files so they can be uploaded to GitHub.
if(! file.exists('data/canFam3.randomSites_annotated.1.RData')) {
  set.seed(42)
  randomSites <- sample(createRandomFragments(refGenome = BSgenome.Cfamiliaris.UCSC.canFam3, n = 5500000, fragWidth = 1), 5000000)
  chromosomes <- unlist(lapply(str_split(names(randomSites), ':'), '[', 1))
  position    <- as.integer(unlist(lapply(str_split(unlist(lapply(str_split(names(randomSites), ':'), '[', 2)), '-'), '[', 1)))
  g <- GRanges(seqnames = chromosomes,  IRanges::IRanges(start = position, end = position), strand = '+')
  g <- addPositionID(g)
  g$refGenome  <- 'canFam3'
  g$patient    <- 'random'
  
  randomSites_annotated <- data.frame(annotateIntSites(g, CPUs = 40))
  
  randomSites_annotated$n <- ntile(1:nrow(randomSites_annotated), 3)
  randomSites_annotated.1 <- subset(randomSites_annotated, n == 1)
  randomSites_annotated.2 <- subset(randomSites_annotated, n == 2)
  randomSites_annotated.3 <- subset(randomSites_annotated, n == 3)
  
  save(list = c('randomSites_annotated.1'), file = 'data/canFam3.randomSites_annotated.1.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
  save(list = c('randomSites_annotated.2'), file = 'data/canFam3.randomSites_annotated.2.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
  save(list = c('randomSites_annotated.3'), file = 'data/canFam3.randomSites_annotated.3.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
} else {
  load('data/canFam3.randomSites_annotated.1.RData')
  load('data/canFam3.randomSites_annotated.2.RData')
  load('data/canFam3.randomSites_annotated.3.RData')
  randomSites_annotated <- bind_rows(randomSites_annotated.1, randomSites_annotated.2, randomSites_annotated.3)
  randomSites_annotated <- randomSites_annotated[! duplicated(randomSites_annotated$posid),]
  rm(randomSites_annotated.1, randomSites_annotated.2, randomSites_annotated.3)
}



# ScanStats -- identify clustering of sites 
#--------------------------------------------------------------------------------------------------
scanStats <-  GenomicRanges::GRangesList(lapply(1:50, function(x){
  set.seed(x)
  a <- makeGRangesFromDataFrame(dplyr::select(expIntSites.noF8, seqnames, start, end, strand, nearestFeature, nearestFeatureDist), keep.extra.columns = TRUE)
  r <- dplyr::sample_n(randomSites_annotated, nrow(expIntSites.noF8))
  b <- makeGRangesFromDataFrame(dplyr::select(r, seqnames, start, end, strand, nearestFeature, nearestFeatureDist), keep.extra.columns = TRUE)
  
  s <- gt23::scanStats(a, b, kvals = '10L:15L', 
                       cutpt.filter.expr = "apply(x, 2, quantile, probs = 0.8, na.rm = TRUE)",
                       cutpt.tail.expr = "critVal.target(k, n, target = 10, posdiff = x)")
  
  d <- dplyr::arrange(data.frame(s),seqnames) %>% 
    dplyr::filter(clusterSource == 'A') %>% 
    dplyr::mutate(widthKB = width / 1000) %>%
    dplyr::select(seqnames, start, end, value1, value2, target.min, widthKB) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(widthKB <= 500, target.min <= 0.20)
  
  if (nrow(d) > 0){
    return(as.data.frame(d) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
  } else {
   return(GRanges()) 
  }
}))

clusters <- GenomicRanges::reduce(unlist(scanStats))
clusters$name <- c('EGR2', 'DUSP1', 'ALB', 'CCND1', 'EGR3')

n <- table(unlist(lapply(scanStats, function(x){
       if(length(x) > 0){
         return(clusters[subjectHits(findOverlaps(x, clusters))]$name)
       } else {
         return(NA)
       }
      })))

clusters$hits <- as.integer(unname(n[match(clusters$name, names(n))]))






# Test for enrichment in TUs compared to random sites
#-------------------------------------------------------------------------------------------------

# Create subsets of the random sites to export to CPU clusters.
randomSites_annotated_lite <- select(randomSites_annotated, posid, nearestFeatureDist, nearestOncoFeatureDist)
randomSites_annotated_lite_inTU <- subset(randomSites_annotated_lite, nearestFeatureDist == 0)


if(file.exists('data/random.tests.rds')){
  random.tests <- readRDS('data/random.tests.rds')
} else {
  cluster <- makeCluster(40)
  expSampleSize <- n_distinct(expIntSites.noF8$posid)
  clusterExport(cluster, c('expSampleSize', 'expIntSites.noF8.percentInTU', 'randomSites_annotated_lite'))
  
  random.tests <- bind_rows(parLapply(cluster, 1:1000000, function(x){
    library(dplyr)
    set.seed(x) 
    s  <- sample_n(randomSites_annotated_lite, expSampleSize, replace = FALSE)
    
    tibble(sample = x, 
           sampledPoistionsInTU =  n_distinct(subset(s, nearestFeatureDist == 0)$posid),
           sampledPoistions = n_distinct(s$posid),
           pInTU =  sampledPoistionsInTU / sampledPoistions,
           randomLower = ifelse(pInTU < expIntSites.noF8.percentInTU, TRUE, FALSE))
  }))
  
  saveRDS(random.tests, file = 'data/random.tests2.rds')
  stopCluster(cluster)
}

inTuPval <- sum(random.tests$randomLower) / nrow(random.tests)

bins <- seq(35, 53, by = .222)/100
random.tests$bin <- cut(random.tests$pInTU, breaks = c(0, bins, Inf), labels = FALSE)
random.tests$binLabel <- sprintf("%.1f%%", bins[random.tests$bin]*100)
random.tests.plotData <- group_by(random.tests, binLabel) %>% summarise(n = n()) %>% ungroup()
random.tests.plotData$binLabel <- factor(as.character(random.tests.plotData$binLabel), levels = sprintf("%.1f%%", bins*100))


inTU.randomTest <- 
  ggplot(random.tests.plotData, aes(binLabel, n)) + 
  geom_col() + 
  theme_bw() + 
  scale_y_continuous(label = scales::comma) +
  scale_x_discrete(drop = FALSE, breaks = sprintf("%.1f%%", bins*100)[seq(1, length(bins), by = 5)])  + 
  #geom_segment(aes(x = 65, y = 10000, xend = 65, yend = 1000), color='black', size=1, arrow = arrow(length = unit(0.025, "npc"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14), axis.title=element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")+
  geom_segment(aes(x = 0, xend = 82, y = 0, yend = 0)) +
  labs(x = 'Percent sites bins', y = 'Simulations') + 
  geom_segment(aes(x = 77, y = 7000, xend = 77, yend = 1000), size = 1.1, color = "dodgerblue1", arrow = arrow(length = unit(0.03, "npc")))

ggsave(inTU.randomTest, file = 'tables_and_figures/inTu_vs_randomSites.pdf')




# Test for enrichment in oncogene TUs compared to random TU sites.
#--------------------------------------------------------------------------------------------------

if(file.exists('data/random.inTU.tests.rds')){
  random.inTU.tests <- readRDS('data/random.inTU.tests.rds')
} else {
  cluster <- makeCluster(40)
  expSampleSize <- n_distinct(expIntSites.noF8.inTU$posid)
  
  clusterExport(cluster, c('randomSites_annotated_lite_inTU', 'expSitesInTU.inOnco', 'expSampleSize'))
  
  random.inTU.tests <- bind_rows(parLapply(cluster, 1:1000000, function(x){
    library(dplyr)
    set.seed(x)
    s  <- sample_n(randomSites_annotated_lite_inTU, expSampleSize, replace = FALSE)
    
    tibble(sample = x, 
           sampledPoistionsInOnco =  n_distinct(subset(s, nearestOncoFeatureDist == 0)$posid),
           sampledPoistions = n_distinct(s$posid),
           pInOnco =  sampledPoistionsInOnco / sampledPoistions,
           randomLower = ifelse(pInOnco < expSitesInTU.inOnco, TRUE, FALSE))
  }))
  
  saveRDS(random.inTU.tests, file = 'data/random.inTU.tests.rds')
  stopCluster(cluster)
}

oncoPval <- sum(random.inTU.tests$randomLower) / nrow(random.inTU.tests)

bins <- seq(9, 22, by = .213)/100  # lowed end bin values.
random.inTU.tests$bin <- cut(random.inTU.tests$pInOnco, breaks = c(0, bins, Inf))
random.inTU.tests$binLabel <- sprintf("%.1f%%", bins[random.inTU.tests$bin]*100)
random.inTU.tests.plotData <- group_by(random.inTU.tests, binLabel) %>% summarise(n = n()) %>% ungroup()
random.inTU.tests.plotData$binLabel <- factor(as.character(random.inTU.tests.plotData$binLabel), levels = sprintf("%.1f%%", bins*100))


random.inTU.tests.plotData$color <- factor(ifelse(random.inTU.tests.plotData$binLabel == '18.6%', 1, 0))
inTU.randomOncoGeneTest <- 
  ggplot(random.inTU.tests.plotData, aes(binLabel, n, fill = color)) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = c('gray25', 'dodgerblue1')) +
  scale_y_continuous(label = scales::comma) +
  scale_x_discrete(drop = FALSE, breaks = sprintf("%.1f%%", bins*100)[seq(1, length(bins), by = 5)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14), axis.title=element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  geom_segment(aes(x = 0, xend = 58, y = 0, yend = 0)) +
  labs(x = 'Percent sites bins', y = 'Simulations')

ggsave(inTU.randomOncoGeneTest, file = 'tables_and_figures/inOncoTu_vs_randomSites.pdf')



# VCN vs intSites / sample mass
#--------------------------------------------------------------------------------------------------

plotData <- group_by(expIntSites, sampleName) %>%
  summarise(VCN = VCN[1], 
            sites.mass = n_distinct(posid) / sampleMass[1],
            subject = gsub('_\\S+$|^p', '', sampleName[1]),
            experimentType = experimentType[1]) %>%
  ungroup() 

plotData.lm <- lm(log10(plotData$VCN) ~ plotData$sites.mass)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

VCNvsInts_rSquared <- summary(plotData.lm)$r.squared
VCNvsInts_modelPval <- lmp(plotData.lm)

VCNvsMassIntSites <-
  ggplot(plotData, aes(log10(VCN), sites.mass, fill = subject, shape = experimentType)) +
  theme_bw() +
  geom_smooth(method = "lm", alpha = 0.1, inherit.aes = FALSE, aes(x = log10(VCN), y = sites.mass)) +
  scale_fill_manual(name = 'Dog', values = colorRampPalette(brewer.pal(12, "Paired"))(6)) +
  scale_shape_manual(name = 'Experiment', values = c(21, 22)) +
  geom_point(size = 3, color = 'black') +
  labs(x = 'log10(VCN)', y = 'Integrations / sample mass') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title=element_text(size = 14), legend.text=element_text(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

ggsave(VCNvsMassIntSites, file = 'tables_and_figures/VCNvsMass_sitesRecovered.pdf')


# IntSite genome map
#--------------------------------------------------------------------------------------------------

chromosomeLengths <- sapply(rev(paste0("chr", c(seq(1:38), "X"))),
                            function(x){length(BSgenome.Cfamiliaris.UCSC.canFam3[[x]])},
                            simplify = FALSE, USE.NAMES = TRUE)

intSiteDistributionPlot <- function (d, chromosomeLengths, alpha = 1) 
{
  library(GenomicRanges)
  library(ggplot2)
  library(gtools)
  library(dplyr)
  d <- d[, c("start", "seqnames", "subject")]
  d <- suppressWarnings(bind_rows(d, bind_rows(lapply(names(chromosomeLengths)[!names(chromosomeLengths) %in% 
                                                                                 unique(as.character(d$seqnames))], function(x) {
                                                                                   data.frame(start = 1, seqnames = x)
                                                                                 }))))
  d$seqnames <- factor(d$seqnames, levels = mixedsort(names(chromosomeLengths)))
  
  d <- lapply(split(d, d$seqnames), function(x) {
    lines <- data.frame(x = rep(x$start, each = 2), 
                        y = rep(c(0, 1), nrow(x)), 
                        g = rep(1:nrow(x), each = 2), 
                        s = x$subject,
                        seqnames = x$seqnames[1])
    box <- data.frame(boxYmin = 0, boxYmax = 1, boxXmin = 1, 
                      boxXmax = chromosomeLengths[[as.character(x$seqnames[1])]], 
                      seqnames = x$seqnames[1])
    list(lines = lines, box = box)
  })
  sites <- do.call(rbind, lapply(d, "[[", 1))
  boxes <- do.call(rbind, lapply(d, "[[", 2))
  
  ggplot() + 
    theme_bw() + 
    scale_color_manual(name = 'Dog', values = colorRampPalette(brewer.pal(12, "Paired"))(6)) +
    geom_line(data = sites, 
              alpha = alpha, 
              size = 0.25,
              aes(x, y, group = g, color = s)) +
    geom_rect(data = boxes, 
              color = "black", 
              size = 0.25,
              alpha = 0,
              mapping = aes(xmin = boxXmin, 
                            xmax = boxXmax, 
                            ymin = boxYmin, 
                            ymax = boxYmax)) + 
    facet_grid(seqnames ~ ., switch = "y") + 
    #scale_x_continuous(expand = c(0, 0)) + 
    labs(x = "Genomic position", y = "") + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), 
          strip.text.y = element_text(size = 12, angle = 180), 
          strip.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

genomeMap <- intSiteDistributionPlot(filter(intSites, ! subject %in% c('HO2', 'M12')), chromosomeLengths, alpha = 1)
ggsave(genomeMap, file = 'tables_and_figures/chromosome_map.pdf', height = 10, width = 8, units = 'in')










# Integrations within vectors and genomic exons.
#--------------------------------------------------------------------------------------------------

library(grid)
library(gridExtra)


createF8exonPlots <- function(sites){
  g <- subset(gt23::canFam3.refSeqGenesGRanges, name2 == 'F8')
  d <- tibble(exonStarts = strsplit(g$exonStarts, ","), exonEnds = strsplit(g$exonEnds, ","))
  d <- unnest(d, exonStarts = d$exonStarts, exonEnds = d$exonEnds)
  
  rev(lapply(1:nrow(d), function(n){
    x <- d[n,]
    x$exonStarts <- as.integer(x$exonStarts)
    x$exonEnds <- as.integer(x$exonEnds)
    o <- subset(sites, start >= x$exonStarts & start <= x$exonEnds)
    
    if(nrow(o) == 0) return(ggplot())
    
    ggplot(o, aes(start, 0, fill = strand)) + 
      theme_bw() +
      scale_fill_manual(values = c('dodgerblue2', 'gold2')) +
      geom_point(size = 3, shape = 21, position = position_jitter(width = 0, height = 0.1)) +
      scale_x_continuous(sec.axis = dup_axis(), expand = c(.1, .1)) +
      scale_y_continuous(sec.axis = dup_axis(), expand = c(.01, .01)) +
      labs(x = '', y = '') +
      theme(legend.position="none",
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
      )
  }))
}

createVectorExonPlots <- function(vectorSitesFile, exons){
  vectorSites <- new.env()
  load(vectorSitesFile, envir = vectorSites)
  vectorSites$sites <- subset(vectorSites$sites, seqnames != 'chrUn')
  
  vectorPlotData <- tibble(subject = sub('^p', '', vectorSites$sites$subject),
                           position = vectorSites$sites$position, 
                           strand = vectorSites$sites$strand)
  
  vectorPlot <- 
    ggplot(vectorPlotData, aes(position, 0, fill = strand)) + 
    theme_bw()+
    scale_fill_manual(values = c('dodgerblue2', 'gold2')) +
    geom_jitter(size = 3, height = 0.75, width = 0, shape = 21, color = 'black') +
    ylim(-3,3) +
    facet_grid(subject~., switch = "y") +
    theme(legend.position="none",
          strip.text.y = element_text(angle = 180),
          strip.background = element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  return(c(list(vectorPlot), lapply(unique(vectorPlotData$subject), function(x){
    gridExtra::arrangeGrob(grobs = createF8exonPlots(subset(intSites, subject == x))[exons],  ncol = 7)  
  })))
}

# Use grid::grid.draw() to view grob layouts.
#
# set.seed(42)
# grobs <- createVectorExonPlots('AAVengeR/outputs/singleChainVector/sites.RData', 11:17)
# single_vectorExon_plots <- arrangeGrob(grobs = grobs, 
#                                        layout_matrix = rbind(c(1,1,1,1,1,1,1),c(1,1,1,1,1,1,1), c(1,1,1,1,1,1,1),  c(1,1,1,1,1,1,1),
#                                                              c(NA,2,2,2,2,2,NA), c(NA,3,3,3,3,3,NA), c(NA,4,4,4,4,4,NA)))
# ggsave(single_vectorExon_plots, file = 'tables_and_figures/singleVector_genomicExonHits.pdf', height = 8, width = 11, units = 'in')
# 
# 
# 
# set.seed(42)
# light_vectorExon_plots <- arrangeGrob(grobs = createVectorExonPlots('AAVengeR/outputs/vectors_light/sites.RData', 11:17), 
#                                       layout_matrix = rbind(c(1,1,1,1,1,1,1),c(1,1,1,1,1,1,1), c(1,1,1,1,1,1,1),
#                                                             c(NA,2,2,2,2,2,NA), c(NA,3,3,3,3,3,NA), c(NA,4,4,4,4,4,NA)))
# ggsave(light_vectorExon_plots, file = 'tables_and_figures/lightVector_genomicExonHits.pdf', height = 6.7, width = 11.25, units = 'in')
# ggsave(light_vectorExon_plots, file = 'tables_and_figures/lightVector_genomicExonHits.svg',  height = 6.7, width = 11.25, units = 'in')
# 
# set.seed(42)
# heavy_vectorExon_plots <- arrangeGrob(grobs = createVectorExonPlots('AAVengeR/outputs/vectors_heavy/sites.RData', 11:17), 
#                                       layout_matrix = rbind(c(1,1,1,1,1,1,1),c(1,1,1,1,1,1,1), c(1,1,1,1,1,1,1),
#                                                             c(NA,2,2,2,2,2,NA), c(NA,3,3,3,3,3,NA), c(NA,4,4,4,4,4,NA)))
# ggsave(heavy_vectorExon_plots, file = 'tables_and_figures/heavyVector_genomicExonHits.pdf',  height = 6.7, width = 11.25, units = 'in')
# ggsave(heavy_vectorExon_plots, file = 'tables_and_figures/heavyVector_genomicExonHits.svg',  height = 6.7, width = 11.25, units = 'in')























# ITR plots.
#--------------------------------------------------------------------------------------------------
# (!) Need ITR names
# 
# ITRplot.shiftVal <- 32
# ITRplot.ITRline  <- data.frame(x1 = 0 , x2 = 98, y = 100)
# ITRplot.ITRticks <- data.frame(x = c(0, 26.5, 37.5, 48.5, 59.5, 69.5, 98), y1 = 98, y2 = 102)
# ITRplot.ITRtips  <- data.frame(x = c(37.5, 59.5), y1 = 0, y2 = 100)
# ITRplot.ITRlabels <- c("A", "C'/B", "C/B'","B/C'", "B'/C", "A'")
# ITRplot.ITRlabelsPos    <- c(14, 32, 44, 55, 65, 83)
# ITRplot.ITRlabelsHeight <- 105
# 
# ITRplotSetup <- 
#   list(theme_bw(),
#        scale_x_continuous(),
#        theme(panel.grid.major.x = element_blank(),
#              panel.grid.minor.x = element_blank(),
#              axis.text.x = element_text(angle = 90, hjust = 1)),
#         geom_segment(data = ITRplot.ITRline,  aes(x = x1, y = y,  xend = x2, yend = y),  inherit.aes = FALSE),
#         geom_segment(data = ITRplot.ITRticks, aes(x = x,  y = y1, xend = x,  yend = y2), inherit.aes = FALSE),
#         geom_segment(data = ITRplot.ITRtips,  aes(x = x,  y = y1, xend = x,  yend = y2), inherit.aes = FALSE, linetype = 'dotted'),
#         annotate("text", x = ITRplot.ITRlabelsPos, y = ITRplot.ITRlabelsHeight, label = ITRplot.ITRlabels),
#         labs(x = 'ITR NTs', y = 'Integrations'))
# 
# ITR.direction.plot <-
#   mutate(expIntSites.noF8, ITRnts = nchar(ltrRepSeq) - ITRplot.shiftVal) %>%
#   select(subject, LTRseqName, estAbund, ITRnts) %>%
#   group_by(LTRseqName, ITRnts) %>%
#   summarise(nSites = n()) %>%
#   ungroup() %>%
#   mutate(LTRseqName = factor(as.character(LTRseqName), levels = c('u', '5p', '3p'))) %>%
#   ggplot(aes(ITRnts, nSites, fill = LTRseqName)) + 
#     ITRplotSetup +
#     geom_col(color = 'black', size = 0.25) +
#     scale_fill_manual(name = 'ITR', values = c('goldenrod2', 'dodgerblue', 'mediumseagreen')) + 
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
#         legend.title=element_text(size = 14), legend.text=element_text(size = 12))
# 
# ggsave(ITR.direction.plot, file = 'tables_and_figures/ITR.direction.plot.pdf')
# 
# ann_text <- data.frame(label      = c("A", "B", "B'", "C'", "C", "A'",    "A", "C'", "C", "B", "B'", "A'",  "A"),
#                        LTRseqName = factor(c("3' direction", "3' direction", "3' direction", "3' direction", "3' direction", "3' direction",  
#                                              "5' direction", "5' direction", "5' direction", "5' direction", "5' direction", "5' direction", "Undetermined"), levels = c("3' direction","5' direction","Undetermined")),
#                        x = c(14, 32, 44, 53, 63, 81, 14, 32, 44, 53, 63, 81, 14), y = 105)
# 
# ITRplot.ITRline$y <- 53
# ITRplot.ITRticks$y1 <- 51
# ITRplot.ITRticks$y2 <- 55
# ITRplot.ITRtips$y2 <- 51
# ann_text$y <- 58
# 
# report$ITR.abundanceDirection.plot <- 
#   mutate(expIntSites.noF8, ITRnts = nchar(LTRseqsRep) - ITRplot.shiftVal) %>%
#   mutate(abundBin = cut(estAbund, breaks = c(0, 5, 10, 50, max(expIntSites.noF8$estAbund)))) %>%
#   select(subject, abundBin, LTRseqName, estAbund, ITRnts) %>%
#   group_by(LTRseqName, abundBin, ITRnts) %>%
#   summarise(nSites = n()) %>%
#   ungroup() %>%
#   mutate(abundBin = recode(abundBin, '(0,5]' = '1 - 5', '(5,10]' = '6 - 10', '(10,50]' = '11 - 50', '(50,131]' = ' 51 - 131')) %>%
#   mutate(LTRseqName = recode(LTRseqName, u = "Undetermined", "5p" = "5' direction", "3p" = "3' direction")) %>%
#   mutate(LTRseqName = factor(as.character(LTRseqName), levels = c('Undetermined', "5' direction", "3' direction"))) %>%
#   ggplot(aes(ITRnts, nSites, fill = abundBin)) +
#   theme_bw() +
#   geom_col(color = 'black', size = 0.25) +
#   scale_fill_manual(name = 'Abundance', values = c('goldenrod2', 'dodgerblue', 'mediumseagreen', 'red')) + 
#   facet_grid(LTRseqName~.) +
#   geom_segment(data = ITRplot.ITRline, aes(x = x1, y = y, xend = x2, yend = y), inherit.aes = FALSE) +
#   geom_segment(data = ITRplot.ITRticks, aes(x = x, y = y1, xend = x, yend = y2), inherit.aes = FALSE) +
#   geom_segment(data = ITRplot.ITRtips, aes(x = x, xend = x, y = y1, yend = y2), inherit.aes = FALSE, linetype = 'dotted') +
#   geom_text(data = ann_text, aes(label = label, x = x, y = y), hjust = 0, inherit.aes = FALSE) +
#   labs(x = 'ITR position', y = 'Number of Integrations') +
#   ylim(c(0, 60)) +
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=14), legend.title =element_text(size=14),
#         strip.text.y = element_text(size = 12, angle = 0),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave(report$ITR.abundanceDirection.plot, file = 'tables_and_figures/ITR.abundanceDirection.plot.pdf')
# 
# 
# report$ITR.abundance.plot <- 
#   mutate(expIntSites.noF8, ITRnts = nchar(LTRseqsRep) - ITRplot.shiftVal) %>%
#   mutate(abundBin = cut(estAbund, breaks = c(0, 5, 10, 50, max(expIntSites.noF8$estAbund)))) %>%
#   select(abundBin, ITRnts) %>%
#   group_by(abundBin, ITRnts) %>%
#   summarise(nSites = n()) %>%
#   ungroup() %>%
#   mutate(abundBin = recode(abundBin, '(0,5]' = '1 - 5', '(5,10]' = '6 - 10', '(10,50]' = '11 - 50', '(50,130]' = ' 51 - 130')) %>%
#   ggplot(aes(ITRnts, nSites, fill = abundBin)) + 
#     ITRplotSetup +
#     geom_col(color = 'black', size = 0.25) +
#     scale_fill_manual(name = 'Abundance', values = c('goldenrod2', 'dodgerblue', 'mediumseagreen', 'red')) + 
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
#         legend.title=element_text(size = 14), legend.text=element_text(size = 12))
# 
# ggsave(report$ITR.abundance.plot, file = 'tables_and_figures/ITR.abundance.plot.pdf')






  






# Local homology tests
#--------------------------------------------------------------------------------------------------


# Create blastn db for single chain vector.
system(paste0(blastBin, '/makeblastdb -in data/sequenceData/singleChainVector.ITR.fasta -dbtype nucl -out singleChainVectorITR.db'), ignore.stderr = TRUE)


# Indentify random sites (5x number of exp sites) which can be flanked +/- 1500 NT without encountering the ends of chromosomes or sequences with Ns.
if(! file.exists('data/localHomology_randomSites.rds')){
  o <- randomSites_annotated[1:(nrow(expIntSites.noF8)*50),]
  o$s <- ntile(1:nrow(o), 40)
  
  cluster <- makeCluster(40)
  localHomology_randomSites <- 
    bind_rows(parLapply(cluster, split(o, o$s), function(y){
      library(BSgenome.Cfamiliaris.UCSC.canFam3)
      library(Biostrings)
      library(dplyr)
    
      bind_rows(lapply(split(y, 1:nrow(y)), function(x){
        if(x$seqnames == 'chrM') return(tibble())
        seq <- tryCatch({
                        as.character(DNAStringSet(getSeq(BSgenome.Cfamiliaris.UCSC.canFam3, x$seqnames, x$start - 2000, x$start + 2000)))
                        }, 
                        error = function(cond) {
                         return('N')
                        })
       if(grepl('N', seq)) return(tribble())
       x
     }))
    }))

  saveRDS(localHomology_randomSites, file = 'data/localHomology_randomSites.rds')
} else {
  localHomology_randomSites <- readRDS('data/localHomology_randomSites.rds')
}

createLocalHomologyDataSet <- function(expSites, randomSitesPool, 
                                       radii = c(25, 50, 100, 250, 500, 1000, 1500, 2000), 
                                       randomSets = 10){
    usedPosids <- vector()
  
    bind_rows(lapply(1:randomSets, function(seed){
      set.seed(seed)
      
      randomSitesPool2 <- subset(randomSitesPool, ! posid %in% usedPosids)
      
      randomSites <- sample_n(randomSitesPool2, nrow(expSites))
      usedPosids <<- c(usedPosids, randomSites$posid)
      
      bind_rows(lapply(radii, function(radius){
  
        message('Starting random set ', seed, ' with radius ', radius)
        
        getGenomeSeq <- function(seqnames, pos, radius){
          tryCatch({
            as.character(getSeq(BSgenome.Cfamiliaris.UCSC.canFam3, seqnames, (pos - radius), (pos + radius)))
          },
          error = function(cond) {
            return(NA)
          })
        }
        
        cluster <- makeCluster(40)
        clusterExport(cluster, varlist = c('getGenomeSeq', 'tmpFile', 'blastBin'), envir = environment())
        clusterEvalQ(cluster, library("dplyr"))
        clusterEvalQ(cluster, library("Biostrings"))
        clusterEvalQ(cluster, library("BSgenome.Cfamiliaris.UCSC.canFam3"))
        
        expSites$s <- ntile(1:nrow(expSites), 40)
        randomSites$s <- ntile(1:nrow(randomSites), 40)
        
        expSites <- bind_rows(parLapply(cluster, split(expSites, expSites$s), function(x){
                       rowwise(x) %>% mutate(seq = getGenomeSeq(seqnames, start, radius)) %>% ungroup()
                    }))
        
         randomSites <- bind_rows(parLapply(cluster, split(randomSites, randomSites$s), function(x){
                          rowwise(x) %>% mutate(seq = getGenomeSeq(seqnames, start, radius)) %>% ungroup()
                        }))
         
       
         i <- which(is.na(expSites$seq) | is.na(randomSites$seq))
         if(length(i) > 0){
           expSites <- expSites[-i,]
           randomSites <- randomSites[-i,]
         }
         
         expSitesSet <- DNAStringSet(expSites$seq)
         names(expSitesSet) <- expSites$posid
         
         randomSitesSet <- DNAStringSet(randomSites$seq)
         names(randomSitesSet) <- randomSites$posid
         
         runBlast <- function(sites){
           bind_rows(parLapply(cluster, as.list(split(sites, ntile(1:length(sites), 40))), function(x){
              f <- tmpFile()
              writeXStringSet(x, file = f)
              
              system(paste0(blastBin, '/blastn -reward 1 -penalty -1 -gapopen 0 -gapextend 2 -evalue 500 -word_size 4 -dust no -soft_masking false -outfmt 6 -query ', 
                            f, ' -db singleChainVectorITR.db -num_threads 10 -out ', 
                            paste0(f, '.blast')))
              blast <- read.table(paste0(f, '.blast'), sep = '\t', header = FALSE)
              invisible(file.remove(list.files(pattern = f)))
              names(blast) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
              group_by(blast, qname) %>% arrange(desc(bitscore)) %>% dplyr::slice(1) %>% ungroup()
            }))
          }
         
         expSitesBlast <- runBlast(expSitesSet)
         randomSitesBlast <- runBlast(randomSitesSet)
         
         stopCluster(cluster)
         
         expSitesBlast <- left_join(expSitesBlast, dplyr::select(expSites, posid, seq), by = c('qname' = 'posid'))
         randomSitesBlast <- left_join(randomSitesBlast, dplyr::select(randomSites, posid, seq), by = c('qname' = 'posid'))
        
         createMissingBlastResults <- function(qnames){
           tibble(qname = qnames,
                  sseqid = 'ITR',
                  pident = 0, length = 0, mismatch = 0, gapopen = 0, qstart = 0,
                  qend = 0, sstart = 0,  send = 0,  evalue = 100, bitscore = 0, seq = NA)
         }
         
         expSitesMissing <- expSites$posid[! expSites$posid %in% expSitesBlast$qname]
         randomSitesMissing <- randomSites$posid[! randomSites$posid %in% randomSitesBlast$qname]
         
         if(length(expSitesMissing) > 0){
           message('  adding ', length(expSitesMissing), ' missing sites to experimental blast results.')
           expSitesBlast <- bind_rows(expSitesBlast, createMissingBlastResults(expSitesMissing))
         }
         
         if(length(randomSitesMissing) > 0){
           message('  adding ', length(randomSitesMissing), ' missing sites to random blast results.')
           randomSitesBlast <- bind_rows(randomSitesBlast, createMissingBlastResults(randomSitesMissing))
         }
         
         tibble(randomSeed = seed,
                siteRadiusDist = radius,
                expSiteBitScores = list(scores = expSitesBlast$bitscore),
                expSiteAlnLengths = list(lengths = expSitesBlast$length),
                expSiteAlnSeqs = list(seqs = expSitesBlast$seq),
                expSiteBitScoreMean = sprintf("%.2f", mean(expSitesBlast$bitscore)),
                expSiteAlignmentLengthMean = sprintf("%.2f", mean(expSitesBlast$length)),
                randSiteBitScores = list(scores = randomSitesBlast$bitscore),
                randSiteAlnLengths = list(lengths = randomSitesBlast$length),
                randSiteAlnSeqs = list(seqs = randomSitesBlast$seq),
                randSiteBitScoreMean = sprintf("%.2f", mean(randomSitesBlast$bitscore)),
                randSiteAlignmentLengthMean = sprintf("%.2f", mean(randomSitesBlast$length)),
                bitScoreWilcoxPval = wilcox.test(expSitesBlast$bitscore, randomSitesBlast$bitscore, alternative = "greater")$p.value,
                alignmentLengthWilcoxPval = wilcox.test(expSitesBlast$length, randomSitesBlast$length, alternative = "greater")$p.value)
       }))
  }))
}


if(! file.exists('data/localHomology.noF8.rds')){
  localHomology_randomSites <- subset(localHomology_randomSites, nearestFeature != 'F8') 
  d <- expIntSites.noF8[! duplicated(expIntSites.noF8$posid),]
  localHomology.noF8 <- createLocalHomologyDataSet(d, localHomology_randomSites)
  saveRDS(localHomology.noF8, file = 'data/localHomology.noF8.rds')
} else {
  localHomology.noF8 <- readRDS('data/localHomology.noF8.rds')
}

set.seed(1)
localHomology_bitScoresVsRadius <- 
  ggplot(subset(localHomology.noF8, bitScoreWilcoxPval < 0.01), aes(factor(siteRadiusDist), bitScoreWilcoxPval , color = factor(randomSeed))) + 
  theme_bw() +
  scale_color_manual(name = 'Random\nsites set', values = createColorPalette(10)) +
  geom_hline(yintercept = 1e-10, color = 'black') +
  geom_hline(yintercept = 1e-3, color = 'red') +
  geom_jitter(height = 0, width = 0.3, size = 2.5) +
  theme(text = element_text(size=15), panel.grid.major.x = element_blank(),  
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) +
  scale_y_continuous() +
  labs(x = 'Site radius', y = 'Alignment bitscore Wilcox p-value') 
  
ggsave(localHomology_bitScoresVsRadius, file = 'tables_and_figures/localHomology_bitScoresVsRadius.pdf', units = 'in', height = 5)


o <- subset(localHomology.noF8, siteRadiusDist == 100)

d <- data.frame(table(unlist(o$expSiteBitScores))) %>% 
     mutate(sites = 'Experimental') %>%
     bind_rows(data.frame(table(unlist(o$randSiteBitScores))) %>% mutate(sites = 'Random')) %>% 
     mutate(Var1 = factor(round(as.numeric(Var1))))

localHomology_randomVsExpAlignmentScores <- 
  ggplot(d, aes(Var1, Freq, fill = sites)) + 
  theme_bw()+
  scale_fill_manual(name = 'Sites', values = c('black', 'gray75')) +
  geom_col(position = 'dodge') + 
  scale_y_continuous(labels = scales::comma) +
  labs(x = 'Alignment bitscore', y = 'Sites') +
  theme(text = element_text(size=15), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) 

ggsave(localHomology_randomVsExpAlignmentScores, 
       file = 'tables_and_figures/localHomology_randomVsExpAlignmentScores.pdf', units = 'in', height = 5)





# Heatmaps
#--------------------------------------------------------------------------------------------------

expIntSites.noF8$refGenome <- 'canFam3'
write(c('sampleName,GTSP,patient', paste0(unique(expIntSites.noF8$subject), ',', unique(expIntSites.noF8$subject), ',', 'xxx', collapse = '\n')), file = 'samples')
write('seqnames,strand,position,sampleName,refGenome', file = 'sites')

set.seed(42)
# /home/everett/manuscript.geneTherapy.HSCdiversity.new/data/tmp
write.table(dplyr::select(subset(expIntSites.noF8, seqnames %in% c(paste0('chr', 1:38), 'chrX', 'chrY')), seqnames, strand, start, subject, refGenome), 
            file = 'sites', col.names = FALSE, row.names = FALSE, sep = ',', append = TRUE, quote = FALSE)


comm <- paste(Rscript_path, './genomicHeatMapMaker/genomic_heatmap_from_file.R samples ', 
              '-c ./genomicHeatMapMaker/INSPIIRED.yml ',
              '-o data/genomicHeatMapMakerOutput ',
              '-f sites ',
              '-r canFam3')

if(! dir.exists('data/genomicHeatMapMakerOutput')) system(comm)

report$genomicHeatmap <- within(
  list(), {
    heatmap_sample_info <- read.csv('samples')
    gen_heatmap <- readRDS('data/genomicHeatMapMakerOutput/roc.res.rds')
    
    heatmap_scale <- seq(0.2, 0.8, 0.1)
    gen_heatmap_colors <- colorspace::diverge_hsv(
      length(heatmap_scale), h = c(240, 0), v = 1, power = 1)
    
    select_gen_features <- row.names(gen_heatmap$ROC)
    select_gen_features <- c(
      "boundary.dist", "start.dist", "general.width", "gene.width",
      "within_refSeq_gene", "refSeq_counts.10k", "refSeq_counts.100k",
      "refSeq_counts.1M", "GC.100", "GC.1k", "GC.10k", "GC.100k", "GC.1M",
      "CpG_counts.1k", "CpG_counts.10k", "CpG_density.10k", "CpG_density.100k",
      "CpG_density.1M")
    
    gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
    
    heatmap_sample_levels <- c("H19", "J60", "Linus", "M06", "M50", "M66")
    
    heatmap_figure_labels <- heatmap_sample_levels
    
    stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
    gen_comp_stats <- structure(cut(
      gen_heatmap$pvalues$op[select_gen_features, 1],
      stat_cuts,
      labels = c("***", " **", " * ", "   "),
      include.lowest = TRUE),
      names = select_gen_features)
    gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
    
    plot_data <- gen_heatmap$ROC %>%
      reshape2::melt() %>%
      mutate(
        feat = Var1,
        comp.sym = gen_comp_stats[Var1],
        Var1 = paste0(Var1, " - ", comp.sym),
        Var1 = factor(Var1, levels = gen_row_names),
        Var2 = factor(Var2, levels = heatmap_sample_levels),
        grp = " ",
        sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
        sym = cut(
          sig, stat_cuts, labels = c("***", " **", " * ", "   "),
          include.lowest = TRUE))
    
    levels(plot_data$Var2) <- heatmap_figure_labels
    
    plot_data$Var1 <- gsub('\\*', '', as.character(plot_data$Var1))
    plot_data$Var1 <- gsub('\\-', '', as.character(plot_data$Var1))
    
    
    levels_var1= gsub('\\*', '',gen_row_names)
    levels_var1= gsub('\\-', '',levels_var1)
    plot_data$Var1=factor(plot_data$Var1,levels=levels_var1)
    
    plot_data$sym  <- ''
    
    gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(
        breaks = c(0.2, 0.4, 0.6, 0.8),
        low = gen_heatmap_colors[1],
        mid = gen_heatmap_colors[round(length(heatmap_scale)/2)],
        high = gen_heatmap_colors[length(heatmap_scale)],
        midpoint = 0.5) +
      guides(fill = guide_colorbar(
        title.position = "left", title.hjust = 0.5,
        direction = "horizontal")) +
      labs(x = NULL, y = NULL, fill = "ROC\nScore") +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
        axis.text.x.top = element_text(
          angle = 90, hjust = 0.5, vjust = 0.5, size = 12),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
  })


ggsave(report$genomicHeatmap$gen_plot, file = 'tables_and_figures/genomicHeatMap.pdf')
file.remove(c('sites_mrcs.gen', 'samples', 'sites'))







# Additional NTs
#--------------------------------------------------------------------------------------------------


parse_cdhitest_output <- function (file) {
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, ">Cluster"))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x) {
    gsub("\\.\\.\\.", "", unlist(lapply(str_match_all(x,
                                                      ">([^\\s]+)"), function(y) {
                                                        y[, 2]
                                                      })))
  })
}


representativeSeq <- function(s){
  if(length(s) == 1 | n_distinct(s) == 1) return(s[1])
  m <- as.matrix(stringdist::stringdistmatrix(s))
  d <- apply(m, 1, sum) 
  s[which(d == min(d))[1]]
}

# Determine the sequences and numbers of additional NTs.
intSites <- rowwise(intSites) %>%
            mutate(LTRseqAdditionNTs = substr(LTRseq2sRep, nchar(LTRseqsRep) + 1, nchar(LTRseq2sRep)),
                   nLTRseqAdditionNTs = nchar(LTRseqAdditionNTs)) %>%
            ungroup()


# Find unique sites allowing for the same site to be seen betwen subjects.
intSites$s <- paste(intSites$posid, intSites$subject, intSites$LTRseq2sRep)
d <- intSites[! duplicated(intSites$s),]

# Create an within / outside F8 flag.
d1 <- data.frame(table(subset(d, nearestFeature == 'F8')$nLTRseqAdditionNTs))
d1$g <- 'Within Factor VIII'
d2 <- data.frame(table(subset(d, nearestFeature != 'F8')$nLTRseqAdditionNTs))
d2$g <- 'Outside Factor VIII'

tab <- bind_rows(d1, d2)
tab$Var1 <- as.integer(as.character(tab$Var1))

# Count the number of integrations with >= 25 additional NTs.
e1 <- sum(subset(tab, Var1 >= 25 & g == 'Within Factor VIII')$Freq)
e2 <- sum(subset(tab, Var1 >= 25 & g == 'Outside Factor VIII')$Freq)

tab <- subset(tab, Var1 < 25) 
tab$Var1 <- as.character(tab$Var1)
tab <- bind_rows(tab, tibble(Var1 = '25+', Freq = e1, g = 'Within Factor VIII'), tibble(Var1 = '25+', Freq = e2, g = 'Outside Factor VIII'))
tab$Var1 <- factor(tab$Var1, levels = c(as.character(0:24), '25+'))
tab$g <- factor(tab$g, levels = c('Within Factor VIII', 'Outside Factor VIII'))

report$additionaNTplot <- 
  ggplot(tab, aes(Var1, Freq, fill = g)) +
  theme_bw() +
  scale_fill_manual(name = '', values = c('gray25', 'gray75')) +
  geom_col(position='dodge') +
  scale_y_continuous(breaks = seq(0, max(tab$Freq), by = 250), label = scales::comma) +
  labs(x = 'Number of additional NTs', y = 'Integrations') +
  theme(panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(report$additionaNTplot, file = 'tables_and_figures/additionalNTdistributions.pdf', height = 5, width = 10, units = 'in')


# Second version with different binning.
tab <- bind_rows(d1, d2)
tab$Var1 <- as.integer(as.character(tab$Var1))

# Count the number of integrations with >= 15 additional NTs.
e1 <- sum(subset(tab, Var1 >= 15 & g == 'Within Factor VIII')$Freq)
e2 <- sum(subset(tab, Var1 >= 15 & g == 'Outside Factor VIII')$Freq)

tab <- subset(tab, Var1 < 15) 
tab$Var1 <- as.character(tab$Var1)
tab <- bind_rows(tab, tibble(Var1 = '15+', Freq = e1, g = 'Within Factor VIII'), tibble(Var1 = '15+', Freq = e2, g = 'Outside Factor VIII'))
tab$Var1 <- factor(tab$Var1, levels = c(as.character(0:14), '15+'))
tab$g <- factor(tab$g, levels = c('Within Factor VIII', 'Outside Factor VIII'))

report$additionaNTplot2 <- 
  ggplot(tab, aes(Var1, Freq, fill = g)) +
  theme_bw() +
  scale_fill_manual(name = '', values = c('gray25', 'gray75')) +
  geom_col(position='dodge') +
  scale_y_continuous(breaks = seq(0, max(tab$Freq), by = 250), label = scales::comma) +
  labs(x = 'Number of additional NTs', y = 'Integrations') +
  theme(panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(report$additionaNTplot2, file = 'tables_and_figures/additionalNTdistributions2.pdf', height = 5, width = 10, units = 'in')






extraNT_Table <- function(d){
  d1 <- tibble(seq = d$LTRseqAdditionNTs)
  d1$s <- paste0('s', 1:nrow(d1))
  write(paste0('>', d1$s, '\n', d1$seq), file = 'extraNTs')
  system('/home/everett/software/cdhit/cd-hit-est -i extraNTs -o extraNTs.out -c 0.95  -aS 0.70 -aL 0.70')
  
  clusters <- parse_cdhitest_output('extraNTs.out.clstr')
  clusters <- clusters[order(unlist(lapply(clusters, length)), decreasing = TRUE)]
  file.remove(c('extraNTs', 'extraNTs.out', 'extraNTs.out.clstr'))
  
  bind_rows(lapply(clusters, function(x){
    o <- subset(d1, s %in% x)
    r <- representativeSeq(o$seq)
    e <- nchar(r) - 25
    r2 <- r
    if(e > 0) r2 <- paste0(substr(r, 1, 15), "... +", e)
    tibble(n = nrow(o), r2 = r2, r = r, l = list(o$seq))
  }))
}

extraNT_Table1 <- extraNT_Table(subset(d, nearestFeature == 'F8' & nLTRseqAdditionNTs >= 10))
extraNT_Table1 <- extraNT_Table1[1:15,]
extraNT_Table1$anno <- c("Region of ITR before first dumbbell tip", "Vector poly-A signal",  "End of promotor, start of F8",
                         "Region of ITR before first dumbbell tip", "F8", "F8", "Vector (between promotor and F8)",
                         "Vector poly-A signal", "Region of ITR after second dumbbell tip", "F8", "F8", "F8", "F8", "F8", "F8")



extraNT_Table2 <- extraNT_Table(subset(d, nearestFeature != 'F8' & nLTRseqAdditionNTs >= 10))
extraNT_Table2 <- extraNT_Table2[1:15,]
extraNT_Table2$anno <- c("Internal edge of 5' ITR + promotor", "Internal edge of 5' ITR", "HC/LC vector between 5' ITR and promotor",
                         "HC/LC - end of vector poly-A signal + start of 3' ITR", "Full - end of vector poly-A signal + start of 3' ITR",
                         "Region of ITR before first dumbbell tip", "Internal edge of ITRs", "Unknown", "Region of ITR before first dumbbell tip",
                         "Region of ITR before first dumbbell tip", "Internal edge of 5' ITR + promotor", "chr25:34548608-34548676",                
                         "Region of ITR before first dumbbell tip", "Region of ITR before first dumbbell tip", "chrX:72159400-72159435")

save.image('data/report.RData')
