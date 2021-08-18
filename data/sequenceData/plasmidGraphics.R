library(tidyverse)
options(stringsAsFactors = FALSE)

d <- read.table('plasmidsFigureFeatures.csv', sep = ',', header = TRUE, quote = '', comment.char = '#')
d$vector <- factor(d$vector)
d$y = as.integer(d$vector)

features <- c("5' ITR", "Spacer", "Promoter", "Intron", "Secretion signal", "B-domain", "F8", "F8 LC", "F8 HC", "PolyA", "3' ITR")
colors   <- c('firebrick1', 'gray80', '#ffc000', 'purple', 'pink', 'orange', '#00af50', '#c5dfb4', '#538235', '#00afef', 'red3') 

d$feature <- factor(d$feature, levels = features)

p <- ggplot(d) + 
     theme_bw() +
     scale_fill_manual(values = colors) +
     geom_rect(mapping=aes(xmin=start, xmax=end, ymin=y-0.25, ymax=y+0.25, fill=feature), color = 'black', size = 1) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(p, file = 'plasmidGraphics.pdf', units = 'in', width = 20)