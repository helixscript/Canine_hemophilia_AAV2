library(IRanges)
library(dplyr)

s <- IRanges(start = c(10, 50, 100), end=c(35, 90, 150), names = c('A', 'B', 'C'))
q <- IRanges(start = c(5, 30, 95, 110, 120), end=c(140, 60, 120, 130, 140))

r <- bind_rows(lapply(split(q, 1:length(q)), function(x){
  o <- IRanges::findOverlaps(x, s)
  d <- data.frame(IRanges::intersect(x, s))
  d$strand <- '*'
  d$feature <- names(s[subjectHits(o)])
  d
}))


df <- data.frame(start = start(v), end = end(v))
ggplot() + geom_segment(aes(x = start, y = 1, xend = end, yend = 1), data = df)