library(BSgenome.Athaliana10.TAIR.TAIR10)
library(rGADEM)
library(hbmcbioR)
peaks <- ChIPseeker::readPeakFile("/wrk/luolin/work/Luolin1_4__araDNase/pre_processed/peaks/A6_6_Ara51216m.bed")
rang <- ranges(peaks)
bed <- data.frame(chr=as.factor(seqnames(peaks)),start=as.numeric(rang@start),end=as.numeric(rang@start + rang@width -1))
bed <- bed[bed$chr!="ChrC",]
bed <- bed[bed$chr!="ChrM",]
seqbio <- Biostrings::DNAStringSet()
for(a in 1:length(levels(seqnames(peaks)))){
  chr <- bed %>% filter(chr==levels(seqnames(peaks))[[a]])
  seqbio <- c(seqbio,Biostrings::DNAStringSet(BSgenome.Athaliana10.TAIR.TAIR10[[levels(seqnames(peaks))[[a]]]], start=chr$start, end=chr$end))
}
gadem <- rGADEM::GADEM(seqbio,verbose=1,genome="/wrk/data/genome/yz_genome_data/aragenome/TAIR10_chr_all.fas")
pwmli <- rGADEM::getPWM(gadem)
pwmlist_t <- TFBSTools::PWMatrixList()
for(i in 1:length(pwmli)){
  t_tmp <- universalmotif::convert_motifs(pwmli[[i]],"TFBSTools-PWMatrix")
  t_tmp@strand <- "+-"
  pwmlist_t[[i]] <- t_tmp %>% hbmcbioR::motif_get_name()
}
write_rds(pwmlist_t,"pwmlist_a6")
ggsave("motif_total_a6.png",view_motifs(pwmlist_t))
