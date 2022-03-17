library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(reshape2)
dele_nu <- function(df){
  for (a in 1:nrow(df)){
    if(a==nrow(df)){
      return(df)
    }
    if(df$start[[a]] < 0){
      df=df[-a,]
    }
    if(a==nrow(df)){
      return(df)
    }
    if(df$chr[[a]]=="Chr1" && df$end[[a]]>30427671){
      df=df[-a,]
    }
    if(a==nrow(df)){
      return(df)
    }
    if(df$chr[[a]]=="Chr2" && df$end[[a]]>19698289){
      df=df[-a,]
    }
    if(df$chr[[a]]=="Chr3" && df$end[[a]]>23459830){
      df=df[-a,]
    }
    if(a==nrow(df)){
      return(df)
    }
    if(df$chr[[a]]=="Chr4" && df$end[[a]]>18585056){
      df=df[-a,]
    }
    if(a==nrow(df)){
      return(df)
    }
    if(df$chr[[a]]=="Chr5" && df$end[[a]]>26975502){
      df=df[-a,]
    }
  }
}

##Download Arabidopsis DHS from DHS database
##5000 bp upstream and downstream of DHS center were taken
dhs <- read.delim("/wrk/yuanzhen/project_bs/openChr_pre/TAIR10_DHSs.gff", header=FALSE)
dhs_score <- read.csv("/wrk/yuanzhen/project_bs//openChr_pre/Download_Ath_DHSs_mean_std_max.csv")
#dhs_bw <- rtracklayer::import("~/openChr_pre/Ath_leaf_DNase (3).bw")
#x_tmp <- str_split_fixed(dhs_score$location,":",n=2) %>% as.data.frame()
#y_tmp <- str_split_fixed(x_tmp$V2,"-",2) %>% as.data.frame()
##identical(dhs$V4 %>% as.numeric(),y_tmp$V1 %>% as.numeric())
dhs_t <- cbind(dhs,dhs_score)
dhs_5000 <- dhs_t %>% arrange(desc(Leaf.mean)) %>% dplyr::slice(1:5000)
dhs_5000 <- dhs_5000 %>% arrange(desc(V5-V4))
dhs_10kb <- data.frame(start=(round((dhs_5000$V4 + dhs_5000$V5) / 2 - 5000)))
dhs_10kb$end <- round((dhs_5000$V4 + dhs_5000$V5) / 2 + 5000)
dhs_10kb$chr <- dhs_5000$V1
dhs_10kb <- dele_nu(dhs_10kb)
dhs_gr <- GenomicRanges::GRanges(seqnames=dhs_10kb$chr,ranges = IRanges::IRanges(start=dhs_10kb$start,end = dhs_10kb$end))
dhs_gr[which(GenomicRanges::width(dhs_gr) < 10001)] <- NULL

##get the top 1000 kmer and 
##Subtract the kmer value randomly disrupted by the original sequence,
##and finally divide it by the total base number of the original sequence
getkmer <- function(seqfile="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/A9_9_IAA__leaf__conc_0.1nM__time_2h__30n.gz",n=1000,kmer_num=9L){
  seqs=fjComm::getSeq_fqfachrFile(seqfile)
  kmer_cnt=fjComm::kmerCntBit(strings =seqs, k = kmer_num, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 10)
  seqs_sf <- stringi::stri_rand_shuffle(seqs)
  kmer_cnt_sf=fjComm::kmerCntBit(strings =seqs_sf, k = kmer_num, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 10) 
  kmer_total <- dplyr::left_join(kmer_cnt,kmer_cnt_sf,by="kmer") %>% filter(!counts.x < counts.y & !counts.x == counts.y)
  kmer_total$kmer_nor <- -log2((kmer_total$counts.x - kmer_total$counts.y) / sum(str_count(seqs)))
  kmer_top1000 <- kmer_total %>% dplyr::arrange(kmer_nor) %>% dplyr::slice(1:n)
  Biostrings::DNAStringSet(kmer_top1000$kmer)
}


##The 10kb regions were divided to 50bp bins, 
## and each bin was then assigned with the mean score of all ten-mers inside the bin
mer_count <- function(genome="/wrk/data/genome/yz_genome_data/aragenome/TAIR10_chr_all.fas",cnt_ob=kmer_stri,dhs_gr=dhs_gr,bins=50,windows=200,by="kmers"){
  genom_ara <- Biostrings::readDNAStringSet(genome)
  dhs_seq <- genom_ara[dhs_gr]
  total_count <- BiocParallel::bplapply(1:length(dhs_seq),function(x){
    seq_50 <- lapply(1:windows,function(y){Biostrings::subseq(dhs_seq[[x]],start = (y-1)*bins+1,end = y*bins)})
    if(by=="kmers"){
      count_m <- lapply(1:windows,function(z){Biostrings::countPDict(cnt_ob,seq_50[[z]])})
    }else if(by=="pwmlist"){
      count_m <- lapply(1:windows,function(z){Biostrings::countPWM(cnt_ob,seq_50[[z]],min.score = "95%")})
    }else if(by=="base"){
      count_m <- lapply(1:windows,function(z){Biostrings::countPattern(cnt_ob,seq_50[[z]])})
    }else{
      stop("please choose correctly form(kmers or pwmlist")
    }
    colSums(do.call("cbind",count_m))
  },BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE)) 
  do.call("rbind",total_count)
}
#melted_cormat <- reshape2::melt(xm)
# kmer_stri <- getkmer()
# kmer_count <- mer_count(cnt_ob=kmer_stri,dhs_gr=dhs_gr)

# write_rds(kmer_stri,"/wrk/yuanzhen/project_bs/openChr_pre/total_kmer_stri.Rds")
# write_rds(kmer_count,"/wrk/yuanzhen/project_bs/openChr_pre/total_kmer_count.Rds")
kmer_stri <- read_rds("/wrk/yuanzhen/project_bs/openChr_pre/total_kmer_stri.Rds")
kmer_count <- read_rds("/wrk/yuanzhen/project_bs/openChr_pre/total_kmer_count.Rds")
kmer_df <- data.frame(num=1:200,count=colSums(kmer_count))

png(paste0("kmerall.png"))
h_line = ComplexHeatmap::HeatmapAnnotation(count=ComplexHeatmap::anno_lines(kmer_df$count,height = unit(2, "cm")))
kmer_pt <- ComplexHeatmap::Heatmap(kmer_count,name="bykmer",
                                   column_title = "Genomic postion(-5kb,5kb)",
                                   column_title_side = "bottom",
                                   row_title = "DHSs in arabidophsis leaf cells",
                                   row_title_gp = gpar(fontsize = 14), #fontface = "bold"),
                                   clustering_distance_rows = function(m) dist(m,method = "maximum"),
                                   use_raster = T,cluster_columns = F,cluster_rows = F,
                                   #raster_by_magick = TRUE, raster_resize_mat = max,
                                   col = circlize::colorRamp2(c(-1, 0, 1), c("red", "white", "blue")),
                                   width = unit(6, "cm"), height = unit(9, "cm"),
                                   top_annotation = h_line)

ComplexHeatmap::draw(kmer_pt)

dev.off()

pwmlist_leaf <- readr::read_rds("/var/www/html/sysbiomotif/yuanzhen/arabidopsis/leaf/motiftotal.Rds")
lapply(1:length(pwmlist_leaf),function(x){
  pwm_count <- mer_count(cnt_ob = pwmlist_leaf[[x]]@profileMatrix,dhs_gr = dhs_gr,by="pwmlist")
  readr::write_rds(pwm_count,paste0(stringr::str_replace_all(pwmlist_leaf[[x]]@ID," ",""),".",x,".Rds"))
  pwm_df <- data.frame(num=1:200,count=colSums(pwm_count))
  readr::write_rds(pwm_df,paste0(stringr::str_replace_all(pwmlist_leaf[[x]]@ID," ",""),"_df",x,".Rds"))
  png(paste0(stringr::str_replace_all(pwmlist_leaf[[x]]@ID," ",""),".png"))
  
  motif_p <- universalmotif::view_motifs(pwmlist_leaf[[1]],show.positions = F)+ylab("")+theme(axis.ticks.y =element_blank(),axis.line.y = element_blank(),axis.text.y = element_blank())
  
  h_line = ComplexHeatmap::HeatmapAnnotation(count=ComplexHeatmap::anno_lines(pwm_df$count,height = unit(2, "cm"),ylim = c(0,30)))
  plotht <- ComplexHeatmap::Heatmap(pwm_count,name="bykmer",
                                    column_title = "Genomic postion(-5kb,5kb)",
                                    column_title_side = "bottom",
                                    row_title = "DHSs in arabidophsis leaf cells",
                                    row_title_gp = gpar(fontsize = 14), #fontface = "bold"),
                                    clustering_distance_rows = function(m) dist(m,method = "maximum"),
                                    use_raster = T,cluster_columns = F,cluster_rows = F,
                                    raster_by_magick = TRUE, raster_resize_mat = max,
                                    col = circlize::colorRamp2(c(-1, 0, 1), c("red", "white", "blue")),
                                    width = unit(6, "cm"), height = unit(9, "cm"),
                                    top_annotation = h_line)
  
  ComplexHeatmap::draw(plotht)
  dev.off()
})




gg <- ggdraw()+draw_plot(y,0,0,1,1)+draw_plot(motif_p,0.1,0.1,0.5,0.5,1)
print(gg)


