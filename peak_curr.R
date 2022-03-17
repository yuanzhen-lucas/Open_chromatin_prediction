library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- GenomicFeatures::makeTxDbFromGFF("/wrk/data/genome/yz_genome_data/aragenome/Athaliana.gff3")
genes <- GenomicFeatures::genes(txdb)
peaks <- ChIPseeker::readPeakFile("/wrk/yuanzhen/dnase_motif/peak_curra/bamfile/root.pos_peaks.narrowPeak")
#peakAnno <- ChIPseeker::annotatePeak(peaks,TxDb=txdb)
peaks_pos <- ChIPseeker::readPeakFile("/wrk/data/ref_data_yz/human/erythroid/bamfile/fseq2test/fseq2_result_peaks.narrowPeak")
ChIPseeker::plotPeakProf2(peak = peaks_pos, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 100,
             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, weightCol = "V5",ignore_strand = F)


genpk <- data.frame(gene_id=peakAnno@anno@elementMetadata$geneId)
genpk$peak_start <- peakAnno@anno@ranges@start
genpk$peak_end <- peakAnno@anno@ranges@start + peakAnno@anno@ranges@width - 1
genpk$peak_len <- peakAnno@anno@ranges@width
genpk$gene_start <- peakAnno@anno@elementMetadata$geneStart
genpk$gene_end <- peakAnno@anno@elementMetadata$geneEnd
genpk$gene_len <- peakAnno@anno@elementMetadata$geneLength
genpk <- dplyr::filter(genpk,peak_len>50)

x_tmp <- function(peak_start,peak_end,gene_start,gene_end,gene_len){
  if(peak_start < gene_start && peak_end < gene_end && peak_end > gene_start){
    if(peak_start <= gene_start-500){
      peak_start = 
    }else{
      peak_start = 500 - (gene_start - peak_start)
    }
    peak_end <- 500 + (gene_len * (peak_end-gene_start) / 1000)
  }else if(peak_start > gene_start && peak_start < gene_end && peak_end > gene_end){
    if(peak_end >= gene_end + 500){
      peak_end = 2000
    }else{
      peak_end = peak_end - gene_end + 1500
    }
    peak_start = 1500 - (gene_len * (gene_end -peak_start) / 1000)
  }else if(peak_start > gene_start && peak_end < gene_end){
    peak_start = 500 + (1000 / gene_len * (peak_start-gene_start))
    peak_end <- 500 + (1000 / gene_len * (peak_end-gene_start))
  }else if(peak_start < gene_start && peak_end < gene_start){
    if(peak_start > gene_start -500 && peak_end > gene_start -500){
      peak_start = 500 - (gene_start - peak_start)
      peak_end = 500 - (gene_start - peak_end)
    }else if(peak_start < gene_start -500 && peak_end > gene_start -500){
      peak_start = 0
      peak_end = 500 - (gene_start - peak_end)
    }else{
      peak_start = 0
      peak_end = 0
    }
    
  }else{
    if(peak_start < gene_end + 500 && peak_end < gene_end + 500){
    peak_start = peak_start - gene_end + 1500
    peak_end = peak_end - gene_end + 1500
  }else if(peak_start < gene_end + 500 && peak_end > gene_end + 500){
    peak_start = peak_start - gene_end + 1500
    peak_end = 2000
  }else{
    peak_start = 2000
    peak_end = 2000
  }
  }

  c(peak_start,peak_end)
  
}
dt_p = t(mapply(x_tmp, genpk$peak_start,genpk$peak_end,genpk$gene_start,genpk$gene_end,genpk$gene_len)) %>% as.data.frame()
dt_p_t = filter(dt_p,V1>0 & V2>0)
dt_p_t = filter(dt_p_t,V1<2000 & V2<2000)

for(m in 1:2000){
  x=0
  for(n in 1:nrow(dt_p_t)){
    if(m>dt_p_t[n,1]&&m<dt_p_t[n,2]){
      x=x+1
    }
    
    if(n==nrow(dt_p_t)){
      cat(m,x,"\n",file = "test.txt",append = T)
    }
  }
  
}

test <- read.table("~/dnase_motif/peak_curra/test.txt", quote="\"", comment.char="")
plot(x=test$V1,y=test$V2)
