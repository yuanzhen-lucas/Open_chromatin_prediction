bam="/wrk/luolin/work/Luolin_lv__h3k27me3/pre_processed/bamfiles/A1_1_zhuwenjun_p6_millipore_h3k27me3.sorted_rmdup.bam"
cov=GenomicRanges::coverage(bam)
txdb <- GenomicFeatures::makeTxDbFromGFF("/wrk/data/genome/yz_genome_data/aragenome/Athaliana.gff3")
trans <- GenomicFeatures::genes(txdb)
seqlevels(trans) %<>% str_replace_all("Chr","chr")
seqlevels(trans,pruning.mode="coarse")=seqlevels(trans)[seqlevels(trans) %>% str_detect("[0-9]$")]
strat <-  promoters(trans,upstream = 0,downstream = 1)
end <- strat
end@ranges@start <- trans@ranges@start + trans@ranges@width -1  %>% as.integer()
#strat@ranges@start <- strat@ranges@start - 1 %>% as.integer()
trim_bin <- function(tss,bin,by="end"){
  tss=resize(tss,bin,fix = by) %>% trim()
  tss=tss[width(tss)==bin]
  rngs=ranges(tss)
  start_=IRanges::start(rngs)
  filter_=start_>0
  tss=tss[filter_]
  return(tss)
}

covmatrix <- function(trans,cov){
  names(cov) %<>% str_replace("Chr","chr")
  cvg=cov[trans]
  cvgs=map(cvg %>% as.list,function(x){approx(x,n=1000)$y})
  cvgs %<>% unlist(use.names = F) %>% as.integer()
  cvg_mat=matrix(cvgs,ncol=1000,byrow = T)
  # flip genes on the reverse strand
  tss_minus=(strand(trans)=="-") %>% as.logical()
  cvg_mat[tss_minus,]= cvg_mat[tss_minus,ncol(cvg_mat):1]
  cvg_acc=colSums(cvg_mat)
  nor_ <- 501:1500
  matcov <- cbind(nor_ %>% as.numeric(),cvg_acc %>% as.numeric())
  return(matcov)
}

tssmat <- function(tss,cov,by="tss"){
  for (j in 1:length(cov)){
    cov[[j]]=c(cov[[j]],Rle(0L,500*2))
  }
  # ajust names according to names in tss
  names(cov) %<>% str_replace("Chr","chr")
  cov=cov[tss]
  cov %<>% unlist(use.names = F) %>% as.integer()
  cov=matrix(cov,ncol=500,byrow = T)
  # flip genes on the reverse strand
  tss_minus=(strand(tss)=="-") %>% as.logical()
  cov[tss_minus,]= cov[tss_minus,ncol(cov):1]
  cvg_acc=colSums(cov)
  if(by=="tss"){
    nor_ <- 1:500
  }else{
    nor_ <- 1501:2000
  }
  
  matcov <- cbind(nor_ %>% as.numeric(),cvg_acc %>% as.numeric())
  return(matcov)
}

start_nor <- trim_bin(strat,500)
end_nor <- trim_bin(end,500,by="start")
strat_ma <- tssmat(start_nor,cov,by="tss")
end_ma <- tssmat(end_nor,cov,by="tts")
gen_ma <- covmatrix(trans,cov)
test <- rbind(strat_ma,gen_ma,end_ma) %>% as.data.frame()
test <- test[order(test$V1),]
ggplot()+geom_line(test,mapping = aes(V1,V2))+scale_x_continuous(c(0,2000))+xlab("")+ylab("coverage")+geom_vline(xintercept=c(500,1500), linetype="dotted")
