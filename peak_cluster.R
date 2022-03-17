library(tidyverse)
library(magrittr)

covmatrix <- function(trans,cov){
  names(cov) %<>% stringr::str_replace("Chr","chr")
  cvg=cov[trans]
  cvg %<>% unlist(use.names = F) %>% as.integer()
  cvg_mat=matrix(cvg,ncol=length(trans),byrow = T)
  # flip genes on the reverse strand
  #tss_minus=(strand(trans)=="-") %>% as.logical()
  #cvg_mat[tss_minus,]= cvg_mat[tss_minus,ncol(cvg_mat):1]
  cvg_acc=colSums(cvg_mat)
  return(cvg_acc)
}

filelist <- Sys.glob("/wrk/luolin/work/Luolin1_4__araDNase/pre_processed/bamfiles/*.bam")
covlist <- BiocParallel::bplapply(filelist,function(x){
  genome_ara <- Biostrings::readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  chrname <- names(genome_ara) %>% str_replace_all("Chr","chr")
  chrlength <- Biostrings::width(genome_ara)
  seq <- lapply(1:5,function(x){
    gourpl <- floor(chrlength[[x]] / 1000)
    test <- data.frame(matrix(NA,gourpl+1,2))
    for(a in 1:gourpl){
      test[a,1] = chrname[[x]]
      test[a,2] = 1000*(a-1)+1
      test[a,3] = 1000 *a
    }
    test[gourpl,1]=chrname[[x]]
    test[gourpl,2]=gourpl*1000
    test[gourpl,3]= chrlength[[x]]
    test
  })
  sqinfo <- do.call("rbind",seq) %>% na.omit()
  gr1 <- GenomicRanges::GRanges(seqnames=sqinfo$X1,ranges=IRanges::IRanges(start=sqinfo$X2, end = sqinfo$V3))
  library(magrittr)
  y=GenomicRanges::coverage(x)
  covmatrix(gr1,y)
  },BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

sapldna <- do.call(cbind,covlist) %>% as.matrix()
colnames(sapldna) <-  paste0(str_sub(basename(filelist),1,str_length(basename(filelist))-17),"_","DNase")


filemn <- Sys.glob("/wrk/wenchenjin/work/Chenjin1_2__MN_ara_b2/pre_processed/bamfiles/*.bam")
covlistmn <- BiocParallel::bplapply(filemn,function(x){
  genome_ara <- Biostrings::readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  chrname <- names(genome_ara) %>% str_replace_all("Chr","chr")
  chrlength <- Biostrings::width(genome_ara)
  seq <- lapply(1:5,function(x){
    gourpl <- floor(chrlength[[x]] / 1000)
    test <- data.frame(matrix(NA,gourpl+1,2))
    for(a in 1:gourpl){
      test[a,1] = chrname[[x]]
      test[a,2] = 1000*(a-1)+1
      test[a,3] = 1000 *a
    }
    test[gourpl,1]=chrname[[x]]
    test[gourpl,2]=gourpl*1000
    test[gourpl,3]= chrlength[[x]]
    test
  })
  sqinfo <- do.call("rbind",seq) %>% na.omit()
  gr1 <- GenomicRanges::GRanges(seqnames=sqinfo$X1,ranges=IRanges::IRanges(start=sqinfo$X2, end = sqinfo$V3))
  library(magrittr)
  y=GenomicRanges::coverage(x)
  covmatrix(gr1,y)
},BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

saplmn <- do.call(cbind,covlistmn) %>% as.matrix()
colnames(saplmn) <- paste0(str_sub(basename(filemn),1,str_length(basename(filemn))-17),"_","MNase")


filepos <- Sys.glob("/wrk/yuanzhen/dnase_motif/peak_curra/bamfile/*.bam")
covlistpos <- BiocParallel::bplapply(filepos,function(x){
  genome_ara <- Biostrings::readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  chrname <- names(genome_ara) %>% str_replace_all("Chr","chr")
  chrlength <- Biostrings::width(genome_ara)
  seq <- lapply(1:5,function(x){
    gourpl <- floor(chrlength[[x]] / 1000)
    test <- data.frame(matrix(NA,gourpl+1,2))
    for(a in 1:gourpl){
      test[a,1] = chrname[[x]]
      test[a,2] = 1000*(a-1)+1
      test[a,3] = 1000 *a
    }
    test[gourpl,1]=chrname[[x]]
    test[gourpl,2]=gourpl*1000
    test[gourpl,3]= chrlength[[x]]
    test
  })
  sqinfo <- do.call("rbind",seq) %>% na.omit()
  gr1 <- GenomicRanges::GRanges(seqnames=sqinfo$X1,ranges=IRanges::IRanges(start=sqinfo$X2, end = sqinfo$V3))
  library(magrittr)
  y=GenomicRanges::coverage(x)
  covmatrix(gr1,y)
},BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

saplpos <- do.call(cbind,covlistpos) %>% as.matrix()
colnames(saplpos) <- paste0(str_sub(basename(filepos),1,str_length(basename(filepos))-17),"_","Positive")

fileatac <- Sys.glob("/wrk/yuanzhen/dnase_motif/peak_curra/raw_data/atacbam/*.bam")
covlistatac <- BiocParallel::bplapply(fileatac,function(x){
  genome_ara <- Biostrings::readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  chrname <- names(genome_ara) %>% str_replace_all("Chr","chr")
  chrlength <- Biostrings::width(genome_ara)
  seq <- lapply(1:5,function(x){
    gourpl <- floor(chrlength[[x]] / 1000)
    test <- data.frame(matrix(NA,gourpl+1,2))
    for(a in 1:gourpl){
      test[a,1] = chrname[[x]]
      test[a,2] = 1000*(a-1)+1
      test[a,3] = 1000 *a
    }
    test[gourpl,1]=chrname[[x]]
    test[gourpl,2]=gourpl*1000
    test[gourpl,3]= chrlength[[x]]
    test
  })
  sqinfo <- do.call("rbind",seq) %>% na.omit()
  gr1 <- GenomicRanges::GRanges(seqnames=sqinfo$X1,ranges=IRanges::IRanges(start=sqinfo$X2, end = sqinfo$V3))
  library(magrittr)
  y=GenomicRanges::coverage(x)
  covmatrix(gr1,y)
},BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

saplatac <- do.call(cbind,covlistatac) %>% as.matrix()
colnames(saplatac) <- paste0(str_sub(basename(fileatac),1,str_length(basename(fileatac))-4),"_","Atac")

filejimpos <- Sys.glob("/wrk/yuanzhen/dnase_motif/peak_curra/raw_data/dhs_jm/bamfile/*.bam")
covlistjimpos <- BiocParallel::bplapply(filejimpos,function(x){
  genome_ara <- Biostrings::readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  chrname <- names(genome_ara) %>% str_replace_all("Chr","chr")
  chrlength <- Biostrings::width(genome_ara)
  seq <- lapply(1:5,function(x){
    gourpl <- floor(chrlength[[x]] / 1000)
    test <- data.frame(matrix(NA,gourpl+1,2))
    for(a in 1:gourpl){
      test[a,1] = chrname[[x]]
      test[a,2] = 1000*(a-1)+1
      test[a,3] = 1000 *a
    }
    test[gourpl,1]=chrname[[x]]
    test[gourpl,2]=gourpl*1000
    test[gourpl,3]= chrlength[[x]]
    test
  })
  sqinfo <- do.call("rbind",seq) %>% na.omit()
  gr1 <- GenomicRanges::GRanges(seqnames=sqinfo$X1,ranges=IRanges::IRanges(start=sqinfo$X2, end = sqinfo$V3))
  library(magrittr)
  y=GenomicRanges::coverage(x)
  covmatrix(gr1,y)
},BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

sapljimpos <- do.call(cbind,covlistjimpos) %>% as.matrix()
colnames(sapljimpos) <- paste0(str_sub(basename(filejimpos),1,str_length(basename(filejimpos))-4),"_","Positive")






saplmn_n <- saplmn %>% sweep(2,colSums(.),FUN = "/")
sapldna_n <- sapldna %>% sweep(2,colSums(.),FUN = "/")
sapldna_l <- -log2(sapldna_n)
saplmn_l <- -log2(saplmn_n)
saplpos_n <- saplpos %>% sweep(2,colSums(.),FUN = "/")
saplpos_l <- -log2(saplpos_n)
saplatac_n <- saplatac %>% sweep(2,colSums(.),FUN = "/")
saplatac_l <- -log2(saplatac_n)
sapljimpos_n <- sapljimpos %>% sweep(2,colSums(.),FUN = "/")
sapljimpos_l <- -log2(sapljimpos_n)

total_data=cbind(sapldna_l,saplmn_l,saplpos_l,saplatac_l,sapljimpos_l)
total_data[is.infinite(total_data)] <- 1
ggcor::quickcor(sapldna_l,saplmn_l,cor.test = TRUE)+geom_square()
ggcor::quickcor(total_data,cor.test = TRUE,cluster=T)+geom_square()+scale_color_gradientn(colours = c("blue","white","red"))

pca_test <- t(total_data) %>% as.data.frame()
pca_test$grp <- c(rep("DNase",8),rep("MNase",17),rep("Positive",2),rep("ATAC",2),rep("Positive",5))
pca_data <- stats::prcomp(pca_test[,-119145])
factoextra::fviz_screeplot(pca_data, addlabels = TRUE)
factoextra::fviz_pca_ind(pca_data, col.ind="cos2", 
                         geom = c("point","text"), # show points only
                         gradient.cols = c("white", "#2E9FDF", "#FC4E07" ))
fviz_pca_ind(pca_data, label="none", habillage=pca_test$grp,
             addEllipses=TRUE, ellipse.level=0.95,
             palette = "Dark2")

fviz_pca_biplot(pca_data, label = "var", habillage=pca_test$grp,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_bw()) +
  theme(panel.grid = element_blank()) +
  scale_color_brewer(palette = "Set1")
