print("============================================================")
print("Description:")
print("CROP-seq analysis: Generate Gene Expression Signature for gRNA-aggregated pseudo-bulk samples")
print("Compulsory Arguments:")
print("--workDir=/full/working/directory/where/results/will/be/generated")
print("--ExpFile=/full/path/to/pre-processed/expression/matrix")
print("--gRNAFile=/full/path/to/gRNA/assignment/file/generated/by/gRNA/assigning/script")
print("--jobName=/name/of/the/CROP-seq/run")
print("Optional Arguments:")
print("--Normalization=/whether/to/use/scran/or/tpm/to/normalize/expression/n/scran/as/default/")
print("--ges_method=/whether/to/use/Perturb-seq/or/Even_pool/method/to/generate/GES/n/Perturb-seq/as/default")
print("--ges_cal_method=/whether/to/use/MAD/or/z(Zscore)/for/Perturb-seq/or/zscore/ttest/or/mean/for/Even_pool/method/to/calculate/GES/n/Zscore/as/default")
print("--bootstrap=/number/of/bootstraps/n/default/is/1000")
print("--pool_size=/number/of/cells/in/a/pool/in/Even_pool/default/is/10")

print("=============================================================")
print("\n")
#manual

args <- commandArgs(trailingOnly = T)
###########################################################################
###/full/working/directory/where/results/will/be/generated
workDir <- args[grep("--workDir=", args)] 
workDir <- substr(workDir, 11, nchar(workDir))
print("working directory:")
print(workDir)
###########################################################################
###/full/path/to/pre-processed/expression/matrix
ExpFile <- args[grep("--ExpFile=", args)] 
ExpFile <- substr(ExpFile, 11, nchar(ExpFile))
print("ExpFile:")
print(ExpFile)
###########################################################################
###/full/path/to/gRNA/assignment/file/generated/by/gRNA/assigning/script
gRNAFile <- args[grep("--gRNAFile=", args)] 
gRNAFile <- substr(gRNAFile, 12, nchar(gRNAFile))
print("gRNAFile:")
print(gRNAFile)
###########################################################################
###/name/of/the/CROP-seq/run
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))
print("job name:")
print(jobName)
##########################################################################
###/whether/to/use/scran/or/tpm/to/normalize/expression
Normalization <- args[grep("--Normalization=", args)] 
Normalization <- substr(Normalization, 17, nchar(Normalization))
if (length(Normalization)==0) Normalization <- "scran"
print("Normalization method:")
print(Normalization)
##########################################################################
###/whether/to/use/Perturb-seq/or/Even_pool/method/to/generate/GES
ges_method <- args[grep("--ges_method=", args)] 
ges_method <- substr(ges_method, 14, nchar(ges_method))
if (length(ges_method)==0) ges_method <- "Perturb-seq"
print("ges_method:")
print(ges_method)
##########################################################################
###/whether/to/use/MAD/or/Zscore/method/to/calculate/GES/n/Zscore/as/default
ges_cal_method <- args[grep("--ges_cal_method=", args)] 
ges_cal_method <- substr(ges_cal_method, 18, nchar(ges_cal_method))
if (length(ges_method)==0) {
  if (ges_method=="Perturb-seq") {
    ges_cal_method <- "z"
  } else {
      ges_cal_method <- "zscore"
      }
}
print("ges_cal_method:")
print(ges_cal_method)
##########################################################################
###/number/of/bootstraps/n/default/is/1000
bootstrap <- args[grep("--bootstrap=", args)] 
bootstrap <- as.numeric(substr(bootstrap, 13, nchar(bootstrap)))
if (length(bootstrap)==0) bootstrap <- 1000
print("bootstrap:")
print(bootstrap)
##########################################################################
###pool_size=/number/of/cells/in/a/pool/in/Even_pool/default/is/10
pool_size <- args[grep("--pool_size=", args)] 
pool_size <- as.numeric(substr(pool_size, 13, nchar(pool_size)))
if (length(pool_size)==0) pool_size <- 10
print("pool_size:")
print(pool_size)
##########################################################################
CellNoThres <- args[grep("--CellNoThres=", args)] 
CellNoThres <- as.numeric(substr(CellNoThres, 15, nchar(CellNoThres)))
if (length(CellNoThres)==0) CellNoThres <- 1
print("CellNoThres:")
print(CellNoThres)
##########################################################################
Multi <- args[grep("--Multi=", args)] 
Multi <- as.logical(substr(Multi, 9, nchar(Multi)))
if (length(Multi)==0) Multi <- F
print("Multi:")
print(Multi)
#arguments


suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(scran))
set.seed(41)

##############################################################################################
### Function: MAD_sig
### Generate gene expression signature by MAD method.
### Return a vector representing gene expression signature.
### Args:
###       sample: a vector representing normalized gene expression data of the sample.
###       ref: a matrix representing normalized gene expression data of the reference set, 
###            rows representing genes and columns representing samples.
MAD_sig <- function(sample,ref) {
  return(sapply(rownames(ref), function(row){
    med <- median(ref[row,])
    mad <- mad(ref[row,])
    return((sample[row]-med)/mad)}))
}
##############################################################################################
### Function: z_sig
### Generate gene expression signature by z-score method.
### Return a vector representing gene expression signature.
### Args:
###       sample: a vector representing normalized gene expression data of the sample.
###       ref: a matrix representing normalized gene expression data of the reference set, 
###            rows representing genes and columns representing samples.
z_sig <- function(sample,ref) {
  return(sapply(rownames(ref), function(row){
    mean <- mean(ref[row,])
    sd <- sd(ref[row,])
    return((sample[row]-mean)/sd)}))
}
##############################################################################################
### Function: ges_perturb
### Generate gene expression signature matrix using Perturb-Seq (Weissman) method.
### Return a matrix representing gene expression signature.
### Args:
###       exprs: a matrix representing normalized gene expression data of the samples including negative controls, 
###              rows representing genes and columns representing samples.
###       gRNA_tags: a vector representing the gRNA tags of the samples in exprs
###       method: method for calculating GES, "MAD" or "z"; default = "MAD"
###       bootstrap: number of bootstraps, default = 1000
ges_perturb <- function(exprs,gRNA_tags,gRNAs,method=c("MAD","z"),bootstrap=1000) {
  negative_group <- grep("CONTROL",gRNA_tags,ignore.case = T)
  ges <- sapply(gRNAs, function(g){
    ### Get cell indexs. ###
    cells <- grep(paste0(g,""),gRNA_tags)
    ### Calculate ges. ####
    if (method=="MAD") {
      if (length(cells)==0) return(rep(0,times=nrow(exprs)))
      else if (length(cells)==1) {
        g_sum <- exprs[,cells]
        NC_pools <- sapply(1:bootstrap, function(i){
          return(exprs[,sample(negative_group,length(cells),replace=T)])
        })
        return(MAD_sig(g_sum,NC_pools))
      }
      else {
        g_sum <- rowSums(exprs[,cells])
        NC_pools <- sapply(1:bootstrap, function(i){
          return(rowSums(exprs[,sample(negative_group,length(cells),replace=T)]))})
        return(MAD_sig(g_sum,NC_pools))
      }
    }
    if (method=="z") {
      if (length(cells)==0) return(rep(0,times=nrow(exprs)))
      else if (length(cells)==1) {
        g_sum <- exprs[,cells]
        NC_pools <- sapply(1:bootstrap, function(i){
          return(exprs[,sample(negative_group,length(cells),replace=T)])
        })
        return(z_sig(g_sum,NC_pools))
      }
      else {
        g_sum <- rowSums(exprs[,cells])
        NC_pools <- sapply(1:bootstrap, function(i){
          return(rowSums(exprs[,sample(negative_group,length(cells),replace=T)]))})
        return(z_sig(g_sum,NC_pools))
      }
    }
  })
  return(ges)
}
##############################################################################################
### Function: ges_even_pool
### Generate gene expression signature matrix using even_pool method.
### Return a matrix representing gene expression signature.
### Args:
###       exprs: a matrix representing normalized gene expression data of the samples including negative controls, 
###              rows representing genes and columns representing samples.
###       gRNA_tags: a vector representing the gRNA tags of the samples in exprs
###       method: method for calculating GES, "zscore", "ttest" or "mean"; default = "zscore"
###       bootstrap: number of bootstraps, default = 1000
###       pool_size: the number of cells in each even pool, default = 10
ges_even_pool <- function(exprs,gRNA_tags,gRNAs,method=c("zscore","ttest","mean"),bootstrap=1000,pool_size=10) {
  negative_group <- grep("CONTROL",gRNA_tags,ignore.case = T)
  ges <- sapply(gRNAs, function(g){
    ### Get cell indexs. ###
    cells <- grep(paste0(g,""),gRNA_tags)
    ### Calculate ges. ####
    if (length(cells)==0) return(rep(0,times=nrow(exprs)))
    else {
      g_sum <- sapply(1:bootstrap, function(i){
        return(rowSums(exprs[,sample(cells,pool_size,replace=T)]))
      })
      NC_pools <- sapply(1:bootstrap, function(i){
        return(rowSums(exprs[,sample(negative_group,pool_size,replace=T)]))
      })
      sig <- viper::viperSignature(g_sum,NC_pools,method,verbose=T,seed=41)$signature
      return(apply(sig, 1, median,na.rm=T))
    }
  })
  return(ges)
}
###############################################################################################
### Function: ges_MASTDetRate
### Generate gene expression signature matrix using MAST with Detection Rate method.
### Return a matrix representing gene expression signature.
### Args:
###       exprs: a matrix representing normalized gene expression data of the samples including negative controls, 
###              rows representing genes and columns representing samples.
###       gRNA_tags: a vector representing the gRNA tags of the samples in exprs
ges_MASTDetRate <- function(exprs,gRNA_tags) {
  suppressPackageStartupMessages(library(MAST))
  suppressPackageStartupMessages(library(edgeR))
  group <- as.character(gRNA_tags)
  group[grep("CONTROL",group,ignore.case = T)] <- "CONTROL"
  group <- gsub("-","_",group)
  dge <- DGEList(exprs, group = group)
  dge <- calcNormFactors(dge)
  cdr <- as.numeric(scale(colMeans(exprs > 0)))
  cpms <- cpm(dge)
  sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                    cData = data.frame(wellKey = colnames(cpms), 
                                       group = group, cdr = cdr))
  zlmdata <- zlm(~0 + cdr + group, sca)
  ges <- sapply(sort(unique(grep("CONTROL",gRNA_tags,invert=T,value=T,ignore.case = T))), function(g){
    mast <- lrTest(zlmdata, Hypothesis(paste0("group",gsub("-","_",g),"-groupCONTROL")))
    mast[,"hurdle",c("lambda","Pr(>Chisq)")]
    contrast0 <- rep(0,times=ncol(coef(zlmdata,"C")))
    names(contrast0) <- colnames(coef(zlmdata,"C"))
    contrast0[grep("CONTROL",names(contrast0),ignore.case = T)] <- 1
    contrast1 <- diag(ncol(coef(zlmdata,"C")))
    rownames(contrast1)<-colnames(contrast1)<-colnames(coef(zlmdata,"C"))
    FC <- getLogFC(zlmdata,contrast0,contrast1)
    FC <- FC[which(FC[,2]==paste0("group",gsub("-","_",g))),"logFC"]
    mast[,"hurdle","lambda"] <- qchisq(mast[,"hurdle","Pr(>Chisq)"],mast[,"hurdle","df"],lower.tail = F)*sign(FC)
  })
  lambda <- ges[1:nrow(exprs),]
  rownames(lambda) <- rownames(exprs)
  PValue <- ges[(nrow(exprs)+1):(2*nrow(exprs)),]
  rownames(PValue) <- rownames(exprs)
  ges <- list(lambda,PValue)
  names(ges) <- c("Signature_Chi","PValue")
  
  return(ges)
}
###############################################################################################
### Function: ges_edgeRQLFDetRate
### Generate gene expression signature matrix using edgeRQLF with Detection Rate method.
### Return a matrix representing gene expression signature.
### Args:
###       exprs: a matrix representing normalized gene expression data of the samples including negative controls, 
###              rows representing genes and columns representing samples.
###       gRNA_tags: a vector representing the gRNA tags of the samples in exprs
ges_edgeRQLFDetRate <- function(exprs,gRNA_tags) {
  suppressPackageStartupMessages(library(edgeR))
  group <- as.character(gRNA_tags)
  group[grep("CONTROL",group,ignore.case = T)] <- "CONTROL"
  dge <- DGEList(exprs, group = group)
  dge <- calcNormFactors(dge)
  cdr <- scale(colMeans(exprs > 0))
  design <- model.matrix(~ 0 + cdr + group)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design=design)
  ges <- sapply(sort(unique(grep("CONTROL",gRNA_tags,invert=T,value=T,ignore.case = T))), function(g){
    contrast <- rep(0,ncol(design))
    contrast[grep("CONTROL",colnames(design),ignore.case = T)] <- -1
    contrast[grep(paste0(g,""),colnames(design),ignore.case = T)] <- 1
    qlf <- glmQLFTest(fit,contrast = contrast)
    tt <- topTags(qlf, n = Inf,sort.by="none")$table[,c("F","logFC","PValue")]
  })
  F <- sapply(ges[1,],function(v){v})
  rownames(F) <- rownames(exprs)
  Sign <- sapply(ges[2,],function(v){v})
  F <- F * sign(Sign)
  PValue <- sapply(ges[3,],function(v){v})
  rownames(PValue) <- rownames(exprs)
  ges <- list(F,PValue)
  names(ges) <- c("Signature_F","PValue")
  return(ges)
}
###############################################################################################
heatmap_all <- function(ges,rownames,row_symbols,ges_method) {
  ges_values <- abs(ges[rownames,])
  limit <- quantile(ges_values[is.finite(ges_values)],probs = 1,na.rm = T)
  pdf(paste0(jobName,ges_cal_method,"_all.pdf"),width = 7,height = 7)
  pheatmap(ges[rownames,],
           color = colorRampPalette(c("blue","white","red"))(100),
           cluster_rows=F,cluster_cols=F,
           labels_row = row_symbols, main=paste0(jobName,"_",ges_method),
           breaks = c(seq(-limit,0,length.out = 51),
                      seq(limit/100,
                          limit,length.out=50))
  )
  dev.off()
}

heatmap_self <- function(ges,rownames,row_symbols,ges_method) {
  pdf(paste0(jobName,ges_cal_method,"_self.pdf"),width = 7,height = 7)
  ges_self <- ges[rownames,colnames(ges)[-grep("Control",colnames(ges),ignore.case = T)]]
  for (i in 1:nrow(ges_self)) {
    for (j in 1:ncol(ges_self)) {
      if (strsplit(colnames(ges_self)[j],split = '-')[[1]][1]==row_symbols[i]) print(i)
      else ges_self[i,j] <- 0
    }
  }
  ges_self_values <- abs(ges_self)
  limit_self <- quantile(ges_self_values[is.finite(ges_self_values)],probs = 1,na.rm = T)
  pheatmap(ges_self, 
           color = colorRampPalette(c("blue","white","red"))(100),
           cluster_rows=F,cluster_cols=F,
           labels_row = row_symbols,main=paste0(jobName,"_",ges_method),
           breaks = c(seq(-limit_self,0,length.out = 51),
                      seq(limit_self/100,
                          limit_self,length.out=50))
  )
  dev.off()
}




setwd(workDir)
### Load exprs from preliminary analysis.
load(ExpFile)
load(gRNAFile)

### Take cells with only 1 gRNA assigned. #######################################################
if (Multi) {
  exprs <- exprs[,sapply(gRNA_assigned,length)>0]
  gRNA_assigned <- gRNA_assigned[sapply(gRNA_assigned,length)>0]
} else {
  exprs <- exprs[,sapply(gRNA_assigned,length)==1]
  gRNA_assigned <- gRNA_assigned[sapply(gRNA_assigned,length)==1]
}





if (Normalization=="tpm") {
  exprs_norm <- apply(exprs,2,function(v){1e6*v/sum(v)})
} else {
  s <- computeSumFactors(exprs,c(20,40,60,80,100))
  exprs_norm <- t(t(exprs)/s*mean(s))
}


gRNAs <- rownames(exprs_gRNA)
gRNAs <- gRNAs[sapply(gRNAs, function(g){length(grep(paste0(g,""),gRNA_assigned))})>=CellNoThres]
# Get cell ids of each MR group.
#groups <- sapply(gRNAs, function(s){
#  grep(s,gRNA_assigned)
#})
#negative_group <- grep("CONTROL",gRNA_assigned,ignore.case = T)
#if(length(negative_group)<100) negative_group <- grep("",gRNA_assigned)


MRs <- unique(sapply(gRNAs[-grep("Control",gRNAs,ignore.case = T)],function(g){strsplit(g,split = '-')[[1]][1]}))
MR_tags_ensg <- mapIds(org.Hs.eg.db,as.character(MRs),column = "ENSEMBL",keytype = "SYMBOL",
                       multiVals=function(n){intersect(n,rownames(exprs))[1]})
inter_rows <- MR_tags_ensg[MR_tags_ensg%in%rownames(exprs)]


if (ges_method=="Perturb-seq") {
  ges_perturb <- ges_perturb(exprs_norm,gRNA_assigned,gRNAs,method=ges_cal_method,bootstrap)
  rownames(ges_perturb) <- rownames(exprs_norm)
  save(ges_perturb,file=paste0(jobName,"_",ges_cal_method,".rda"))
  ### GES_perturb ###
  heatmap_all(ges_perturb,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  heatmap_self(ges_perturb,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  
  
} else if (ges_method=="Even_pool") {
  ges_even_pool <- ges_even_pool(exprs_norm,gRNA_assigned,gRNAs,method=ges_cal_method,bootstrap,pool_size)
  rownames(ges_even_pool) <- rownames(exprs_norm)
  save(ges_even_pool,file=paste0(jobName,"_",ges_cal_method,".rda"))
  ### GES_even_pool ###
  heatmap_all(ges_even_pool,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  heatmap_self(ges_even_pool,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  
  
} else if (ges_method=="MASTDetRate") {
  ges_MASTDetRate <- ges_MASTDetRate(exprs,gRNA_assigned)
  save(ges_MASTDetRate,file=paste0(jobName,"_",ges_cal_method,".rda"))
  ges_MASTDetRate <- ges_MASTDetRate$Signature_lambda
  ### ges_MASTDetRate ###
  heatmap_all(ges_MASTDetRate,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  heatmap_self(ges_MASTDetRate,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  
  
} else if (ges_method=="edgeRQLFDetRate") {
  ges_edgeRQLFDetRate <- ges_edgeRQLFDetRate(exprs,gRNA_assigned)
  save(ges_edgeRQLFDetRate,file=paste0(jobName,"_",ges_cal_method,".rda"))
  ges_edgeRQLFDetRate <- ges_edgeRQLFDetRate$Signature_F
  ### ges_edgeRQLFDetRate ###
  heatmap_all(ges_edgeRQLFDetRate,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
  heatmap_self(ges_edgeRQLFDetRate,inter_rows,MRs[MR_tags_ensg%in%rownames(ges_perturb)],ges_method = ges_method)
}
