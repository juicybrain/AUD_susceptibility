
library("tidyverse")
library("ggbiplot")
library(DESeq2)
library("VennDiagram")

setwd("C:/Users/Yx Li/Documents/02RNAseqHisat2_FC/rawdata")

files <- list.files()
 dat <- read.table(list.files()[1],sep='\t',header = TRUE,skip = 1)
for (i in 2:length(files))
     {
       da <- read.csv(files[i],sep='\t',header = TRUE,skip = 1)
       dat <- merge(dat,da, by="Geneid")
}

#colnames(dat)[2:8]=str_split_fixed(colnames(dat)[2:8],"\\.",6)[,1]
#samples <- cbind(colnames(dat)[2:8],str_split_fixed(colnames(dat)[2:8],"_",3)[,2])

# DEseq2 allsamples


spl_group <- c("ct","ko","ct","ko","ct","ko","ct","ko","ct","ko","ct","ko","ct","ko","ko","ct")
cond=factor(spl_group)

datf_matrix <- as.matrix(dat[,-1])
rownames(datf_matrix) <- dat$Geneid

dds <- DESeqDataSetFromMatrix(datf_matrix,DataFrame(cond),~cond)
p1 <- plotPCA(rlog(dds),ntop=500, intgroup="cond")+theme_bw()
p1
dDEG <- DESeq(dds)

res <- results(dDEG)
res.ord <- res[order(res$padj),]
res.ord.sgnft <- subset(res.ord,padj <=1)


gEx <- as.data.frame(res.ord.sgnft)

gEx$threshold <- ifelse(gEx$padj<0.05 & abs(gEx$log2FoldChange)>=1,ifelse(gEx$log2FoldChange>1,"up","down"),"NOT")

dds <- estimateSizeFactors(dds)
normed_counts <- as.data.frame(counts(dds, normalized=TRUE))


p2 <- plotCounts(dDEG, gene=which.min(res$padj),intgroup = "cond",pch=2, col=dDEG$cond,cex=3,transform = T)

p2

library(stringr)
# PCA
#plotPCA(rlog(dds), intgroup="condition")+theme_bw()
alltable <- normed_counts
 inputmatrix <- filter(alltable,rowMeans(alltable)>1)
inputmatrix_1 <- alltable[names(which(apply(alltable,1,mean)>1)),]


allpca <- prcomp(t(inputmatrix),scale=TRUE)
#allpca_1 <- prcomp(t(inputmatrix_1),scale=T)

g1 <- ggbiplot(allpca, obs.scale = 2, var.scale = 1,var.axes = F,
    groups = spl_group, ellipse = TRUE, circle = TRUE) +
    scale_color_discrete(name = '') +geom_text(aes(label=c(seq(1:length(spl_group)))),hjust=-0.5,vjust=-0.5,size=3)+
    theme(legend.direction = 'horizontal', legend.position = 'top')
g1

#ggbiplot(allpca_1, obs.scale = 2, var.scale = 1,var.axes = F,groups = c("-","-","+","-","+","+","+"), ellipse = TRUE, circle = TRUE) +scale_color_discrete(name = '') +theme(legend.direction = 'horizontal', legend.position = 'top')

probe_var <- apply(inputmatrix,1,var)
order_var <- order(probe_var,decreasing = TRUE)[1:200]
pca_var <- prcomp(t(inputmatrix[order_var,]),scale=TRUE)
g2 <- ggbiplot(pca_var, obs.scale = 2, var.scale = 1,var.axes = F, groups = spl_group, ellipse = TRUE, circle = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', legend.position = 'top') +geom_text(aes(label=c("A1","B1","A2","B2","A3","B3","A4","B4","A5","B5","A6","B6","A7","B7","B8","A8")),hjust=-0.5,vjust=-0.5,size=3)+ggtitle("PCA of TET1 D1 male cell RNAseq Data")

g2


tet <- c('ENSMUSG00000047146.17','ENSMUSG00000040943.12','ENSMUSG00000034832.15')
#TET_Ex <- dat %>% filter(Geneid %in% tet )
TET_Ex <- alltable[tet,]

probesetvar = apply(TET_Ex[,],1,var)
ord = order(probesetvar,decreasing=TRUE) 

 
 pca_TET_ex <- prcomp(t(TET_Ex[ord,]),scale=TRUE)
 pca_TET_ex_dat <- as.data.frame(pca_TET_ex$x[,1:2])
 pca_TET_ex_dat$name <- row.names(pca_TET_ex_dat)
 pca_TET_ex_dat$group <- spl_group
 
 
 g3 <- ggplot(pca_TET_ex_dat,aes(x=PC1,y=PC2),label=pca_TET_ex_dat$name) + geom_point(aes(color=group)) + geom_text(aes(label=name),hjust=-0.5,vjust=-0.5,size=3)
 g3
g4 <- ggbiplot(pca_TET_ex, obs.scale = 2, var.scale = 1,var.axes = F,
    groups = spl_group, ellipse = TRUE, circle = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', legend.position = 'top') +geom_text(aes(label=c("A1","B1","A2","B2","A3","B3","A4","B4","A5","B5","A6","B6","A7","B7","B8","A8")),hjust=-0.5,vjust=-0.5,size=3)+ggtitle("PCA of TET1 D1 male cell RNAseq Data")
g4


## select ko (4,5,7,8) and all ct

dat_DE <- dat[,-c(1,3,5,7,13)]
row.names(dat_DE) <- dat$Geneid
dat_DE <- dat_DE[which(apply(dat_DE,1,mean)>1),]
#alltable[names(which(apply(alltable,1,mean)>1)),]



cond=factor(c("ct","ct","ct","ct","ko","ct","ko","ct","ct","ko","ko","ct"))

datf_matrix <- as.matrix(dat_DE)
rownames(datf_matrix) <- rownames(dat_DE)

dds <- DESeqDataSetFromMatrix(datf_matrix,DataFrame(cond),~cond)

#PCA
p5 <- plotPCA(rlog(dds),ntop=2000, intgroup="cond")+theme_bw()
p5

dds_pca <- estimateSizeFactors(dds)
normed_counts <- as.data.frame(counts(dds_pca, normalized=TRUE))

probesetvar = apply(normed_counts,1,var)
ord = order(probesetvar,decreasing=TRUE)[1:2000] 

 
 pca_TET_ex <- prcomp(t(datf_matrix[ord,]),scale=TRUE)
 pca_TET_ex_dat <- as.data.frame(pca_TET_ex$x[,1:2])
 pca_TET_ex_dat$name <- row.names(pca_TET_ex_dat)
 pca_TET_ex_dat$group <- c("ct","ct","ct","ct","ko","ct","ko","ct","ct","ko","ko","ct")

PCA_selected <- ggbiplot(pca_TET_ex, obs.scale = 1,groups= c("ct","ct","ct","ct","ko","ct","ko","ct","ct","ko","ko","ct"), var.scale = 1,var.axes = F, ellipse = TRUE, circle = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', legend.position = 'top') +geom_text(aes(label=c("A1","A2","A3","A4","B4","A5","B5","A6","A7","B7","B8","A8")),hjust=-0.5,vjust=-0.5,size=3)+ggtitle("PCA of TET3 KO  RNAseq Data")

PCA_selected




dDEG <- DESeq(dds)

res <- results(dDEG)
res.ord <- res[order(res$padj),]
res.ord.sgnft <- subset(res.ord,padj <=1)

#file1 <- "C:/Users/Yx Li/Documents/01HXuRNAseq/D1_12vs367.txt"
#write.table(res.ord.sgnft,file=file1,sep="\t",quote=F)

gEx_1 <- as.data.frame(res.ord.sgnft)

gEx_1$threshold <- ifelse(gEx_1$padj<0.05 & abs(gEx_1$log2FoldChange)>=0.2,ifelse(gEx_1$log2FoldChange>=0.2,"up","down"),"NOT")



###Volcano###


p6 <- ggplot(data = gEx_1, aes(x = log2FoldChange, y = -log10(padj), color=threshold)) +
  geom_point(alpha=0.8, size = 1) +
  xlim(c(-7, 7)) +
  ylim(c(0, 5)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(
    legend.position="right",
    panel.grid=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_color_manual(name = "", values = c("blue", "red", "grey"), limits = c("up", "down", "NOT"))+
  labs(x="log2 (Fold Change)",y="-log10 (P-adjusted value)",title="NAc Tet3 WT vs KO")


p6
# DESeq2 inneit function
lfcshrink_res_ddds <- lfcShrink(dDEG, coef=2, res=res)


## Sort the fold changes by adjusted p-value

lfcshrink_res_ddds_ordered=lfcshrink_res_ddds[order(lfcshrink_res_ddds$padj),]

## Now plot the fold changes. 
p7 <- plotMA(lfcshrink_res_ddds_ordered)
p7


res.ord.sgnft <- subset(lfcshrink_res_ddds_ordered,padj <=1)

#file1 <- "C:/Users/Yx Li/Documents/01HXuRNAseq/D1_12vs367.txt"
#write.table(res.ord.sgnft,file=file1,sep="\t",quote=F)

gEx <- as.data.frame(res.ord.sgnft)

gEx$threshold <- ifelse(gEx$padj<0.05 & abs(gEx$log2FoldChange)>=0.2,ifelse(gEx$log2FoldChange>0.2,"up","down"),"NOT")


p8 <- ggplot(data =gEx, aes(x = log2FoldChange, y = -log10(padj), color=threshold)) +
  geom_point(alpha=0.8, size = 1) +
  xlim(c(-1, 1)) +
  ylim(c(0, 6)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(
    legend.position="right",
    panel.grid=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_color_manual(name = "", values = c("blue", "red", "grey"), limits = c("up", "down", "NOT"))+
  labs(x="log2 (Fold Change)",y="-log10 (P-value)",title="Drd1-WT vs KO")


p8


#library("TxDb.Mmusculus.UCSC.mm10.knownGene")

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#gene_id <- read.table("C:/Users/Yx Li/Documents/mm10_genecode_id")

#colnames(gene_id) <- c("V2","V1","V3","V4","V5")
#gEx$V1 <- row.names(gEx)
 
#u_d_genes <- gEx %>% filter(gEx$threshold!="NOT")

#u_d_genes_ids <- merge(u_d_genes,gene_id,by="V1",all=FALSE)


#dat$V1 <- rownames(dat)

#u_d_genes_ids_ex <- merge(dat,u_d_genes_ids,by="V1",all=FALSE)

#file3 <- "C:/Users/Yx Li/Documents/01HXuRNAseq/D1_12vs367_passed_0.2_0.05_ids.txt"

#write.table(u_d_genes_ids_ex,file=file3,sep="\t",quote=F,row.names = FALSE)




# annotation

#library(annotables)
library(annotables)
library("ggrepel")

gEx_1$V1 <- row.names(gEx_1)
gEx_1$gene <- str_split_fixed(gEx_1$V1,"\\.",2)[,1]
# also can use: substr(gEx$V1,1,18)  
test <- gEx_1[which(gEx_1$threshold!="NOT"),] %>% 
  arrange(padj) %>% 
#  head(50) %>% 
  inner_join(grcm38, by=c("gene"="ensgene")) 
 # dplyr::select(gene,entrez, log2FoldChange, padj, symbol, description,threshold) 

upgenes <- gEx_1[which(gEx_1$threshold=="up"),] %>% 
  arrange(padj) %>% 
#  head(50) %>% 
  inner_join(grcm38, by=c("gene"="ensgene")) 


downgenes <- gEx_1[which(gEx_1$threshold=="down"),] %>% 
  arrange(padj) %>% 
#  head(50) %>% 
  inner_join(grcm38, by=c("gene"="ensgene")) 




gEx$V1 <- row.names(gEx)
gEx$gene <- str_split_fixed(gEx$V1,"\\.",2)[,1]
test2 <- gEx[which(gEx$threshold!="NOT"),] %>% 
  arrange(padj) %>% 
#  head(50) %>% 
  inner_join(grcm38, by=c("gene"="ensgene")) %>% 
  dplyr::select(gene,entrez, log2FoldChange, padj, symbol, description,threshold) 

p8_plus <- p8+ geom_text_repel(data=head(test2,5), aes(label=symbol))


library("clusterProfiler") 
library("org.Mm.eg.db") 

#convert genesymbol to Entrez geneid
#x <- org.Mm.egSYMBOL2EG
#mapped_genes <- mappedkeys(x)
#xx <- as.list(x[mapped_genes])
#genename <- u_d_genes_ids_ex_2$V3
#u_d_genes_ids_ex_2$geneid <- as.character(genename)%>% map(~xx[[.x]])
#gene_down <- u_d_genes_ids_ex_2%>% filter(u_d_genes_ids_ex_2$threshold=="up")


DMRgo_up <- enrichGO(gene = na.omit(test$entrez), OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

p9 <- dotplot(DMRgo_up,title="Tet3 ko upregulated genes GO terms") 
p9
compKEGG_up <- compareCluster(geneCluster = test, fun = "enrichKEGG", organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH") 
p10 <- dotplot(compKEGG_up, showCategory = 20, title = "KEGG Enrichment Analysis for Tet3 ko upregulated genes")
p10

DMRgo_down <- enrichGO(gene = na.omit(downgenes$entrez), OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

p11 <- dotplot(DMRgo_down,title="Tet3 ko Downregulated genes GO terms") 
p11
compKEGG_down <- compareCluster(geneCluster = downgenes , fun = "enrichKEGG", organism = "mouse", pvalueCutoff = 0.05, pAdjustMethod = "BH") 
p12 <- dotplot(compKEGG_down, showCategory = 20, title = "KEGG Pathway Enrichment in Downregulated Genes in tet3 KO NAc tissue")
p12



dat_go <- as.data.frame(DMRgo_up)
dat_keg <- as.data.frame(compKEGG_up)
view(dat_go)
view(dat_keg)


#prep normalized reads count

normed_counts$V1 <- rownames(normed_counts)
all_result <- merge(normed_counts,test,by="V1",all=F)

write.table(all_result, "Tet3_KO_RNAseq_results.txt",row.names = F)
