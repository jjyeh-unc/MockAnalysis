DIR = "~/Google Drive File Stream/My Drive/postdoc/JenJenYeh/analysis/MockAnalysis/"
RES = "~/Google Drive File Stream/My Drive/postdoc/JenJenYeh/analysis/RESULTS/"

library(M3C)
library(gplots)
library(NMF)
library(DESeq2)

options(stringsAsFactors = FALSE)

# input and pre-process data
load(paste0(DIR,"TCGA_PAAD_plus.RData"))
sample = as.character(TCGA_PAAD_plus$sampInfo$Tumor.Sample.ID[which(TCGA_PAAD_plus$sampInfo$Grade !="")])

Moffitt_basal25 <- read.table(paste0(DIR,"Moffitt-basal25.txt"), header = F)
Moffitt_classical25 <- read.table(paste0(DIR,"Moffitt-classical25.txt"), header = F)
AURKA_sig_gene <- read.table(paste0(DIR,"AURKA_sig_genes.txt"), header = F)
TCGA_PAAD_plus.cnt <- read.table(paste0(DIR,"TCGA_PAAD_plus.cnt.txt"), header = F, sep = "\t", row.names = 1)

#########
# Q1.
#########
data.Q1 = TCGA_PAAD_plus$ex[which(as.character(TCGA_PAAD_plus$featInfo$SYMBOL) %in% c(Moffitt_basal25[,1],Moffitt_classical25[,1])),sample]
data.Q1 = log2(data.Q1+1)
rownames(data.Q1) = as.character(TCGA_PAAD_plus$featInfo$SYMBOL)[which(as.character(TCGA_PAAD_plus$featInfo$SYMBOL) %in% c(Moffitt_basal25[,1],Moffitt_classical25[,1]))]

# run consensus clustering
res.Q1 <- M3C(data.Q1, cores = 4, seed = 123, maxK = 2, removeplots = TRUE)
saveRDS(res.Q1, file = paste0(RES,"res_Q1.rds"))

res.Q1 <- readRDS(file = paste0(RES,"res_Q1.rds"))
annon = res.Q1$realdataresults[[2]]$ordered_annotation 
colnames(annon) = "subtype"

# not sure whether colnames changed, need manually check
rownames(annon) = gsub("\\.","-",rownames(annon))
rownames(annon)[which(!rownames(annon) %in% colnames(data.Q1))]
colnames(data.Q1)[which(!colnames(data.Q1) %in% rownames(annon))]
rownames(annon)[which(!rownames(annon) %in% colnames(data.Q1))] = c("TCGA-XN-A8T5-01A", "TCGA-IB-A7LX-01A", "TCGA-LB-A7SX-01A", "TCGA-XD-AAUI-01A", "TCGA-3A-A9IX-01A",
                                                                    "TCGA-XD-AAUH-01A", "TCGA-F2-A7TX-01A", "TCGA-XN-A8T3-01A", "TCGA-XD-AAUG-01A", "TCGA-HV-AA8X-01A",
                                                                    "TCGA-XD-AAUL-01A")
all(rownames(annon) %in% colnames(data.Q1))

# plot heatmap
data.Q1.order <- res.Q1$realdataresults[[2]]$ordered_data 
colnames(data.Q1.order) = rownames(annon)
data.Q1.order = t(scale(t(data.Q1.order))) 
data.Q1.order = apply(data.Q1.order, 2, function(x) ifelse(x > 4, 4, x)) 
data.Q1.order = apply(data.Q1.order, 2, function(x) ifelse(x < -4, -4, x)) 

pdf(file = paste0(RES,"Q1.pdf"))
ann_colors = rev(colorRampPalette(c('orange','blue'))(2))
pal = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
aheatmap(data.Q1.order, annCol = annon, Colv = NA, cexCol = 0.4, distfun = 'pearson', 
              color = gplots::bluered(256), annColors = list(subtype=ann_colors))
dev.off()


#########
# Q2.
#########
data.Q2 = TCGA_PAAD_plus$ex[which(as.character(TCGA_PAAD_plus$featInfo$SYMBOL) %in% AURKA_sig_gene[,1]),sample]
data.Q2 = log2(data.Q2+1)
rownames(data.Q2) = as.character(TCGA_PAAD_plus$featInfo$SYMBOL)[which(as.character(TCGA_PAAD_plus$featInfo$SYMBOL) %in% AURKA_sig_gene[,1])]

sample.classical = rownames(annon)[which(as.character(annon$subtype) == 1)]
sample.basal = rownames(annon)[which(as.character(annon$subtype) == 2)]

# t-test
p.val <- NULL
for (i in 1:nrow(data.Q2)) {
  p.val[i] = t.test(data.Q2[i,which(colnames(data.Q2) %in% sample.basal)],
                    data.Q2[i,which(colnames(data.Q2) %in% sample.classical)])$p.value
}
p.adj = round(p.adjust(p.val, method = "BH"),2)

# boxplot
pdf(file = paste0(RES,"Q2.pdf"))
par(mfrow = c(2,2))
for (i in 1:nrow(data.Q2)) {
  boxplot(data.Q2[i,which(colnames(data.Q2) %in% sample.basal)],
          data.Q2[i,which(colnames(data.Q2) %in% sample.classical)], 
          col = c("orange","blue"), names = c("basal-like","classical"),
          ylab = "expression", main = rownames(data.Q2)[i])
  legend("topright", inset =.02, 
         legend = paste0("p.adj=",p.adj[i]), cex=0.8)
}
dev.off()


#########
# Q3.
#########
data.Q3 = TCGA_PAAD_plus.cnt[c(3:20533), sapply(lapply(sample, function(x) grep(x,TCGA_PAAD_plus.cnt[1,])), function(x) x[1])]
data.Q3 = round(mapply(data.Q3, FUN=as.numeric))
colnames(data.Q3) = sample
rownames(data.Q3) = as.character(TCGA_PAAD_plus$featInfo$SYMBOL)
  
coldata = annon
# match data names
all(rownames(coldata) %in% colnames(data.Q3))
all(rownames(coldata) == colnames(data.Q3))
data.Q3 <- data.Q3[, rownames(coldata)]
all(rownames(coldata) == colnames(data.Q3))

# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = data.Q3,
                              colData = coldata,
                              design= ~ subtype)
dds = DESeq(dds)
resultsNames(dds) 
res.Q3 = results(dds, name="subtype_2_vs_1", lfcThreshold=1, alpha = 0.05)
saveRDS(res.Q3, file = paste0(RES,"res_Q3.rds"))

res.Q3 <- readRDS(file = paste0(RES,"res_Q3.rds"))
pdf(file = paste0(RES,"Q3.pdf"))
# MA plot
plotMA(res.Q3, ylim=c(-2,2), alpha = 0.05, main = "MA plot")

# volcano plot
cols = densCols(res.Q3$log2FoldChange, -log10(res.Q3$pvalue))
plot(res.Q3$log2FoldChange, -log10(res.Q3$padj), col = cols, panel.first = grid(),
     main = "Volcano plot", xlab = "Effect size: log2(fold-change)", ylab = "-log10(adjusted p-value)",
     pch = 20, cex = 0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(0.05), col="brown")

gn.selected = res.Q3$padj < 0.05 
text(res.Q3$log2FoldChange[gn.selected], -log10(res.Q3$padj)[gn.selected],
     lab = rownames(res.Q3)[gn.selected ], cex = 0.6)
dev.off()


#########
# Q4.
#########
data.Q4 = data.frame(sample = TCGA_PAAD_plus$sampInfo$Tumor.Sample.ID,
                     age = TCGA_PAAD_plus$sampInfo$Age.at.initial.pathologic.diagnosis,
                     gender = TCGA_PAAD_plus$sampInfo$Gender,
                     race = TCGA_PAAD_plus$sampInfo$Race)
data.Q4 = data.Q4[which(data.Q4$sample %in% sample),]
rownames(data.Q4) = data.Q4$sample

data.Q4 = merge(data.Q4, annon, by=0)

# chi-sequre test for subtype vs. gender
p.gender = chisq.test(data.Q4$subtype, data.Q4$gender)$p.value

# Fisher's exact test for subtype vs. race
p.race = fisher.test(data.Q4$subtype, data.Q4$race)$p.value

# Kruskal-Wallis test for subtype vs. age
p.age = kruskal.test(age~subtype, data = data.Q4)$p.value

print(data.frame(geneder = p.gender,
                 race = p.race,
                 age = p.age, row.names = "p-value"))


