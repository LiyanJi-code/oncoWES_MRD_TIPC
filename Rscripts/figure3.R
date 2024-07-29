download_work_dir <- getwd()
setwd(download_work_dir)


library(atable)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(ggnewscale)
library(fst)
library(gtable)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library("cowplot")
### read in data
file.in = "./sourcedata/TableS3.clinical_sample_information.xls"
clin <- read.delim(file.in, stringsAsFactors = F)

uclin <- readxl::read_xlsx("./sourcedata/TableS1.Chracterization_features_of_CRC_patients.xlsx",skip=1)
uclin <- as.data.frame(uclin, stringsAsFactors = F)

### summary table

the_table <- atable::atable(uclin,
                            target_cols = c("Age", "Gender", "Stage"),
                            group_col = "TIPC_ctDNA")



### swimmers plot
label_size <- 3      
max_time <- ceiling(max(max(clin$relapse_to_surgery), max(clin$days_post_surgery) )/100) *100
p_swimmer <- ggplot(uclin, aes(x = relapse_to_surgery, y = PID)) + 
  geom_segment(aes(x = -5, xend = relapse_to_surgery, y = PID, yend = PID), size = 0.5, color="grey") + # linetype = 'dotted') +
  new_scale_colour() +
  geom_point(data = clin, aes(x = days_post_surgery, y = PID,  fill=`TIPC_ctDNA_status`, shape=`TIPC_ctDNA_status`), size = 3) +
  scale_fill_manual(name='TIPC_ctDNA_status', values=c("Positive"="darkred", "Negative"="#bababa"),na.value="#bababa") + 
      #new_scale_fill() +
  geom_point(data = uclin[uclin$Relapse %in% 1,], aes(x = relapse_to_surgery, y = PID), size = 3, shape = 'x') + 
  geom_point(data = uclin[uclin$Relapse %in% 0,], aes(x = relapse_to_surgery, y = PID), size = 3, shape = '>') + 
  scale_shape_manual(values=c(1, 21, 2, 6, 8))+     ### define the shape of geom_point 
  ggtitle("PATIENT ID") +
  ylab('Coloretal cancer patients') + xlab('Days post surgery') +
  scale_x_continuous(expand = c(0,0), limits = c(-5,max_time), position = 'bottom') +
  theme_classic() + theme(axis.line.y = element_blank(), 
                          strip.background = element_rect(fill = '#d9d9d9'), 
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                         # axis.title = element_blank(), 
                          plot.title.position = 'plot', 
                          plot.title = element_text(size = label_size + 1, vjust = -1)) 
p_swimmer
p_leng <- ggplot(uclin, aes(x = relapse_to_surgery, y = PID)) + 
  geom_point(aes(shape = factor(Relapse)), size = 2) + 
  scale_shape_manual(name = 'Relapse', values = c('1' = 'x', '0' = '>'),
                     labels = c('1' = 'Death/Recurrence', '0' = 'Censored')) +
  theme_bw()

p_ctdna <- ggplot(uclin, aes(x = 1, y = PID, fill = TIPC_ctDNA)) + 
  geom_tile(size = 0) + 
  scale_fill_manual(name = 'TIPC_ctDNA', values = c('Negative' = '#d1e5f04D', 'Positive' = 'grey33' ),  na.value="#bababa") + #'#e6f5d0'
  scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'TIPC_ctDNA') +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())

p_ctdna2 <- ggplot(uclin, aes(x = 1, y = PID, fill = TIO_ctDNA)) + 
  geom_tile(size = 0) + 
  scale_fill_manual(name = 'TIO_ctDNA', values = c('Negative' = '#d1e5f04D', 'Positive' = 'grey33' ),  na.value="#bababa") + #'#e6f5d0'
  scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'TIO_ctDNA') +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())




#Extract the legend from a ggplot object with this function:
g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) == 0){
    return(NULL)
  }
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:
align.legends.fun  <- function(legend.list){
  aligned.leg <- legend.list[[1]]
  for(i in 2:length(legend.list)){
    leg1        <- legend.list[[i]]$grobs[[1]]
    leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))
    leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)
    aligned.leg <- leg_final
  }
  return(aligned.leg)
}

#align legends
legends         <- list(
  g_legend(p_swimmer),
  g_legend(p_ctdna),  
  g_legend(p_ctdna2)   
)

aligend.legends <- align.legends.fun(legends)
aligend.legends1 <- align.legends.fun(legends[c(1:3)]) 

#exclude legends from original plots
p_swimmer      <- p_swimmer + theme(legend.position = 'none')
p_ctdna      <- p_ctdna + theme(legend.position = 'none')
p_ctdna2      <- p_ctdna2 + theme(legend.position = 'none')

#set margins
p_swimmer      <- ggplotGrob(p_swimmer + theme(plot.margin = unit(c(0.5, 0, 0.5, 0.1), "cm")))
p_ctdna        <- ggplotGrob(p_ctdna + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
p_ctdna2      <- ggplotGrob(p_ctdna2 + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))


plot_list <- list(p_swimmer, p_ctdna, p_ctdna2) 
all_heights <- lapply(plot_list, function(x) {x$heights})
plot_list_alignedHeights <- lapply(plot_list, function(x){
  x$heights <- do.call(unit.pmax, all_heights)
  return(x)
})



full_list <- c(rev(plot_list_alignedHeights))
widths    <- c(0.02, 0.02, 0.6) #  / 2
xpos      <- c(0, cumsum(widths)) + 0.01
full_plot <- cowplot::ggdraw()
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0, width = widths[x], height = 0.98)
}
plot(full_plot)
plot_list_legend <- list(plot = plot_list_alignedHeights, legend = aligend.legends)



pdf("./Figures/singleFigures/figure4a.pdf", width=18,height=9)
plot(full_plot)
dev.off()
pdf("./Figures/singleFigures/figure4a_legend.pdf", width=6,height=6)
plot(plot_list_legend$legend)
dev.off()

### stage and ctdna detection

cols <- c("#317ec2","#c03830", "#825ca6", "#c43d96",  "#e7872b", "#f2b342", "#5aaa46")
cols_light <- colorspace::lighten(cols, amount = 0.6, space = "HLS")
g1 <- ggplot(data = uclin,aes(x= Stage, fill = TIPC_ctDNA)) +
  geom_bar(position="fill", stat="count")+ 
  scale_y_continuous(labels = scales::percent) + 
  theme_pubr(legend="right") +
  scale_fill_manual(values=cols_light) +
  ylab("Percentage of patients") + 
  geom_text(
    aes(label = after_stat(
      scales::percent(
        ave(count, x, FUN = function(x) x / sum(x))
      )
    )),
    stat = "count", position = "fill"
  )

g2 <- ggplot(data = uclin,aes(x= Stage, fill = TIO_ctDNA)) +
  geom_bar(position="fill", stat="count")+ 
  scale_y_continuous(labels = scales::percent) + 
  theme_pubr(legend="right") +
  scale_fill_manual(values=cols_light) +
  ylab("Percentage of patients") +geom_text(
    aes(label = after_stat(
      scales::percent(
        ave(count, x, FUN = function(x) x / sum(x))
      )
    )),
    stat = "count", position = "fill"
  )

df2 <- as.data.frame(table(uclin$Stage), stringsAsFactors = F)
df2$value <- df2$Freq/sum(df2$Freq) *100

df2$labs <- paste0(df2$Var1, "\n(", round(df2$value,2), "%)")


fig4a <- ggdonutchart(data = df2, x = "value", label = "labs", color = "white",
                      lab.font = "white",  fill = "Var1",lab.pos = "in", 
                      orientation = "horizontal")+
#        ggsci::scale_fill_jama()+
  scale_fill_manual(values=cols_light) +
        theme(legend.position = "none")+
  annotate(geom = 'text', x = 0.5, y = 0, label = "n=59")

df3 <- reshape2::melt(uclin, id="Stage", measure=c("TIPC_ctDNA","TIO_ctDNA"))
df3.count <- as.data.frame(table(df3$Stage, df3$variable, df3$value))
df3.count <- df3.count[df3.count$Var3=="Positive",,drop=T]
#df3.count <- df3.count[df3.count$Var1!="I",,drop=T]
g3 <- ggbarplot(data = df3.count, x= "Var1", y="Freq", 
                color = "Var2", group="Var2", position = position_dodge(0.9)) +
                scale_fill_manual(values=cols_light) +
                stat_compare_means(method="kruskal.test") + 
                xlab("Stage")+ylab("Number of patients")

res <- ggdraw() +
  draw_plot(fig4a, x = .01, y = .7, width = .4, height = .3) +
  draw_plot(g1, x = .01, y = .35, width = .5, height = .3) +
  draw_plot(g2, x = .01, y = 0, width = .5, height = .3) +
  draw_plot_label(label = c("A","B", "C"), size = 12,
                  x = c(.01, .01, 0.01), y = c(1, .7, .35))

ggsave("./Figures/singleFigures/figure3a-c.pdf", res, width=9, height = 9)


### mutational spectrum
file.in = "./sourcedata/TableS2.Somatic_mutations_CRC_byWES.xls"
tissue.wes <- read.delim(file.in,stringsAsFactors = F)

cnv_events = c(c("Amp", "Del"), unique(as.character(tissue.wes$Function) ) )

oncomat <- reshape2::dcast(data = tissue.wes, formula = Gene.Symbol ~ PID,
                                 fun.aggregate = function(x, cnv = cnv_events){
                                   #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                   x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                   xad = x[x %in% cnv]
                                   xvc = x[!x %in% cnv]
                                   
                                   if(length(xvc)>0){
                                     xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                   }
                                   
                                   x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                   x = gsub(pattern = ';$', replacement = '', x = x)
                                   x = gsub(pattern = '^;', replacement = '', x = x)
                                   x = gsub(pattern = '"', replacement = '', x = x)
                                   return(x)
                                 } , value.var = 'Function', fill = '', drop = FALSE)

#convert to matrix
data.table::setDF(oncomat)
rownames(oncomat) = oncomat$Gene.Symbol
oncomat = as.matrix(oncomat[,-1, drop = FALSE])

variant.classes = as.character(unique(tissue.wes$Function))
variant.classes = c('',variant.classes, 'Multi_Hit')
names(variant.classes) = 0:(length(variant.classes)-1)

#Complex variant classes will be assigned a single integer.
vc.onc = unique(unlist(apply(oncomat, 2, unique)))
vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
variant.classes2 = c(variant.classes, vc.onc)

oncomat.copy <- oncomat
#Make a numeric coded matrix
for(i in 1:length(variant.classes2)){
  oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
}

#convert from character to numeric
mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
rownames(mdf) = rownames(oncomat.copy)
mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
  length(x[x != "0"])
}))
#Sort by total variants
mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
#colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
nMut = mdf[, ncol(mdf)]

mdf = mdf[, -ncol(mdf)]

mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
tmdf = t(mdf) #transposematrix
mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
mdf = mdf.temp.copy

#organise original character matrix into sorted matrix
oncomat.copy <- oncomat.copy[,colnames(mdf)]
oncomat.copy <- oncomat.copy[rownames(mdf),]
oncomat.copy[oncomat.copy=="."] <- ""

tmb <- as.data.frame(table(tissue.wes$PID))
clonal <- as.data.frame(table(tissue.wes$PID[tissue.wes$PYVI_cluster==1]))
subclonal <- as.data.frame(table(tissue.wes$PID[tissue.wes$PYVI_cluster!=1]))
tnb <- as.data.frame(table(tissue.wes$PID[tissue.wes$isTNB=="TNB"]))

### oncoplot

anno_df <- uclin[,c(9,3:5,7,11)]

rownames(anno_df) <- anno_df$PID
anno_df <- anno_df[,-1]

anno_df <- anno_df[colnames(oncomat.copy), ,drop=F]

ha <- HeatmapAnnotation(
  df = anno_df, col = list(
    Gender = c("Male" = "#542788", "Female"= "#d8daeb"),
    Age = circlize::colorRamp2(c(19,84), c( "white", "#313695")),
    Stage = c("I" = "#d1e5f0", "II"= "#92c5de","III" = "#2166ac", "IV"= "#053061"), 
    TIPC_ctDNA = c("Positive" = "#67001f", "Negative"= "grey66"),
    TIO_ctDNA = c("Positive" = "#7f3b08", "Negative"= "grey66")
  ))

cbPalette <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")            
cbPalette <- colorspace::lighten(cbPalette, amount = 0.2, space = "HLS")
names(cbPalette) <- c("missense","nonsense","frameshift","span","cds-del","splice-3","splice-5","cds-ins","cds-indel")
excludegene <- c("TTN", "MUC16", "SYNE1", "NEB", "MUC19", "CCDC168", "FSIP2", "OBSCN", "GPR98","USH2A")

### ctdna
file.in = './sourcedata/TableS4.ctDNA_mutations_TIPC.xls'
ctdna.tipc <- read.delim(file.in , stringsAsFactors = F)

ctdna.tipc$flag <- paste(ctdna.tipc$PID, ctdna.tipc$MutationName)

tissue.wes$flag <- paste(tissue.wes$PID, tissue.wes$MutationName)

ctdna.tipc$clonality <- NA
ctdna.tipc$clonality[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$PYVI_cluster!=1]] <- "tSubclonal"
ctdna.tipc$clonality[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$PYVI_cluster==1]] <- "tClonal"

ctdna.tipc$group <- NA
ctdna.tipc$group[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$PYVI_cluster!=1]] <- "tSubclonal"
ctdna.tipc$group[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$isTNB=="TNB"]] <- "tTNB"
ctdna.tipc$group[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$caseAF > 0.2]] <- "tVAFh"
ctdna.tipc$group[ctdna.tipc$flag %in% tissue.wes$flag[tissue.wes$PYVI_cluster==1]] <- "tClonal"


clonal.ctdna <- as.data.frame(table(ctdna.tipc$PID, ctdna.tipc$clonal))
clonal.ctdna.mat <- reshape2::dcast(clonal.ctdna, Var1 ~ Var2, value.var="Freq")
clonal.ctdna.mat$ctdna_clonal <- clonal.ctdna.mat$tClonal/(clonal.ctdna.mat$tClonal + clonal.ctdna.mat$tSubclonal)

ctdna.type <- as.data.frame(table(ctdna.tipc$PID, ctdna.tipc$group))
ctdna.type.mat <- reshape2::dcast(ctdna.type, Var1 ~ Var2, value.var="Freq")
#ctdna.type.mat$ctdna_clonal <- ctdna.type.mat$tClonal/sum(ctdna.type.mat$tClonal,ctdna.type.mat$tSubclonal,ctdna.type.mat$tTNB,ctdna.type.mat$tVAFh)

add.neg.pts <- data.frame(Var1 = uclin$PID[!uclin$PID %in% ctdna.type.mat$Var1],
                          tClonal = rep(0,13),
                          tSubclonal= rep(0,13),
                          tTNB = rep(0,13),
                          tVAFh = rep(0,13))
ctdna.type.mat.all <- rbind(ctdna.type.mat, add.neg.pts)
rownames(ctdna.type.mat.all) <- ctdna.type.mat.all$Var1
ctdna.type.mat.all <- ctdna.type.mat.all[,-1]

ctdna.type.mat.all <- ctdna.type.mat.all[colnames(oncomat.copy),,drop=T]
stat_ctdna <- ctdna.type.mat.all
stat_ctdna$clonal_rate <- stat_ctdna$tClonal/rowSums(stat_ctdna)
mean(stat_ctdna$clonal_rate,na.rm=T)

###heatmap
bar_anno <- merge(clonal, subclonal, by="Var1", all=T, no.dups=T)
bar_anno <- merge(bar_anno, clonal.ctdna.mat[,c("Var1","ctdna_clonal")] , by="Var1", all=T, no.dups =T)
bar_anno <- merge(bar_anno, tnb , by="Var1", all=T, no.dups =T)
colnames(bar_anno) <- c("Var1","clonal","subclonal","ctdna_clonal","TNB")
rownames(bar_anno) <- bar_anno$Var1
bar_anno <- bar_anno[,-1]
bar_anno <- bar_anno[colnames(oncomat.copy),,drop=T]
bar_anno$log_Clonal <- log2(bar_anno$clonal)
bar_anno$log_subClonal <- log2(bar_anno$subclonal)
bar_anno$log_tnb <- log2(bar_anno$TNB)

#col_fun = circlize::colorRamp2(c(0, 5, 10), c("#313695", "white", "#a50026"))
column_ha = HeatmapAnnotation(Subclonal = anno_barplot(bar_anno$log_subClonal, gp=gpar(fill = "#CCCCCC80", col="#CCCCCC80") ) ,
                              Clonal = anno_barplot(bar_anno$log_Clonal, gp=gpar(fill = "#a50026", col="#CCCCCC80") ), 
                              TNB = anno_barplot(bar_anno$log_tnb, gp=gpar(fill = "#313695", col="#CCCCCC80")  ),
                              ctdna_Clonal=anno_barplot(as.matrix(ctdna.type.mat.all), gp=gpar(fill = cbPalette, col= cbPalette)  )
                              )

oncomat.dat <- oncomat.copy[!rownames(oncomat.copy) %in% excludegene,]

pdf("./Figures/singleFigures/figure3.pdf", width=8,height=10)

Heatmap(oncomat.dat[1:20,],col = cbPalette, na_col = "grey99",
                        show_column_names = F,row_names_side = "left",
                        top_annotation = ha,
                        bottom_annotation = column_ha
        )
dev.off()

### TIO

file.in = './sourcedata/TableS5.ctDNA_mutations_TIO.xls'
ctdna.tio <- read.delim(file.in , stringsAsFactors = F)

### paired TIPC and TIO

ctdna.tio$flag <- paste(ctdna.tio$publicationID, ctdna.tio$MutationName)
ctdna.tipc$tag <- paste(ctdna.tipc$publicationID, ctdna.tipc$MutationName)
ctdna <- merge(ctdna.tipc, ctdna.tio[,c("flag","caseAF")], by.x="tag", by.y="flag", all=T, no.dups = T)
colnames(ctdna) <- gsub("caseAF.y","TIO_caseAF", colnames(ctdna))
colnames(ctdna) <- gsub("caseAF.x","TIPC_caseAF", colnames(ctdna))

ctdna.df2 <- reshape2::melt(ctdna, id="publicationID", measure=c( "TIPC_caseAF", "TIO_caseAF"))
ctdna.df2$value <- as.numeric(ctdna.df2$value)
ctdna.df2$variable <- gsub("_caseAF","", ctdna.df2$variable)
### 
p1 <- ggpaired(ctdna.df2,  x = "variable", y = "value",
                         color = "variable", line.color = "gray", line.size = 0.4,
                         palette = "npg", width = 0.4, xlab="",
                      ylab="Variant allele frequency", legend="none") +
        stat_compare_means(
        method="wilcox",aes(label=paste0("p = ",..p.format..)), label.x=1.5) 

ctdna.df3_1 <- as.data.frame(table(ctdna$publicationID[!is.na(ctdna$TIPC_caseAF)]))
ctdna.df3_1$group <- "TIPC"
ctdna.df3_2 <- as.data.frame(table(ctdna$publicationID[!is.na(ctdna$TIO_caseAF)]))
ctdna.df3_2$group <- "TIO"
ctdna.df3 <- rbind(ctdna.df3_1, ctdna.df3_2)


p2 <- ggboxplot(ctdna.df3,  x = "group", y = "Freq",
                add="jitter",
               color = "group", line.color = "gray", line.size = 0.4,
               palette = "npg", width = 0.4, xlab="",
               ylab="Number of SNVs", legend="none") +  
    scale_y_continuous(trans = 'log10', limits = c(1,100),
                     breaks = c(1,2.5,5,10,25,50,100),
                     labels = c(1,2.5,5,10,25,50,100))+
    stat_compare_means(
    method="wilcox",aes(label=paste0("p = ",..p.format..)), label.x=1.5) 
res <- ggarrange(p1, p2, align = "hv")
outdir <- './Figures/singleFigures/'
ggsave(paste0(outdir,'figureS2a-b.pdf'), res, width=8, height= 4 )


