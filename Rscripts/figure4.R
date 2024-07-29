download_work_dir <- getwd()
setwd(download_work_dir)

### libraries
library(ggpubr)
library(ggnewscale)
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
library(sciClone)
library(eclipse)
### read in data
file.in = "./sourcedata/TableS3.clinical_sample_information.xls"
clin <- read.delim(file.in, stringsAsFactors = F, fileEncoding = "GBK")


uclin <- readxl::read_xlsx("./sourcedata/TableS1.Chracterization_features_of_CRC_patients.xlsx",skip=1)
uclin <- as.data.frame(uclin, stringsAsFactors = F)
 
### tissue snv for four patients with multiple sampling time points

### P45,  P46, P47, P56

### eclapse
file.in = "./sourcedata/TableS7.ECLIPSE_input.xls"
eclipse.input <- read.delim(file.in, stringsAsFactors = F, check.names = F)
### run ESCLIPSE https://github.com/amf71/ECLIPSE

eclipse.input <- eclipse.input[!is.na(eclipse.input$ddp),,drop=T]
eclipse.input <- as.data.table(eclipse.input)

hard_filters <- c( "primer_abundance_filter", "primer_strand_bias","forward_strand_bias", 'sequence_strand_bias', 'dro_cutoff', 'dao_imbalance' )
mrd_filters <- c( "tnc_error_rate" )

eclipse.input[, hard_filtered := grepl( paste( hard_filters, collapse = "|" ), `failed filters` ) | is.na(tnc_error_rate) | ddp == 0 ]
eclipse.input[, mrd_filtered := grepl( paste( mrd_filters, collapse = "|" ), `failed filters` ) | is.na(tnc_error_rate) | ddp == 0 ]


normsd <- extract_normalised_sd(eclipse.input, sample_id_col = 'sample_id', hard_filters = hard_filters, 
                                tumour_totcn_col = 'tumour_total_cpn', normal_totcn_col = 'normal_total_cpn')

eclipse.res <- clonal_deconvolution(data = eclipse.input, normalisedSD_max = normsd['hci95'], sample_id_col ='sample_id', niose_col= 'tnc_error_rate', 
                               chromosome_col = 'chromosome', position_col = 'position', alt_base_col = 'alternate', 
                               varcount_col = 'dao', depth_col = 'ddp', clone_col = 'PyCloneCluster_SC', 
                               is_clonal_col = 'is_clonal', multiplicity_col = 'mean_multiplicity',
                               tumour_totcn_col = 'tumour_total_cpn', normal_totcn_col = 'normal_total_cpn', 
                               hard_filtered_col = 'hard_filtered', mrd_filtered_col = 'mrd_filtered')
eclipse.res <- eclipse.res[order(eclipse.res$days_post_surgery),, drop=T]


### fishplot
library(fishplot)
library(clonevol)

for (patientID in sel.pts){
#  patientID <- 'P46'
#  tmp <- ctdna.tipc[ctdna.tipc$PID==patientID,,drop=T]
  tmp <- eclipse.res[grep(patientID, eclipse.res$sample_id),,drop=T]
  tmp$PID <- gsub("_.*","", tmp$sample_id)
    ### fishplot
  #provide a list of timepoints to plot
  #You may need to add interpolated points to end up with the desired
  #visualization. Example here was actually sampled at days 0 and 150
  
  timepoints = as.numeric(unique(tmp$days_post_surgery))

  tmp.mat <- reshape2::dcast(MutationName + PyCloneCluster_SC + Gene.Symbol ~ sample_id, data= tmp, value.var="clone_ccf", fun=sum)
  tmp.mat <- as.data.frame(tmp.mat, stringsAsFactors = F)
  
  vaf.col.names <- grep('^P[0-9][0-9]',  colnames(tmp.mat), value = T)

  sample.names <-  grep('^P[0-9][0-9]',  colnames(tmp.mat), value = T)
  
  
  # prepare sample grouping
#  sample.groups <- gsub("^P*.*_","", vaf.col.names)
#  sample.groups <- gsub(" ","_", sample.groups)
  sample.groups <- sample.names
  names(sample.groups) <- vaf.col.names
  
  
  # setup the order of clusters to display in various plots (later)

  tmp.mat <- tmp.mat[!is.na(tmp.mat$PyCloneCluster_SC),,drop=T]
  tmp.mat <- tmp.mat[tmp.mat$PyCloneCluster_SC !=".",,drop=T]
  tmp.mat <- tmp.mat[order(tmp.mat$PyCloneCluster_SC),]
  tmp.mat$PyCloneCluster_SC[which(tmp.mat$PyCloneCluster_SC==0)] <- 1
  
  tmp.mat[,c(4:ncol(tmp.mat))] <- apply(tmp.mat[,c(4:ncol(tmp.mat))],2, function(x){ifelse(is.na(as.numeric(x)), 0, x*100)} )
  
  #### https://gist.github.com/chrisamiller/f4eae5618ec2985e105d05e3032ae674
  #clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e',"#9370DB", "#E066FF")
  
  clone.colors <- c("#ece2f0","#a6bddb","#1c9099","grey88")
  
  #clone.colors <- NULL
  colnames(tmp.mat)[2] <- 'cluster'
  res = infer.clonal.models(variants = tmp.mat,
                            cluster.col.name = 'cluster',
                            ccf.col.names = vaf.col.names,
                            sample.groups = sample.groups,
                            cancer.initiation.model='polyclonal',#polyclonal monoclonal
                            subclonal.test = 'bootstrap',
                            subclonal.test.model = 'non-parametric',
                            num.boots = 1000,
                            founding.cluster = 1,
                            cluster.center = 'mean',
                            ignore.clusters = NULL,
                            clone.colors = clone.colors,
                            min.cluster.vaf = 0,
                            # min probability that CCF(clone) is non-negative
                            sum.p = 0.05,
                            # alpha level in confidence interval estimate for CCF(clone)
                            alpha = 0.05)
  
  # map driver events onto the trees
  res <- transfer.events.to.consensus.trees(res,
                                            tmp.mat,
                                            cluster.col.name = 'cluster',
                                            event.col.name = 'Gene.Symbol')
  
  # prepare branch-based trees
  res <- convert.consensus.tree.clone.to.branch(res, branch.scale = 'sqrt')
  y <- res
  
  #create a fish object
  #fish = createFishObject(frac.table,parents,timepoints=timepoints)
  
  ## create a list of fish objects - one for each model (in this case, there's only one)
  res = generateFishplotInputs(results=res)
  fishes = createFishPlotObjects(res)
  
  
  
  ## plot each of these with fishplot
  sample.groups <- gsub("^P*.*_","", vaf.col.names)
  outdir <- "./fishPlot/"
  pdf(paste0(outdir,patientID,'_fish.pdf'), width=8, height=4)
  for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,res$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm=patientID, cex.title=1,
#             vlines=seq(1, length(sample.names)), vlab=sample.names, 
             vlines = timepoints/max(timepoints)+1, vlab=as.numeric(sample.groups), 
             pad.left=0.5,
             bg.type = "gradient", bg.col = c("white","white","white"))
  }
  dev.off()

  ### https://github.com/hdng/clonevol/blob/master/examples/aml1.r
  plot.clonal.models(y,
                     # box plot parameters
                     box.plot = TRUE,
                     fancy.boxplot = TRUE,
                     fancy.variant.boxplot.highlight = 'is.driver',
                     fancy.variant.boxplot.highlight.shape = 21,
                     fancy.variant.boxplot.highlight.fill.color = 'red',
                     fancy.variant.boxplot.highlight.color = 'black',
                     fancy.variant.boxplot.highlight.note.col.name = 'gene',
                     fancy.variant.boxplot.highlight.note.color = 'blue',
                     fancy.variant.boxplot.highlight.note.size = 2,
                     fancy.variant.boxplot.jitter.alpha = 1,
                     fancy.variant.boxplot.jitter.center.color = 'grey50',
                     fancy.variant.boxplot.base_size = 12,
                     fancy.variant.boxplot.plot.margin = 1,
                     fancy.variant.boxplot.vaf.suffix = '.VAF',
                     # bell plot parameters
                     clone.shape = 'bell',
                     bell.event = TRUE,
                     bell.event.label.color = 'blue',
                     bell.event.label.angle = 60,
                     clone.time.step.scale = 1,
                     bell.curve.step = 2,
                     # node-based consensus tree parameters
                     merged.tree.plot = TRUE,
                     tree.node.label.split.character = NULL,
                     tree.node.shape = 'circle',
                     tree.node.size = 30,
                     tree.node.text.size = 0.5,
                     merged.tree.node.size.scale = 1.25,
                     merged.tree.node.text.size.scale = 2.5,
                     merged.tree.cell.frac.ci = FALSE,
                     # branch-based consensus tree parameters
                     merged.tree.clone.as.branch = TRUE,
                     mtcab.event.sep.char = ',',
                     mtcab.branch.text.size = 1,
                     mtcab.branch.width = 0.75,
                     mtcab.node.size = 3,
                     mtcab.node.label.size = 1,
                     mtcab.node.text.size = 1.5,
                     # cellular population parameters
                     cell.plot = TRUE,
                     num.cells = 100,
                     cell.border.size = 0.25,
                     cell.border.color = 'black',
                     clone.grouping = 'horizontal',
                     #meta-parameters
                     scale.monoclonal.cell.frac = TRUE,
                     show.score = FALSE,
                     cell.frac.ci = TRUE,
                     disable.cell.frac = FALSE,
                     # output figure parameters
                     out.dir = paste0('./fishPlot/',patientID),
                     out.format = 'pdf',
                     overwrite.output = TRUE,
                     width = 8,
                     height = 4,
                     # vector of width scales for each panel from left to right
                     panel.widths = c(3,4,2,4,2),
                     max.num.models.to.plot=3)
  # plot trees only  
  pdf(paste0(outdir,patientID,'_tree.pdf'), width=18, height=10)
  plot.all.trees.clone.as.branch(y, branch.width = 0.5,
                                 node.size = 1, node.label.size = 0.5)
  dev.off()
}

### dynamic change curves with time
eclipse.res.df <- aggregate(eclipse.res, clonal_purity_mut ~ sample_id, max)
eclipse.res.df <- tidyr::separate(eclipse.res.df, col="sample_id", into=c("PID","days_post_surgery"), sep="_")
p1 <- ggplot(eclipse.res.df, aes(days_post_surgery,clonal_purity_mut, group=PID )) +
        geom_point(color="darkred") +
        geom_line(color="grey") +
        facet_wrap(~ PID, nrow =1, scales = "free_x", strip.position = "top") +
        theme_pubr() +
        theme(strip.background = element_rect(colour="black", fill="white"), strip.placement = "outside")+
        xlab("Days post surgery") + ylab("Clonal ctDNA levels")
outdir <- './Figures/singleFigures/'
#ggsave(paste0(outdir,'figure4c_dynamic_ctDNA.pdf'), p1, width=8, height=2.5)


### TIO

file.in = './sourcedata/TableS5.ctDNA_mutations_TIO.xls'
ctdna.tio <- read.delim(file.in , stringsAsFactors = F)

### TIPC
file.in = './sourcedata/TableS4.ctDNA_mutations_TIPC.xls'
ctdna.tipc <- read.delim(file.in , stringsAsFactors = F)

### paired TIPC and TIO

ctdna.tio$flag <- paste(ctdna.tio$publicationID, ctdna.tio$MutationName)
ctdna.tipc$tag <- paste(ctdna.tipc$publicationID, ctdna.tipc$MutationName)
ctdna <- merge(ctdna.tipc, ctdna.tio[,c("flag","caseAF")], by.x="tag", by.y="flag", all=T, no.dups = T)
colnames(ctdna) <- gsub("caseAF.y","TIO_caseAF", colnames(ctdna))
colnames(ctdna) <- gsub("caseAF.x","TIPC_caseAF", colnames(ctdna))
ctdna.tio$publicationID <- gsub(" ","_",ctdna.tio$publicationID)
ctdna.tio$flag <- paste(ctdna.tio$publicationID, ctdna.tio$MutationName)
### add TIO
eclipse.res$flag <- paste(eclipse.res$sample_id, eclipse.res$MutationName)
eclipse.res.tio <- eclipse.res[eclipse.res$flag %in% ctdna.tio$flag,,drop=T]

eclipse.res.df.tipc <- aggregate(eclipse.res, clonal_purity_mut ~ sample_id, max)
eclipse.res.df.tipc$group <- "TIPC"

eclipse.res.df.tio <- aggregate(eclipse.res.tio, clonal_purity_mut ~ sample_id, max)
eclipse.res.df.tio$group <- "TIO"
eclipse.res.df.tio.add <- data.frame(sample_id = eclipse.res.df.tipc$sample_id[!eclipse.res.df.tipc$sample_id %in% eclipse.res.df.tio$sample_id],
                                     clonal_purity_mut = 0,
                                     group = "TIO")

eclipse.res.df <- rbind(eclipse.res.df.tipc, eclipse.res.df.tio, eclipse.res.df.tio.add)
eclipse.res.df <- tidyr::separate(eclipse.res.df, col="sample_id", into=c("PID","days_post_surgery"), sep="_")

p2 <- ggplot(eclipse.res.df, aes(days_post_surgery,clonal_purity_mut, group=PID )) +
  geom_point(aes(color=factor(group)) ) +
  geom_line(aes(group=group),color = "grey") +
  facet_wrap(~ PID, nrow =1, scales = "free_x", strip.position = "top") +
  theme_pubr() +
  theme(strip.background = element_rect(colour="black", fill="white"), strip.placement = "outside")+
  xlab("Days post surgery") + ylab("Clonal ctDNA levels") +
  scale_color_manual(values = c("#0072B2","darkred"))
 
outdir <- './Figures/singleFigures/'
ggsave(paste0(outdir,'figure4c_dynamic_ctDNA_v2.pdf'), p2, width=10, height=4.5)
