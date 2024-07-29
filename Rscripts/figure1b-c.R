download_work_dir <- getwd()
setwd(download_work_dir)

cap_eff <- readxl::read_xlsx("./WES-MRD文章原始数据/不同大小核心panel捕获效率.xlsx")
cap_eff$id <- rownames(cap_eff)
cap_eff <- reshape2::melt(cap_eff, id="id")
cap_eff <- cap_eff[!is.na(cap_eff$value),,drop=T]
cap_eff$panel_size <- gsub("K","",cap_eff$variable)
cap_eff$panel_size <- as.numeric(cap_eff$panel_size)
cap_eff$variable <- paste0(cap_eff$variable,"b")
cap_eff$variable <- factor(cap_eff$variable, levels = c("6Kb","12Kb","20Kb","22Kb","24Kb","30Kb"))
 

p4 <- ggscatter(cap_eff, x = "panel_size", y = "value", 
                add.params = list(color = "#8c2d04", fill = "grey88"), 
                color="variable", palette=cols_light,legend="right") + 
                ylab("Capture efficiency (%)") + xlab("Panel size (Kb)") +
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE, col="grey") +
    geom_vline(xintercept=26, linetype="dashed", color = "grey")  


outdir <- './Figures/singleFigures/'
ggsave(paste0(outdir,'figure1b.pdf'), p4, width=12, height= 6 )


healthy <- readxl::read_xlsx('./WES-MRD文章原始数据/健康阴性标准品突变检出结果.xlsx',sheet="SNV")
fp <- healthy[healthy$CRI == "H0",,drop = T]
fp$group <- gsub("阴性标准品1","NA12878_TIPC",fp$Mark...5 )
fp$group <- gsub("阴性标准品","SW480",fp$group )
fp$group <- gsub("健康样本","Healthy_cfDNA_TIPC",fp$group )
fp$Mutation <- paste0(fp$Chr...25, ":",fp$vcfPos...26)
fp$background_error_rate <- fp$caseAltDep/fp$caseDep

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

fp2 <- data_summary(fp, varname="background_error_rate", 
                    groupnames=c("SampleNo", "group"))



fp2$group <- factor(fp2$group, levels = c('SW480','NA12878_TIPC','Healthy_cfDNA_TIPC'))
fp2 <- fp2[order(fp2$group),,drop=T]
fp2$SampleNo <- factor(fp2$SampleNo,levels=unique(fp2$SampleNo) )
###
library(ggpubr)
library(scales)

p1 <- ggplot(fp2, aes(x = SampleNo, y = background_error_rate, color = '#999999')) +
  geom_point(size=3,color='#999999') + 
  scale_y_continuous(trans = 'log10', limits = c(0.00001,0.001),
                     breaks = c(0.00001,0.0001,0.001),
                     labels = c(0.00001,0.0001,0.001)) +
  geom_errorbar(aes(ymin=background_error_rate-sd, ymax=background_error_rate+sd), width=.2,
                position=position_dodge(0.05)) +
  theme_pubr(legend="none") + 
  ylab("Background error rate") + xlab("") +
  scale_color_manual(values=c("#8c2d04")) +
#  theme(axis.text.x = element_text(angle=60,hjust=1,size = 2)) 
  theme(axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank())


cols_light <- c("#AACCE9", "#CDBEDB", "#E7B1D5", "#E9ACA9", "#F5CFAA", "#FAE1B3", "#BBE0B3")

p1 <- p1 + annotation_logticks(side="l", outside = TRUE )  + coord_cartesian(clip = "off")  
pc <- p1 + geom_rect(aes(xmin=0, xmax=19, ymin=0, ymax=0.001), alpha=0.3, fill="#AACCE9",color="#AACCE9") +
  geom_rect(aes(xmin=20, xmax=28, ymin=0, ymax=0.001), alpha=0.3, fill="#CDBEDB",color="#CDBEDB") + 
  geom_rect(aes(xmin=29, xmax=36, ymin=0, ymax=0.001), alpha=0.3, fill="#E7B1D5",color="#E7B1D5")

outdir <- './Figures/singleFigures/'
ggsave(paste0(outdir,'figure1c.pdf'), pc, width=6, height=2.5)


