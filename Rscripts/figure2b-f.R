download_work_dir <- getwd()
setwd(download_work_dir)

### load packages
library(ggpubr)
library(scales)

file.in <- "./WES-MRD文章原始数据/分析性能验证结果.xlsx"
dat1 <- readxl::read_xlsx(file.in, sheet="60ng-样本维度灵敏度特异性-与菁良MRD标准品性能一致" )

spec <- data.frame(Number_of_variants=as.numeric(dat1[1,-1]),
                   specificity=as.numeric(dat1[10,-1])
                   )

spec2 <- dat1[c(1,11:12),]
colnames(spec2) <- spec2[1,]
spec2 <- spec2[-1,]
colnames(spec2)[1] <- "Group"
spec.df <- reshape2::melt(spec2, id="Group")
spec.df$type <- gsub("20标准品","NA12878_standards", spec.df$Group)
spec.df$type <- gsub("44健康人","The_healthy_cfDNA", spec.df$type)
cols <- c("#fee090","#91bfdb")
f2e <- ggplot(spec.df, aes(variable, value, group="type") ) +
  geom_point(aes(colour = factor(type)), size=4) +
  geom_smooth(aes(group = type), method=loess,  se=FALSE)+ 
  # geom_smooth(aes(x = variable, y = value, group=type), method = "nls", se = FALSE, 
  #                    method.args = list(formula = y ~ (1/(1 + exp(-b*x + c))), start = list(b= 100, c=0))) +
  theme_pubr(legend="none") + #border() +
  geom_hline(yintercept=95, linetype="dashed", color = "grey")  +
  xlab("Monitoring variants") + ylab("Specificity (%)") +ylim(0,100) +
  scale_color_manual(values=cols) +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  ggbreak::scale_y_break(c(70, 75), scales = 0.2)


dat1 <- dat1[1:8,]
colnames(dat1) <- dat1[1,]
dat1 <- dat1[-1,]
colnames(dat1)[1] <- "VAF"
dat.df <- reshape2::melt(dat1, id="VAF")
colnames(dat.df)[1] <- "VAF(%)"


cols <- c(
  "#feedde",
  "#fdd0a2",
  "#fdae6b",
  "#fd8d3c",
  "#f16913",
  "#d94801",
  "#8c2d04")

f2f <- ggplot(dat.df, aes(variable, value, color = `VAF(%)`, group=`VAF(%)`) ) +
#  geom_point(size=4) +
#  geom_smooth(method=lm, formula = y ~ poly(x,3),  se=FALSE)+ 
#  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE) +
  geom_smooth(method=loess,  se=FALSE)+ 
  theme_pubr(legend="right") + #border() +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  geom_hline(yintercept=95, linetype="dashed", color = "grey")  +
#  geom_vline(xintercept="50", linetype="dashed", color = "grey")  +
  scale_color_manual(values=cols) + xlab("Number of variants") + ylab("Sensitivity (%)")
  


dat2 <- readxl::read_xlsx(file.in, sheet="ctDNA含量一致性(50)" )
colnames(dat2) <- c("Theoretical_ctDNA_content","Actural_ctDNA_content")
f2d <- ggscatter(dat2, x="Theoretical_ctDNA_content",y="Actural_ctDNA_content",
          color="grey", add="reg.line", add.params = list(color="#8c2d04"),
          conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab="Expected ctDNA fraction", ylab="Estimated ctDNA fraction"
          ) + 
              scale_y_continuous(labels = percent) + 
              scale_x_continuous(labels = percent) #+ border()

### 
input_depth <- readxl::read_xlsx("./WES-MRD文章原始数据/不同投入量下的有效深度统计.xlsx",sheet="Sheet1")
input_depth$id <- rownames(input_depth)
input_depth.df <- reshape2::melt(input_depth, id ="id")
input_depth.df <- input_depth.df[!is.na(input_depth.df$value),,drop=T]
input_depth.df$variable <- gsub(" ng","", input_depth.df$variable)
library(ggbeeswarm)
#compare to jitter

cbPalette <- c("#a6cee3","#386cb0","#fb9a99","#fdb462","#ef3b2c","#662506","#7fc97f","#984ea3","#ffff33")
                       
cbPalette <- colorspace::lighten(cbPalette, amount = 0.2, space = "HLS")                             
f2b <- ggplot(input_depth.df,aes(variable, value, color=factor(variable), group="variable")) + 
        geom_beeswarm(cex = 2.5, method = "swarm", aes(colour = factor(variable))) + theme_pubr(legend="none") + xlab("DNA input (ng)")+
        ylab("Effective depth") + 
        scale_color_manual(values=cbPalette) +
  stat_summary(fun.min=function(y) {mean(y) - qt((1-0.95)/2, length(y) - 1) * sd(y) / sqrt(length(y) - 1)}, fun.max=function(y) {mean(y) + qt((1-0.95)/2, length(y) - 1) * sd(y) / sqrt(length(y) - 1)}, geom='errorbar', width=0.3, size =1) +
  stat_summary(fun.min=mean, fun.max=mean, geom='errorbar', width=0.6, size =1) +
  geom_hline(yintercept=3000, linetype="dashed", color = "grey")  +
  ggbreak::scale_y_break(c(2800, 3500), scales = 6) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE, col="grey")


#https://rforbiochemists.blogspot.com/2015/06/plotting-two-enzyme-plots-with-ggplot.html

nonhotLoD <- readxl::read_xlsx("./WES-MRD文章原始数据/不同投入量下的有效深度统计.xlsx",sheet="Sheet2")
nonhotLoD <- nonhotLoD[,c(1,6)]
nonhotLoD$VAF_100 <- 100*nonhotLoD$VAF
nonhotLoD$specificity_100 <- 100*nonhotLoD$specificity
f2c <- ggplot(nonhotLoD, aes(VAF_100, specificity_100 ,color="#8c2d04") ) +
  geom_point(size=4) +
  geom_smooth(method=loess,  se=FALSE)+
  theme_pubr(legend="none") + #border() +
  geom_hline(yintercept=95, linetype="dashed", color = "grey")  +
  xlab("VAF (%)") + ylab("Sensitivity (%)") +ylim(0,100)


library("cowplot")
res <- ggdraw() +
  draw_plot(f2b, x = .39, y = .5, width = .28, height = .5) +
  draw_plot(f2c, x = .68, y = .5, width = .27, height = .5) +
  draw_plot(f2d, x = .01, y = 0, width = .4, height = .5) +
  draw_plot(f2e, x = .40, y = 0, width = .27, height = .5) +
  draw_plot(f2f, x = .68, y = 0, width = .35, height = .5) +
  draw_plot_label(label = c("B", "C","D", "E", "F"), size = 15,
                  x = c(.39, .67, 0, 0.39, 0.67), y = c(1, 1, .5, .5, .5))

ggsave("./Figures/singleFigures/figure2b-f.V2.pdf", res, width=18, height = 9)
