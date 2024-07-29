download_work_dir <- getwd()
setwd(download_work_dir)


uclin <- readxl::read_xlsx("./sourcedata/TableS1.Chracterization_features_of_CRC_patients.xlsx",skip=1)
uclin <- as.data.frame(uclin, stringsAsFactors = F)
### survival curve
library(survival)
library(survminer)
cols <- c("#317ec2","#c03830", "#825ca6", "#c43d96",  "#e7872b", "#f2b342", "#5aaa46")
cols_light <- colorspace::lighten(cols, amount = 0.6, space = "HLS")
xtitle = "Time in days"
ytitle = "Relapse free survival"

pl <- list()
fit <- survfit(Surv(relapse_to_surgery, Relapse) ~ TIPC_ctDNA, data = uclin)

column = "TIPC"
legend.labs <- names(fit$strata)
legend.labs = gsub(paste0(column,"_ctDNA="),"",legend.labs)

data.survdiff <- survival::survdiff(Surv(relapse_to_surgery, Relapse)~ TIPC_ctDNA, data = uclin)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("HR = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = "-"), sep = "")

pl[[column]] <- survminer::ggsurvplot(
  fit, 
  data = uclin, 
  size = 1,                 
  break.time.by = 180,	  # x轴间隔天数/月数
  xlab = xtitle,   # customize X axis label.
  ylab = ytitle,
  conf.int = F,          # Add confidence interval
  palette = cols_light,
   pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                              paste("p = ",round(p.val,3), sep = "")),
                HR, CI, sep = "\n") ,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.title = "No. at risk",
  tables.theme = theme_cleantable(), # clean theme for tables
  tables.y.text = FALSE, # hide tables y axis text 
  tables.height = 0.3, 
  tables.y.text.col = T, 
  tables.col = "strata",
  legend.labs = legend.labs,    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  title = column,
  legend.title="",
  ggtheme = theme_pubr()      # Change ggplot2 theme
)


fit <- survfit(Surv(relapse_to_surgery, Relapse)~ TIO_ctDNA, data = uclin)

column = "TIO"

legend.labs <- names(fit$strata)
legend.labs = gsub(paste0(column,"_ctDNA="),"",legend.labs)

data.survdiff <- survival::survdiff(Surv(relapse_to_surgery, Relapse)~ TIO_ctDNA, data = uclin)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("HR = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = "-"), sep = "")

pl[[column]] <- survminer::ggsurvplot(
  fit, 
  data = uclin, 
  size = 1,                 
  break.time.by = 180,	  # x轴间隔天数/月数
  xlab = xtitle,   # customize X axis label.
  ylab = ytitle,
  conf.int = F,          # Add confidence interval
  palette = cols_light,
  pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                             paste("p = ",round(p.val,3), sep = "")),
               HR, CI, sep = "\n") ,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.title = "No. at risk",
  tables.theme = theme_cleantable(), # clean theme for tables
  tables.y.text = FALSE, # hide tables y axis text 
  tables.height = 0.3, 
  tables.y.text.col = T, 
  tables.col = "strata",
  legend.labs = legend.labs,    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  title = column,
  legend.title="",
  ggtheme = theme_pubr()      # Change ggplot2 theme
)

res <- arrange_ggsurvplots(pl, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.2)
ggsave("./Figures/singleFigures/figure5a.pdf", res, width=12, height = 4.5)

### stratified by stage
model <- coxph( Surv(relapse_to_surgery, Relapse)~  Age + Gender + Stage  + TIPC_ctDNA , data = uclin )
ggforest(model)
model2 <- coxph( Surv(relapse_to_surgery, Relapse)~  Age + Gender + Stage  + TIO_ctDNA , data = uclin )
ggforest(model2)
pdf("./Figures/singleFigures/figure5b.pdf", width=5, height = 4.5)
ggforest(model)
dev.off()
pdf("./Figures/singleFigures/figureS2.pdf", width=5, height = 4.5)
ggforest(model2)
dev.off()
### sensitivity and specificity


calculate_statistics <- function(tbl) {
  P <- sum(tbl[,2])
  N <- sum(tbl[,1])
  TP <- tbl[2,2]
  TN <- tbl[1,1]
  FP <- tbl[1,2]
  FN <- tbl[2,1]
  TPR <- TP / P
  FPR <- FP / N
  FDR <- FP / (TP + FP)
  NPV <- TN / (TN + FN)
  data.frame(P = P, N = N, TP = TP, TN = TN, FP = FP, 
             FN = FN, PPV = TPR, NPV = NPV, FPR = FPR, FDR = FDR)
}
tbl <- table(uclin$TIPC_ctDNA,uclin$Relapse)
tipc_accurary <- calculate_statistics(tbl)
tbl <- table(uclin$TIO_ctDNA,  uclin$Relapse)
tio_accurary <- calculate_statistics(tbl)

tipc_accurary <- as.data.frame(t(tipc_accurary), stringsAsFactors = F)
colnames(tipc_accurary) <- "TIPC"
tio_accurary <- as.data.frame(t(tio_accurary), stringsAsFactors = F)
colnames(tio_accurary) <- "TIO"

acc <- cbind(tipc_accurary, tio_accurary)

plot.acc <- acc[rownames(acc) %in% c("PPV","NPV"),,drop=T]

plot.acc$group <- rownames(plot.acc)

plot.acc.df <- reshape2::melt(plot.acc, id="group")
plot.acc.df$variable <- factor(plot.acc.df$variable, levels = c("TIO", "TIPC"))

###
validate.pts <- c('P55','P56','P57','P58','P59')
uclin2 <- uclin[!uclin$PID %in% validate.pts,,drop=T]
tbl <- table(uclin2$TIPC_ctDNA,uclin2$Relapse)
tipc_accurary <- calculate_statistics(tbl)
tbl <- table(uclin2$TIO_ctDNA,  uclin2$Relapse)
tio_accurary <- calculate_statistics(tbl)

tipc_accurary <- as.data.frame(t(tipc_accurary), stringsAsFactors = F)
colnames(tipc_accurary) <- "TIPC"
tio_accurary <- as.data.frame(t(tio_accurary), stringsAsFactors = F)
colnames(tio_accurary) <- "TIO"

acc <- cbind(tipc_accurary, tio_accurary)

plot.acc <- acc[rownames(acc) %in% c("PPV","NPV"),,drop=T]

plot.acc$group <- rownames(plot.acc)

plot.acc.df <- reshape2::melt(plot.acc, id="group")
plot.acc.df$variable <- factor(plot.acc.df$variable, levels = c("TIO", "TIPC"))
fig4c <- ggbarplot(plot.acc.df, x = "group", y="value", 
                    fill="variable", palette = cols_light, label = T,
                    position = position_dodge(.5), width = .4,
                    xlab="", ylab="Ratio of patients",legend="top")
###

uclin3 <- uclin[uclin$PID %in% validate.pts,,drop=T]
tbl3 <- table(uclin3$TIPC_ctDNA,uclin3$Relapse)

### plot figure 5d
file.in <- "./sourcedata/TableS6.WESMRDpostivie_TIOnegative_n9.xlsx"
dat2 <- readxl::read_xlsx(file.in, sheet="Sheet1",skip=4)
colnames(dat2) <- c('Cohort','monitoringNo_TIPC','TIPC_MRD_status','TIPC_ctDNA_fraction','monitoringNo_TIO','TIO_MRD_status','Reason')
dat2$TIO <- 0
dat.df <- reshape2::melt(dat2, id="PID",measure.vars = c("TIPC_ctDNA_fraction","TIO") )
dat.df$variable <- factor(dat.df$variable, levels = c("TIO","TIPC_ctDNA_fraction"))
fig4d <- ggscatter(dat.df, x="PID", y="value", fill="variable", color="variable",
          palette = cols,xlab="", ylab="Estimated ctDNA fraction",
#          title="TIPC in TIO negative patients",
          legend="top")

### TIO negative muts and TIPC positive muts 

tioneg_tipcpos <- uclin[uclin$TIPC_ctDNA=="Positive" & uclin$TIO_ctDNA=="Negative",,drop=T]

tioneg_tipcpos$group <- "Off-target by TIO"
tioneg_tipcpos$group[tioneg_tipcpos$PID %in% dat2$PID[grep("Undetectable", dat2$Reason)]] <- "Under limitation"

df2 <- as.data.frame(table(tioneg_tipcpos$group), stringsAsFactors = F)
df2$value <- df2$Freq/sum(df2$Freq) *100

df2$labs <- paste0(df2$Var1, " (", round(df2$value,2), "%)")

### comparison in TIO and TIPC TIPC positive but TIO negative mutation: reason (off-target; low frequency) venn plot figure 5e
library(ggpubr)
fig4e <- ggdonutchart(data = df2, x = "value", label = "labs", color = "white",
                     lab.font = "white",  fill = "Var1",lab.pos = "in", 
                     orientation = "horizontal")+
#                  ggsci::scale_fill_jama()+
  scale_fill_manual(values=cols_light) +
                  theme(legend.position = "none")


library("cowplot")
fig5 <-  ggdraw() +
  draw_plot(fig4c, x = .01, y = 0, width = .33, height = .5) +
  draw_plot(fig4d, x = .34, y = 0, width = .37, height = .5) +
  draw_plot(fig4e, x = .71, y = 0, width = .33, height = .5) +
  draw_plot_label(label = c("C","D","E"), size = 15,
                  x = c(0, 0.33, 0.73), y = c( .5, .5, .5))

ggsave("./Figures/singleFigures/figure5c-e.pdf", fig5, width=18, height = 9)
