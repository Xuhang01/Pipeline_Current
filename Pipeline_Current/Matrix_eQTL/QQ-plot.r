#!/usr/bin/env Rscript
# gg_qqplot.R
# Kamil Slowikowski
# February 16, 2014

# Create a quantile-quantile plot with ggplot2.
#
# Assumptions:
#   - Expected P values are uniformly distributed.
#   - Confidence intervals assume independence between tests,
#     so they are conservative for GWAS with SNPs with LD.
#
# xs is a vector of p-values
# ci is the size of the confidence interval
library(ggplot2)
gg_qqplot = function(xs, ci=0.95) {
  N = length(xs)
  df = data.frame(observed=-log10(sort(xs)),
                  expected=-log10(1:N / N),
                  cupper=-log10(qbeta(ci,     1:N, N - 1:N + 1)),
                  clower=-log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
  log10Pe = expression(paste("Expected -log"[10], plain(P)))
  log10Po = expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), size=1.5, col="grey") +
    geom_abline(intercept=0, slope=1, alpha=0.4,size=0.8) +
    geom_line(aes(expected, cupper), col="red", alpha=0.4,size=0.8,linetype=2) +
    geom_line(aes(expected, clower), col="green", alpha=0.4,size=0.8,linetype=2) +
    xlab(log10Pe) +
    ylab(log10Po)+
    #theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
    theme_bw()+theme(panel.grid=element_blank())+
    ggtitle("QQ-plot for cis-eQTL p-values")+
    theme(plot.title=element_text(hjust=0.5))
}
work_dir <- commandArgs(T)[1]
output <- read.table(file.path(work_dir, "cis_eqtl.txt"),header=T)
colnames(output) <- c("SNP","gene","beta","t_stat","p_value","FDR")
pdf(file.path(work_dir,"QQ-plot.cis_eqtl.p_value.pdf"))
gg_qqplot(output$p_value)
dev.off()
