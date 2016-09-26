#!/usr/bin/env Rscript

require(ggplot2)

args = commandArgs(trailingOnly=TRUE)

dat<-read.table(args[1])

tdat<-dat[dat$V16 > args[2],]

vdat<-aggregate(dat$V5, by=list(dat$V4), max)

thePlot<-ggplot(dat, aes(y=V2/1e6, yend=V3/1e6, x=V5/1e6, xend=V6/1e6))+geom_segment()
thePlot<-thePlot+labs(x="GRCh38", y="Clint")
thePlot<-thePlot+geom_point(data=tdat, aes(x=V5/1e6, y=V2/1e6), size=0.5)
thePlot<-thePlot+geom_point(data=tdat, aes(x=V6/1e6, y=V3/1e6), size=0.5)
thePlot<-thePlot+theme_classic()
thePlot<-thePlot+geom_vline(xintercept=vdat$x/1e6, linetype="longdash", alpha=0.2)

ggsave(filename="layout.pdf", width=10, height=10)

