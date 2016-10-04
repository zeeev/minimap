#!/usr/bin/env Rscript
library("optparse")
require(ggplot2)

option_list = list(
            make_option(c("-f", "--file"),
            type="character",
            default=NULL,
            help="layout file",
            metavar="character"
            ),
            make_option(c("-p", "--pdf"),
            type="character",
            default=NULL,
            help="output file name (supports png or pdf)",
            metavar="character"
            ),
            make_option(c("-l", "--length"),
            type="numeric",
            default=10000,
            help="Filter alignments less than XXXBp [10000]",
            metavar="numeric"
            ),
            make_option(c("-d", "--dots"),
            type="numeric",
            default=20000,
            help="Filter dots less than XXXBp [20000]",
            metavar="numeric"
            ),
            make_option(c("-x", "--xlab"),
            type="character",
            default="target genome (Mb)",
            help="the xlab title",
            metavar="character"
            ),
            make_option(c("-y", "--ylab"),
            type="character",
            default="query genome (Mb)",
            help="the ylab title",
            metavar="character"
            ),
            make_option(c("-w", "--width"),
            type="numeric",
            default=10,
            help="The width of the plot",
            metavar="numeric"
            ),
            make_option(c("-i", "--height"),
            type="numeric",
            default=10,
            help="The height of the plot",
            metavar="numeric"
            ),
            make_option(c("-t", "--title"),
            type="character",
            default="layout",
            help="title of plot",
            metavar="character"
            ),
            make_option(c("-b", "--bed"),
            type="character",
            default=NULL,
            help="a bed file to paint on",
            metavar="character"
            )
);

options(error=traceback)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

dat<-read.table(opt$file)
dat<-dat[dat$V17 > opt$length, ]
tdat<-dat[dat$V17 > opt$dots,  ]

vdat<-aggregate(dat$V5, by=list(dat$V4), max)

thePlot<-ggplot()+geom_segment(data=dat,mapping=aes(y=V2/1e6, yend=V3/1e6, x=V5/1e6, xend=V6/1e6, colour=V7))
thePlot<-thePlot+theme_classic()
thePlot<-thePlot+theme(legend.position="none")

thePlot<-thePlot+labs(x=opt$xlab, y=opt$ylab, title=opt$title)
thePlot<-thePlot+geom_point(data=tdat, aes(x=V5/1e6, y=V2/1e6, colour=V7), size=0.5)
thePlot<-thePlot+geom_point(data=tdat, aes(x=V6/1e6, y=V3/1e6, colour=V7), size=0.5)

thePlot<-thePlot+geom_vline(xintercept=vdat$x/1e6, linetype="longdash", alpha=0.2)

if (!is.null(opt$bed)){
   bedDat<-read.table(opt$bed)
   bedDat<-bedDat[ bedDat$V1 %in% tdat$V4,]

   thePlot<-thePlot+geom_rect(data=bedDat, mapping=aes(ymin=0, ymax=max(dat$V3/1e6), xmin=V2/1e6, xmax=V3/1e6), alpha=0.25)
}


ggsave(filename=opt$pdf, width=opt$height, height=opt$height)

