require(ggplot2)

dat<-read.table("~/Desktop/test.txt")

tdat<-dat[dat$V16 > 10e3,]

vdat<-aggregate(dat$V5, by=list(dat$V4), max)

ggplot(dat, aes(y=V2/1e6, yend=V3/1e6, x=V5/1e6, xend=V6/1e6))+geom_segment()+labs(x="GRCh38", y="Clint")+geom_point(data=tdat, aes(x=V5/1e6, y=V2/1e6), size=0.5)+geom_point(data=tdat, aes(x=V6/1e6, y=V3/1e6), size=0.5)+theme_classic()+geom_vline(xintercept=vdat$x/1e6, linetype="longdash", alpha=0.2)

