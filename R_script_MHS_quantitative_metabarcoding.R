setwd("set your working directory")
library(ggplot2)
library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(ggplot2)
library(plyr)
library(estimatr)
library(vegan)
library(mgcv)
library(ggnewscale)
library(directlabels)
library(patchwork)
#===========Data Collation================
dat<-read.csv("zotu_table.csv")
dat1<-subset(dat,YN=="Y")
dat.std<-subset(dat1,name=="STD")
copy<-c(40,20,10,5,1)
sample<-colnames(dat.std)[4:ncol(dat.std)]
result<-c()
dat.res<-c()
for (i in sample) {
  mod1<-lm(dat.std[,i]~copy+0)
  mod2<-lm_robust(dat.std[,i]~copy+0)
  resultb<-data.frame(sample=i,
                      slope=summary(mod1)$coef[1,1],
                      sd=summary(mod1)$coef[1,2],
                      p=summary(mod1)$coef[1,4],
                      r.sq=summary(mod1)$r.squared,
                      adj.r.sq=summary(mod1)$adj.r.squared,
                      se.rob=summary(mod2)$coef[1,2],
                      p.rob=summary(mod2)$coef[1,4])
  result<-rbind(result,resultb)
  
  dat.resb<-data.frame(sample=i,x=ifelse(length(dat.std[,i]==5),
                                         c(40,20,10,5,1),c(40,20,10,5)),
                       res=mod1[["residuals"]],
                       fit=mod1[["fitted.values"]])
  dat.res<-rbind(dat.res,dat.resb)
}
hist(result$adj.r.sq,
     xlab = "Adjusted R-squared",
     family='serif',
     main=NULL)
write.csv(result,"MHS_qmeta_STD_result.csv",row.names = F)

#residual~fit
relm<-c()
for (i in unique(dat.res$sample)) {
  lm.resi<-lm(data=subset(dat.res,sample==i),res~fit+I(fit^2))
  relmb<-data.frame(sample=i,
                    inter=summary(lm.resi)$coef[1,1],
                    inte.p=summary(lm.resi)$coef[1,4],
                    slo=summary(lm.resi)$coef[2,1],
                    slo.p=summary(lm.resi)$coef[2,4],
                    slo2=summary(lm.resi)$coef[3,1],
                    slo.p2=summary(lm.resi)$coef[3,4])
  relm<-rbind(relm,relmb)
}
relm<-transform(relm,p=ifelse(slo.p2<=0.05|slo.p<=0.05,"Y","N"))

library(ggplot2)
ggplot(dat.res,aes(x=fit,y=res,alpha=0.1))+
  geom_point(aes(group=sample),size=0.5)+
  theme_bw()+xlab("Fitted values")+ylab("Residuals")+
  theme(legend.position ='none',
        panel.grid=element_blank(),
        text=element_text(size=20,family="serif")) +
  geom_line(data=subset(dat.res,sample%in%relm[relm$p=='N','sample']),
            aes(group=sample),stat='smooth',method='lm',
            formula=y~x,se=F,linewidth=0.5,alpha=0.2)+
  geom_line(data=subset(dat.res,sample%in%relm[relm$p=='Y','sample']),
            aes(group=sample),stat='smooth',method='lm',
            formula=y~x+I(x^2),se=F,linewidth=0.5,alpha=0.2)+
  geom_hline(yintercept = 0,linetype=2)+scale_x_log10()

hist(dat.res$res,breaks = 76,xlab="Residuals",
     main="Histogram of residuals",family='serif')
qqnorm(dat.res$res,main = "Normal Q-Q Plot of residuals",
       family="serif");qqline(dat.res$res)

#Converting reads to copy
dat.sam<-subset(dat1,name!="STD")
dat.copy<-data.frame()
for (i in sample) {
  dat.copyb<-data.frame(sample=i,
                        species=dat.sam$species,
                        copy=dat.sam[,i]/result[result$sample==i,"slope"])
  dat.copy<-rbind(dat.copy,dat.copyb)
}
dat.info<-read.csv("MHS-info.csv")
dat2<-merge(dat.info,dat.copy,by='sample')

#Subtract PCR negative control
dat.sn<-subset(dat2,type!='pcrnc')
dat.pnc<-subset(dat2,type=='pcrnc')
for (i in unique(dat.pnc$pcr)) {
  for (sp in unique(dat.pnc[dat.pnc$pcr==i,'species'])) {
    dat.sn[dat.sn$pcr==i&dat.sn$species==sp,'copy']<-
      dat.sn[dat.sn$pcr==i&dat.sn$species==sp,'copy']-
      dat.pnc[dat.pnc$pcr==i&dat.pnc$species==sp,'copy']
  }
}
dat.sn[dat.sn$copy<0,"copy"]<-0

#Calculate the ratio for filtering negative control
dat.fnc<-subset(dat.sn,type=='filternc')
dat.s<-subset(dat.sn,type=='sample')
fnc.ratio<-c()
for (i in unique(dat.fnc$date)) {
  fnc.ratio<-c(fnc.ratio,sum(dat.fnc[dat.fnc$date==i,'copy'])/
                 sum(dat.s[dat.s$date==i,'copy']))
}
mean(fnc.ratio)

#Calculate eDNA concentration (copy/L)
dat.s<-transform(dat.s,copy=copy*100*1000/filter)
dat.s<-dat.s[,colnames(dat.sn)!='filter'&colnames(dat.sn)!='pcr'&
               colnames(dat.sn)!='type']

#Copy < 1 is recorded as 0
dat.s[dat.s$copy<1,"copy"]<-0

#Fish species with detection count <2 were removed
library(plyr)
dat.species<-ddply(dat.s,.(species),summarise,number=sum(copy>0))
dat.species2<-dat.species[dat.species[,2]>1,1]
dat.s<-dat.s[dat.s$species %in% dat.species2,]
write.csv(dat.s,"MHS_qmeta_copy.csv",row.names = F)

#zotu table base on copy
rm(list = ls())
dat<-read.csv("MHS_qmeta_copy.csv")
dat2<-data.frame(sample=unique(dat$sample))
for (i in unique(dat$species)) {
  dat2b<-subset(dat,species==i)
  dat2b[,i]<-dat2b[,'copy']
  dat2b<-dat2b[,c('sample',i)]
  dat2<-merge(dat2,dat2b,by='sample')
}
write.csv(dat2,"zotu_table_copy.csv",row.names = F)

#======phylogenetic tree based on the neighbor-joining method===
library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(ggplot2)
test.dna <- read.dna('fish_12s_seq.fasta', format = 'fasta')
test.phyDat <- phyDat(test.dna, type = 'DNA', levels = NULL)
dna.dist <- dist.ml(test.phyDat,model = 'JC69')
tree.nj <- NJ(dna.dist)
p.nj<-ggtree(tree.nj,size=0.1)+
  geom_tiplab(align = T,fontface=3,size=2.9,
              linesize=0.1,family='serif')+
  xlim(NA,0.73)
p.nj

library(plyr)
dat.e<-read.csv("zotu_table_copy.csv",row.names = 1)
dat.ei<-read.csv("MHS-info.csv")
dat.ei$year<-substr(dat.ei$date,1,4)
dat.c<-read.csv("MHS_capture.csv",row.names = 1)

dat.e2<-data.frame()
for (i in c(2019,2020)) {
  dat.eb<-dat.e[row.names(dat.e) %in% dat.ei[dat.ei$year==i,'sample'],]
  dat.ec<-data.frame(Species=colnames(dat.eb),
                     Ratio=colSums(dat.eb)/sum(dat.eb),
                     Type=paste0("eDNA",i))
  dat.e2<-rbind(dat.e2,dat.ec)
}

dat.c2<-data.frame()
for (i in row.names(dat.c)) {
  dat.cb<-data.frame(Species=colnames(dat.c),
                     Ratio=colSums(dat.c[i,])/sum(dat.c[i,]),
                     Type=i)
  dat.c2<-rbind(dat.c2,dat.cb)
}
dat.m<-rbind(dat.c2,dat.e2)
for (i in c(1:nrow(dat.m))) {
  sp<-dat.m[i,"Species"]
  if (length(grep("_sp.",sp))==1) {
    dat.m[i,'Species']<-gsub("_"," ",sp)
  }else{
    dat.m[i,'Species']<-paste0(substr(sp,1,1),". ",
                               strsplit(sp,"_")[[1]][2])}
}
dat.m[dat.m$Species=="R. lagowskii","Species"]<-"R. steindachneri"
ttla<-data.frame(Type=c("y2007","y2008","y2009","y2010","y2019",
                        "eDNA2019","eDNA2020"),
                 New=c("T\n2007",
                       "T\n2008",
                       "T\n2009",
                       "T\n2010",
                       "T\n2019",
                       "e\n2019",
                       "e\n2020"))

for (i in unique(dat.m$Type)) {
  dat.m[dat.m$Type==i,"Type"]<-ttla[ttla$Type==i,"New"]
}

dat.m[dat.m$Ratio==0,'Ratio']<-NA
dat.m$Type<-factor(dat.m$Type,levels = c("T\n2007",
                                         "T\n2008",
                                         "T\n2009",
                                         "T\n2010",
                                         "T\n2019",
                                         "e\n2019",
                                         "e\n2020"))
treefish<-na.omit(data.frame(Fish=p.nj[["data"]][["label"]],
                             y=p.nj[["data"]][["y"]]))
treefish2<-treefish[order(-treefish$y),"Fish"]
treefish<-treefish[order(treefish$y),"Fish"]


dat.m$Species<-factor(dat.m$Species,levels = treefish)

p.ce<-ggplot(data=dat.m,aes(x=Type,y=Species,size=Ratio))+
  geom_point(na.rm=T)+
  scale_size_continuous(breaks = c(0.05,0.15,0.25,0.35),
                        labels = c("5%","15%","25%","35%"))+
  theme(axis.text.y = element_blank(),
        text=element_text(size=10,family="serif"))+
  xlab(NULL)+ylab(NULL)+scale_fill_grey(start = 1,end = 0)
p.ce
p.tree<-p.nj+p.ce
p.tree


#eDNA ~ traditional surveys=====2019 only
rm(list = ls()[ls()!='treefish2'])
dat.e<-read.csv("zotu_table_copy.csv",row.names = 1)
dat.ei<-read.csv("MHS-info.csv")
dat.ei$year<-substr(dat.ei$date,1,4)
dat.c<-read.csv("MHS_capture.csv",row.names = 1)
dat.c['y2019','Carassius_cuvieri']<-0
dat.c['y2019','Misgurnus_anguillicaudatus']<-0

dat.e2019<-data.frame(Species=colnames(dat.e),
                      copy=colSums(
                        dat.e[row.names(dat.e) %in% 
                                dat.ei[dat.ei$year==2019,'sample'],]))
dat.t2019<-data.frame(Species=colnames(dat.c),
                      indi=colSums(dat.c["y2019",]))
dat.ci2019<-merge(dat.e2019,dat.t2019,by="Species",all=T)
dat.ci2019[is.na(dat.ci2019)]<-0
dat.ci2019<-dat.ci2019[rowSums(dat.ci2019[,c(2,3)])!=0,]
library(estimatr)
mod.lm3<-lm_robust(data=dat.ci2019,(copy)~(indi))
summary(mod.lm3)
ggplot(data=dat.ci2019,aes(x=(indi),y=(copy)))+
  stat_smooth(method="lm_robust",formula=y~x,col="black")+
  geom_point()+
  theme(text=element_text(size=15,family="serif"))+
  ylab("eDNA survey (copies/L)")+
  xlab("Traditional survey (individuals)")

#eDNA ~ traditional surveys=====multiple years
dat.e.all<-data.frame(Species=colnames(dat.e),copy=colSums(dat.e))
dat.t.all<-data.frame(Species=colnames(dat.c),indi=colSums(dat.c))
dat.et.all<-merge(dat.e.all,dat.t.all,by="Species",all=T)
dat.et.all[is.na(dat.et.all)]<-0
mod.lm4<-lm_robust(data=dat.et.all,(copy)~(indi))
summary(mod.lm4)
ggplot(data=dat.et.all,aes(x=(indi),y=(copy)))+
  stat_smooth(method="lm_robust",formula=y~x,col="black")+
  geom_point()+
  theme(text=element_text(size=15,family="serif"))+
  ylab("eDNA survey (copies/L)")+
  xlab("Traditional survey (individuals)")

#=============Beta diversity===============
rm(list = ls()[ls()!='treefish2'])
dat<-read.csv("zotu_table_copy.csv",row.names = 1)
dat<-log(dat+1)
dat.s<-read.csv("MHS_qmeta_copy.csv")
dat.s<-dat.s[,c(1:6)]
dat.s<-unique(dat.s)
Sys.setlocale("LC_TIME", "English")
dat.s[,'mon']<-months(as.Date(dat.s$date))
dat.s$mon<-factor(dat.s$mon,
                  levels = c("March","April","May",
                             "June","July","August"))
dat.s[,'Year']<-substr(dat.s$date,1,4)
row.names(dat.s)<-dat.s$sample

#pcoa
library(vegan)
library(ggplot2)
bray_dist<-vegdist(dat,method = "bray")
pcoa1<-cmdscale(bray_dist,k=2,eig=T)
ggplot(data.frame(eig=pcoa1$eig),aes(x=1:141,y=eig))+
  geom_point()+geom_line()+ylab("Eigenvalues")+xlab("")+
  scale_x_continuous(breaks = c(1,50,100),
                     labels = paste("PCoA",c(1,50,100)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.1, vjust = 0.5),
        text=element_text(size=15,family="serif"))

df.pcoa<-data.frame(pcoa1$point)[1:2]
df.pcoa$sample<-row.names(df.pcoa)
df.plot<-merge(df.pcoa,dat.s,by='sample',all.x=T)
df.plot.sp<-data.frame(wascores(pcoa1$points[,1:2], dat,expand = T))
for (i in c(1:nrow(df.plot.sp))) {
  sp<-row.names(df.plot.sp)[i]
  if (length(grep("sp.",sp))==1) {
    df.plot.sp[i,'species']<-gsub("_"," ",sp)
  }else{
    df.plot.sp[i,'species']<-paste0(substr(sp,1,1),". ",
                                    strsplit(sp,"_")[[1]][2])}
}

pcoa_exp <- summary(eigenvals(pcoa1))['Proportion Explained',1:2]
axis1 <- 'PCoA axis1'
axis2 <- 'PCoA axis2'

p1<-ggplot(data = df.plot,aes(X1, X2)) +
  geom_point(aes(shape = site,color = mon),size=3)+
  stat_ellipse(aes(fill=mon,color = mon),
               geom="polygon",level=0.8,alpha=0.1)+
  labs(color="Month",fill="Month",shape="Site",
       x = axis1, y = axis2, 
       title = 'a. PCoA based on Bray-Curtis dissimilarity')+
  theme_bw()+
  theme(text=element_text(size=20,family="serif"),
        panel.grid=element_blank())+
  geom_segment(data = df.plot.sp, aes(x = 0, xend = X1, y=0, yend = X2), 
               size = 0.3, 
               arrow = arrow(length = unit(0.01, "npc")), alpha = 0.5) +
  geom_text(data = df.plot.sp,
            aes(x = X1+c(0,0,0.11,0,0.14,0,0.11,0,0,0.07,0.095,0.01,0,0,0),
                y = X2+c(0,0,0,-0.008,0,0,0,0,0,-0.01,0,0,0,0,0), 
                label = species,
                hjust = "inward", vjust = "outward"),
            color = "black", size = 6,family="serif",fontface = "italic")
p1

#permanover
veg.bray <- vegdist(dat,method="bray")
per1<-adonis2(veg.bray~mon+Year+site,data=dat.s,
              permutations = 9999,by='margin')
per1

group_name <- levels(dat.s$mon)

dat.paried <- data.frame()
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(dat.s, mon %in% 
                         c(as.character(group_name[i]), 
                           as.character(group_name[j])))
    otu_ij <- dat[group_ij$sample, ]
    adonis_result_otu_ij <- adonis2(otu_ij~mon, group_ij, 
                                    permutations = 9999, 
                                    distance = 'bray')
    res.temp <- as.data.frame(adonis_result_otu_ij)[1,]
    rownames(res.temp) <- paste(as.character(group_name[i]),'-',
                                as.character(group_name[j]))
    
    dat.paried <- rbind(dat.paried,res.temp)
  }
}
write.csv(dat.paried,file = 'sup_permanover_paired.csv')

#betadisper
disper<-betadisper(veg.bray,dat.s$mon)
permutest(disper,permutations = 9999)
dat.dis<-data.frame(dis=disper$distances,Month=dat.s$mon)
ggplot(dat.dis,aes(x=Month,y=dis))+geom_boxplot()+
  ylab("Distance to centroid")+
  geom_text(data=aggregate(dat.dis,dis~Month,max),
            aes(x=Month,y=dis+0.05),size=10,family="serif",
            label=c('a','a','b','b','b','b'))+
  theme_bw()+
  theme(text=element_text(size=20,family="serif"),
        panel.grid=element_blank())

dat.tukey<-round(TukeyHSD(disper)[['group']],4)
write.csv(dat.tukey,file = 'sup_betadisper_paired.csv')

disper<-betadisper(veg.bray,dat.s$Year)
permutest(disper,permutations = 9999)

disper<-betadisper(veg.bray,dat.s$site)
permutest(disper,permutations = 9999)

#ordisurf
library(ggnewscale)
library(directlabels)
#Water temperature
nm.ordi<-ordisurf(pcoa1~dat.s$Tem,npoints=100,method='REML')
nm.grid<-nm.ordi$grid
nm.ordi.coord<-expand.grid(x=nm.grid$x,y=nm.grid$y)
nm.ordi.coord$z<-as.vector(nm.grid$z)
nm.ordi.coord<-data.frame(na.omit(nm.ordi.coord))
lab.nm<-paste0("DE = ",round(summary(nm.ordi)$dev.expl,3)*100,"%")
p2<-ggplot(data = df.plot,aes(X1, X2)) +
  geom_point(aes(shape = site,color = mon),size=3)+
  labs(color="Month",fill="Month",shape="Site",
       x = 'PCoA axis1', y = 'PCoA axis2',title = 
         expression(paste('b. ',"Water temperature (°C), ", 
                          italic('p'), "< 0.001, DE = 51.7%")))+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(text=element_text(size=20,family="serif"),
        legend.position = 'none')+
  new_scale_color()+
  stat_contour(data=nm.ordi.coord,aes(x,y,z=z,
                                      colour=..level..))+
  scale_colour_gradient(low = "green",high = "red",name="Tem")
p2<-direct.label(p2,dl.combine(list("first.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black"),
                               list("last.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black")))
#pH
nm.ordi<-ordisurf(pcoa1~dat.s$pH,npoints=100,method='REML')
nm.grid<-nm.ordi$grid
nm.ordi.coord<-expand.grid(x=nm.grid$x,y=nm.grid$y)
nm.ordi.coord$z<-as.vector(nm.grid$z)
nm.ordi.coord<-data.frame(na.omit(nm.ordi.coord))
lab.nm<-paste0("DE = ",round(summary(nm.ordi)$dev.expl,3)*100,"%")
p3<-ggplot(data = df.plot,aes(X1, X2)) +
  geom_point(aes(shape = site,color = mon),size=3)+
  labs(color="Month",fill="Month",shape="Site",
       x = 'PCoA axis1', y = 'PCoA axis2',
       title = expression(paste('c. ',"pH, ",
                                italic(p), "< 0.001, DE = 24.9%")))+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(text=element_text(size=20,family="serif"),
        legend.position = 'none')+
  new_scale_color()+
  stat_contour(data=nm.ordi.coord,aes(x,y,z=z,
                                      colour=..level..))+
  scale_colour_gradient(low = "green",high = "red",name="Tem")
p3<-direct.label(p3,dl.combine(list("first.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black"),
                               list("last.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black")))
#EC
nm.ordi<-ordisurf(pcoa1~dat.s$EC,npoints=100,method='REML')
nm.grid<-nm.ordi$grid
nm.ordi.coord<-expand.grid(x=nm.grid$x,y=nm.grid$y)
nm.ordi.coord$z<-as.vector(nm.grid$z)
nm.ordi.coord<-data.frame(na.omit(nm.ordi.coord))
lab.nm<-paste0("DE = ",round(summary(nm.ordi)$dev.expl,3)*100,"%")
p4<-ggplot(data = df.plot,aes(X1, X2)) +
  geom_point(aes(shape = site,color = mon),size=3)+
  labs(color="Month",fill="Month",shape="Site",
       x = 'PCoA axis1', y = 'PCoA axis2',
       title = expression(paste('d. ',"EC (mS/m), ",
                                italic(p)," = 0.002, DE = 14.5%")))+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(text=element_text(size=20,family="serif"),
        legend.position = 'none')+
  new_scale_color()+
  stat_contour(data=nm.ordi.coord,aes(x,y,z=z,
                                      colour=..level..))+
  scale_colour_gradient(low = "green",high = "red",name="Tem")
p4<-direct.label(p4,dl.combine(list("first.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black"),
                               list("last.points",cex=1.2,
                                    fontface="plain",fontfamily="serif",
                                    color="black")))

library(patchwork)
pm<-(p1+p2)/(p3+p4)+plot_layout(guides = "collect")
pm

#Hierarchical generalized additive model (HGAM)
rm(list = ls()[ls()!='treefish2'])
library(mgcv)
dat<-read.csv("MHS_qmeta_copy.csv")
dat<-cbind(dat[,c(1:2)],Year=substr(dat$date,1,4),dat[,c(3:8)])
colnames(dat)<-c("Sample","Date","Year","Site",
                 "Tem","pH","EC","Fish","Copy")
dat$Site<-as.factor(dat$Site)
dat$Fish<-as.factor(dat$Fish)
dat$Year<-as.factor(dat$Year)

copy.fish<-aggregate(dat,Copy~Fish,sum)
copy.fish<-copy.fish[order(-copy.fish$Copy),]
dat.8f<-subset(dat,Fish%in%copy.fish[1:8,'Fish'])
dat.8f$Date<-as.Date(dat.8f$Date)
dat.8f$Fish<-factor(dat.8f$Fish,levels=copy.fish[1:8,'Fish'])
Sys.setlocale("LC_TIME", "English")

palab<-c("Tem"="Water temperature (℃)",
         "pH"="pH",
         "EC"="EC (mS/m)")
filab<-c("Carassius_cuvieri"="C. cuvieri",
         "Carassius_sp."="Carassius sp.",
         "Ctenopharyngodon_idella"="C. idella",
         "Cyprinus_carpio"="C. carpio",
         "Gnathopogon_elongatus_elongatus"="G. elongatus",
         "Hemibarbus_barbus"="H. barbus",
         "Nipponocypris_temminckii"="N. temminckii",
         "Pseudorasbora_parva"="P. parva",
         "Zacco_platypus"="Z. platypus",
         "Misgurnus_anguillicaudatus"="M. anguillicaudatus",
         "Misgurnus_sp._Clade_A"="Misgurnus sp. Clade A",
         "Tribolodon_hakonensis"="T. hakonensis",
         "Lefua_echigonia"="L. echigonia",
         "Gymnogobius_urotaenia"="G. urotaenia",
         "Rhinogobius_sp."="Rhinogobius sp.",
         "Micropterus_salmoides"="M. salmoides",
         "Lepomis_macrochirus"="L. macrochirus",
         "Silurus_asotus"="S. asotus",
         "Anguilla_japonica"="A. japonica")

ggplot(subset(dat.8f,Copy>0),aes(x=Date,y=Copy,col=Site))+
  geom_point(size=3)+
  facet_grid(Fish~Year,scales = "free",
             labeller = labeller(Fish=filab))+
  scale_y_log10(labels = scales::comma)+
  scale_x_date(breaks = "1 month",date_labels = "%B")+
  theme_bw()+
  xlab("Date")+ylab("eDNA concentration (Copy/L)")+
  theme(text=element_text(size=20,family="serif"),
        strip.text.y = element_text(face = "italic"))

modGS<-gam(data=subset(dat,Copy>=1),family='gaussian',
           log(Copy)~
             s(Tem,m=2)+s(Tem,Fish,bs='fs',m=2)+
             s(pH,m=2)+s(pH,Fish,bs='fs',m=2)+
             s(EC,m=2)+s(EC,Fish,bs='fs',m=2)+
             s(Site,Year,bs='re'),method='REML')
save(modGS,file="MHS_modGS.rda",compress = 'xz')
load("MHS_modGS.rda")
summary(modGS)
gam.check(modGS)

te<-data.frame(res=modGS$residuals,
               lip=modGS$linear.predictors,
               response=modGS$y)
ggplot(te,aes(x=lip,y=res))+geom_point(alpha=0.2)+
  geom_smooth(method='lm',formula=y~x)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_bw()+
  xlab('Linear predictors')+ylab('Residuals')+labs(title = 'a')+
  theme(panel.grid=element_blank(),
        text=element_text(size=20,family="serif"))

ggplot(te,aes(x=res))+
  geom_histogram(fill='grey',stat='bin',binwidth = 0.25,
                 col='black',alpha=0.5)+
  theme_bw()+xlab('Residuals')+ylab('Frequency')+labs(title = 'b')+
  theme(panel.grid=element_blank(),
        text=element_text(size=20,family="serif"))

ggplot(te,aes(sample=res))+geom_qq(shape=1,alpha=0.2)+geom_qq_line()+
  theme_bw()+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  labs(title = 'c')+
  theme(panel.grid=element_blank(),
        text=element_text(size=20,family="serif"))

ggplot(te,aes(x=lip,y=response))+geom_point(alpha=0.2)+
  geom_abline(slope=1,linetype=2)+
  theme_bw()+xlab('Fitted Values')+ylab('Response')+labs(title = 'd')+
  theme(panel.grid=element_blank(),
        text=element_text(size=20,family="serif"))

#partial residual of HGAM
lev1.5<-levels(modGS[["model"]][["Fish"]])
lev2<-c()
for (i in lev1.5) {
  lev2<-c(lev2,rep(i,100))
}

dat.plot<-plot(modGS,select=1,residuals = T,se=1.96)
dat.p.raw<-data.frame()
dat.p.fit<-data.frame()
for (i in c(1,3,5)) {
  pa<-c("Tem","X","pH","X","EC")[i]
  dat.p.rawb<-data.frame(rawx=dat.plot[[i]]$raw,
                         PA=pa,
                         Fish=na.omit(dat)[na.omit(dat)$Copy>=1,"Fish"],
                         resi1=dat.plot[[i]]$p.resid,
                         resi2=dat.plot[[i+1]]$p.resid)
  dat.p.raw<-rbind(dat.p.raw,dat.p.rawb)
  dat.p.fitb<-data.frame(fit1=dat.plot[[i]]$fit,
                         fit2=dat.plot[[i+1]]$fit,
                         fitx=dat.plot[[i]]$x,
                         Fish=lev2,
                         se=dat.plot[[i]]$se,
                         PA=pa)
  dat.p.fit<-rbind(dat.p.fit,dat.p.fitb)
}
dat.p.raw$Fish<-as.character(dat.p.raw$Fish)
for (i in c(1:nrow(dat.p.raw))) {
  sp<-dat.p.raw[i,"Fish"]
  if (length(grep("_sp.",sp))==1) {
    dat.p.raw[i,'Fish']<-gsub("_"," ",sp)
  }else{
    dat.p.raw[i,'Fish']<-paste0(substr(sp,1,1),". ",
                                strsplit(sp,"_")[[1]][2])}
}
dat.p.fit$Fish<-as.character(dat.p.fit$Fish)
for (i in c(1:nrow(dat.p.fit))) {
  sp<-dat.p.fit[i,"Fish"]
  if (length(grep("_sp.",sp))==1) {
    dat.p.fit[i,'Fish']<-gsub("_"," ",sp)
  }else{
    dat.p.fit[i,'Fish']<-paste0(substr(sp,1,1),". ",
                                strsplit(sp,"_")[[1]][2])}
}


dat.p.raw$Fish<-factor(dat.p.raw$Fish,levels = treefish2)
dat.p.raw$PA<-factor(dat.p.raw$PA,levels = c("Tem","pH","EC"))
dat.p.fit$Fish<-factor(dat.p.fit$Fish,levels = treefish2)
dat.p.fit$PA<-factor(dat.p.fit$PA,levels = c("Tem","pH","EC"))

ggplot(data=dat.p.raw,aes(x=rawx,y=(resi1+resi2)))+
  facet_grid(Fish~PA,scales = "free",
             labeller = labeller(PA=palab))+
  geom_line(data=dat.p.fit,aes(x=fitx,y=(fit1+fit2)))+
  geom_ribbon(data=dat.p.fit,
              aes(ymin=(fit1+fit2 - se), 
                  ymax=(fit1+fit2 + se), 
                  x=fitx),alpha=0.3,inherit.aes=FALSE)+
  theme(text=element_text(size=10,family="serif"),
        strip.text.y = element_text(face = "italic"))+
  ylab("Component smooth function")+
  xlab("Environment variable")

ggplot(data=dat.p.raw,aes(x=rawx,y=(resi1+resi2)))+
  facet_grid(Fish~PA,scales = "free",
             labeller = labeller(PA=palab))+
  geom_point(size=0.1)+
  geom_line(data=dat.p.fit,aes(x=fitx,y=(fit1+fit2)))+
  geom_ribbon(data=dat.p.fit,
              aes(ymin=(fit1+fit2 - se), 
                  ymax=(fit1+fit2 + se), 
                  x=fitx),alpha=0.3,inherit.aes=FALSE)+
  theme(text=element_text(size=10,family="serif"),
        strip.text.y = element_text(face = "italic"))+
  ylab("Component smooth function")+
  xlab("Environment variable")

#================spawning activity================
rm(list = ls()[ls()!='treefish2'])
dat<-read.csv('MHS_qmeta_copy.csv')
dat<-cbind(dat[,c(1:2)],Year=substr(dat$date,1,4),dat[,c(3:8)])
colnames(dat)<-c("Sample","Date","Year","Site",
                 "Tem","pH","EC","Fish","Copy")
dat[(dat$Copy)==0,"Copy"]<-NA
dat<-transform(dat,Con.c=log(Copy))

datb<-data.frame()
fishname<-unique(dat$Fish)
sitename<-unique(dat$Site)
yearname<-unique(dat$Year)
for (fish in fishname) {
  for (site in sitename) {
    for (year in yearname) {
      dat1<-subset(dat,Fish==fish&Site==site&Year==year)
      q3=quantile(dat1$Copy,na.rm = T)[4]
      q1=quantile(dat1$Copy,na.rm = T)[2]
      out=q3+1.5*(q3-q1)
      dat1$spawning<-ifelse(dat1$Copy>out,ifelse(dat1$Tem<30,1,0),0)
      datb<-rbind(datb,dat1)
    }
  }
}

datb$Fish<-as.factor(datb$Fish)
library(plyr)
Fishspawning<-ddply(datb,.(Fish),summarise,spyn=sum(na.omit(spawning)))
Fishspawning<-Fishspawning[Fishspawning$spyn>0,'Fish']
datb<-datb[datb$Fish %in% Fishspawning,]

mod.sp2<-gam(data=datb,family = 'binomial',method = 'REML',
             spawning~Con.c+Con.c*Fish+
               s(Tem,m=2)+s(Tem,Fish,bs='fs',m=2))
summary(mod.sp2)

ree<-summary(mod.sp2)
ree.est<-data.frame(est=ree$p.coeff)
ree.se<-data.frame(se=ree$se)
ree.z<-data.frame(z.value=ree$p.t)
ree.p<-data.frame(p=ree$p.pv)

ree.m<-cbind(ree.est,se=ree.se[1:26,],ree.z,ree.p)
write.csv(ree.m,file='sup_result of modsp2.csv')

pred.re.sp<-cbind(dat[dat$Fish %in% Fishspawning,],
                  pro=predict(mod.sp2,type='response',
                              dat[dat$Fish %in% Fishspawning,],
                              se.fit = F))
pred.re.sp$YearSite<-paste0(pred.re.sp$Year,"-",pred.re.sp$Site)
pred.re.sp$DateNew<-pred.re.sp$Date
pred.re.sp[pred.re.sp$Date=="2020/4/27","DateNew"]<-"2020/4/28"
pred.re.sp[pred.re.sp$Date=="2020/5/6","DateNew"]<-"2020/5/5"
pred.re.sp[pred.re.sp$Date=="2020/8/17","DateNew"]<-"2020/8/18"
pred.re.sp[pred.re.sp$Year=="2020",
           "DateNew"]<-gsub("2020","2019",
                            pred.re.sp[pred.re.sp$Year=="2020",
                                       "DateNew"])
pred.re.sp$DateNew<-as.Date(pred.re.sp$DateNew)
pred.re.sp[pred.re.sp$Year=="2020",
           "DateNew"]<-pred.re.sp[pred.re.sp$Year=="2020",
                                  "DateNew"]+2
for (i in c(1:nrow(pred.re.sp))) {
  sp<-pred.re.sp[i,"Fish"]
  if (length(grep("_sp.",sp))==1) {
    pred.re.sp[i,'Fish']<-gsub("_"," ",sp)
  }else{
    pred.re.sp[i,'Fish']<-paste0(substr(sp,1,1),". ",
                                 strsplit(sp,"_")[[1]][2])}
}

pred.re.sp[is.na(pred.re.sp$pro),"pro"]<-0
pred.re.sp$Fish<-factor(pred.re.sp$Fish,levels = treefish2)
ggplot(data=pred.re.sp,aes(x=DateNew,y=YearSite))+
  facet_grid(Fish~.)+theme_bw()+
  geom_tile(aes(fill=pro))+theme(panel.grid = element_blank())+
  scale_fill_gradient(name="Probability",
                      low = "white", 
                      high = "black")+
  scale_x_date(date_breaks = '1 month',date_labels = "%m.%d")+
  theme(text=element_text(size=10,family="serif"),
        strip.text.y = element_text(face = "italic"))+
  xlab("Date")+ylab("Year & Site")
