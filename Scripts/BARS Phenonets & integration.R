require(igraph);require(PCIT);require(qgraph);require(psych);require(plotrix);require(assortnet);require(ggplot2); require(gtools); require(devtools);require(PHENIX);require(ggrepel)
#

####>>>>>>>>>>>>>>>>>>>>>
## Make custom plot function
Q<-function(COR,lab.col,lab.scale,lab.font,lay,...){
  if(missing(lab.col)){lab.col="black"}
  if(missing(lab.scale)){lab.scale=T}
  if(missing(lab.font)){lab.font=2}
  if(missing(lay)){lay="spring"}
  G<-qgraph(COR,diag=F,fade=F,label.color=lab.col,label.font=lab.font,label.scale=lab.scale,label.norm="0000",layout=lay,mar=c(4,7,7,4),...)
return(G)}
#<<<<<<<<<<<<<<<


#END FUNCTION DEFINITIONS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Read in my custom function definitions
source_url("https://dl.dropboxusercontent.com/s/b32xvn9x18gahc6/bioworkflow.R") #bioworkflow.R script (includes NA_outliers function)
D.raw<-read.csv("data/8pop.MF.TCS.csv")
#get rid of values for which sex is uncertain
D.raw$Sex[which(D.raw$Sex=="?"|D.raw$Sex=="M?"|D.raw$Sex=="F?"|D.raw$Sex=="F? "|D.raw$Sex=="N"|D.raw$Sex=="U"|D.raw$Sex=="TS"|D.raw$Sex=="")]<-NA


D.raw$Pop[which(D.raw$Pop=="Colorado")]<-"CO"
D.raw$Pop<-droplevels(D.raw$Pop)
ggplot(D.raw,aes(x=Pop,y=T_Avg.Brightness))+geom_boxplot()
#There are still a TON of Colorado outliers for several color axes
NA_outliers(subset(D.raw,Pop=="CO"),c(.01,.99),id="band")
#get rid of these bastards
COfixed<-NA_outliers(subset(D.raw,Pop=="CO"),id="band")$newdata
COfixed<-NA_outliers(subset(COfixed,Pop=="CO"),id="band")$newdata #run it again after getting rid of the extremely wrong vals
#rejoin with D.raw
toi<-c('band','Pop','Year','bandyear','Sex','CI1','RS','mRWL','mLTS','mRTS','maxTS','Mass','Date','T_Avg.Brightness','T_Hue','T_Chrom','R_Avg.Brightness','R_Hue','R_Chrom','B_Avg.Brightness','B_Hue','B_Chrom','V_Avg.Brightness','V_Hue','V_Chrom')
D<-rbind(subset(D.raw[,toi],Pop!="CO"),COfixed[,toi])
ggplot(D,aes(x=Pop,y=T_Avg.Brightness))+geom_boxplot()

pops<-levels(D$Pop)
NA_outliers(subset(D,Pop==pops[2]),c(.01,.99),id="band")#manually go thru all pops...look good enough, except Turkey


NA_outliers(subset(D,Pop=="TU"),c(.01,.9),id="band")
TUfixed<-NA_outliers(subset(D,Pop=="TU"),c(.01,.99),id="band")$newdata
#rejoin with D
D<-rbind(subset(D,Pop!="TU"),TUfixed)
NA_outliers(subset(D,Pop=="TU"),c(.01,.99),id="band")
#Leave this one alone

#Now do the whole dataset
NA_outliers(D,c(.01,.99),id="band",ignore="CI1")
#gonna leave these alone, since they aren't outliers in their respective populations

# #Get rid of outlier values
# D<-NA_outliers(D,c(.01,.99),id="band",ignore="CI1")$newdata
# NA_outliers(D,id="band",ignore="CI1") #outliers switched to NA

#OOOO K! A clean dataset.

### * Something weird happened with the color for those CO individuals...but now it's fixed

# toi<-c('band','Pop','Year','bandyear','Sex','CI1','mRWL','maxTS','Mass','T_tcs.h.theta','T_tcs.h.phi','T_tcs.r.achieved','T_sum.B2','R_tcs.h.theta','R_tcs.h.phi','R_tcs.r.achieved','R_sum.B2','B_tcs.h.theta','B_tcs.h.phi','B_tcs.r.achieved','B_sum.B2','V_tcs.h.theta','V_tcs.h.phi','V_tcs.r.achieved','V_sum.B2')

D2<-D
D2$Year<-as.factor(D2$Year)
head(D2)
#check that there are no issues with data entry
check_class(D2)
class(D2$RS)<-"numeric"

#************************************************************************
###### D2 contains repeated individuals; D3 will just have first capture
#Remove duplicates
D2$popID<-paste(D2$Pop,D2$band)
sum(duplicated(D2$popID)) #623 duplicates
#Order D2 by Population, then Year
D2<-D2[order(D2$Pop,D2$Year,D2$band),]
D2$PopYr<-factor(paste0(D2$Pop,sapply(D2$Year,function(x) substr(x,3,4))))

D3<-subset(D2,!duplicated(D2$popID))
sum(duplicated(D3$popID)) #No repeated individuals

###################################
#define traits of interest (toi2)
toi
(toi2<-toi[-c(1:13)])#list of color measures

shps<-c(rep("triangle",12))#"square","square",

#****** Get rid of rows w/ incomplete data*********
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
D3<-D3[which(complete.cases(D3[,c(toi2,"Sex")])),]
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
D3$Sex<-factor(D3$Sex)

(goodpops<-tapply(D3$band,D3$PopYr,length))#sample sizes for each pop (including both sexes)
(Ns<-ftable(Pop~Sex,data=D3))# sample sizes


# #remove pops w incomplete data
# D2<-subset(D2,PopYr=="CO08"|PopYr=="CO09"|PopYr=="CO10"|PopYr=="CO11"|PopYr=="CO12"|PopYr=="CO13")
# D2$PopYr<-factor(D2$PopYr)

#######################################
##########################
##Sep sexes
machos<-subset(D3,Sex=="M") #subset D3, which has no repeated individuals (only values for first capture)
hembras<-subset(D3,Sex=="F")
########################

#One last test to look for outliers by sex
testmachos<-lapply(levels(machos$Pop),function(x) NA_outliers(subset(machos,Pop==x),id="band")$changelog)
names(testmachos)<-levels(machos$Pop)
testmachos #these CO traits don't affect any color metrics

testhembras<-lapply(levels(hembras$Pop),function(x) NA_outliers(subset(hembras,Pop==x),id="band")$changelog)
names(testhembras)<-levels(hembras$Pop)
testhembras #Actually looks pretty good. Turkey has one outlier, but that could be real

####  Make list w/ each country as a separate dataframe
#Lumped across years
DlistM<-lapply(levels(machos$Pop),function(x) (subset(machos,Pop==x)))
(names(DlistM)<-levels(machos$Pop))
# Females
DlistF<-lapply(levels(hembras$Pop),function(x) (subset(hembras,Pop==x)))
(names(DlistF)<-levels(hembras$Pop))


#Separated by Years
# Males
DlistMyr<-lapply(levels(machos$PopYr),function(x) (subset(machos,PopYr==x)))
(names(DlistMyr)<-levels(machos$PopYr))
# Females
DlistFyr<-lapply(levels(hembras$PopYr),function(x) (subset(hembras,PopYr==x)))
(names(DlistFyr)<-levels(hembras$PopYr))
#NY02 & CZ10 don't have females, so take these out
DlistFyr[[which(names(DlistFyr)=="NY02")]]<-NULL
DlistFyr[[which(names(DlistFyr)=="CZ10")]]<-NULL

########################
#### Make list of correlations for each country (for traits of interest)
#----

##Lumped across years
# Males
ClistM<-lapply(names(DlistM), function(x) cor(as.matrix(DlistM[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistM)<-names(DlistM))
# Females
ClistF<-lapply(names(DlistF), function(x) cor(as.matrix(DlistF[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistF)<-names(DlistF))

#----
##Separated by Year
# Males
ClistMyr<-lapply(names(DlistMyr), function(x) cor(as.matrix(DlistMyr[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistMyr)<-names(DlistMyr))
# Females
ClistFyr<-lapply(names(DlistFyr), function(x) cor(as.matrix(DlistFyr[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistFyr)<-names(DlistFyr))


#igraph test
gCO10<-graph_from_adjacency_matrix(ClistMyr$CO10,"undirected",weighted=T,diag=F)
plot(gCO10,vertex.size=40,layout=layout.circle)

######## Plot networks
#######


#Unfiltered networks!
# png("figs/EvoFig_unfiltered.png")
#  qgraph(Clist$CO,labels=F,color=rep("gray40",14),shape=shps,vsize=8,vize2=8,edge.color="white",border.color="black",layout="spring",bg="transparent",diag=F)
#  dev.off()
 
# gCR<-Q(Clist$CR)
# gIA<-Q(Clist$IA)
# gIL<-Q(Clist$IL)
# gRM<-Q(Clist$RM)
# gTR<-Q(Clist$TR)
# gTW<-Q(Clist$TW)
# gUK<-Q(Clist$UK)

##### Filter ze networks!!! (PCIT)
# Get list of PCIT-filtered correlation matrices
#lumped across years
## Males
  Clist_pcitM<-lapply(names(ClistM),function(x) {
  old<-ClistM[[x]] #current dataset is xth element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcitM)<-names(ClistM)

## Females
 Clist_pcitF<-lapply(names(ClistF),function(x) {
  old<-ClistF[[x]] #current dataset is xth element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcitF)<-names(ClistF)

# #Plot M nets
# par(mfrow=c(2,4))
# for (i in 1: length(Clist_pcitM)){
#   Q(Clist_pcitM[[i]],title=names(Clist_pcitM)[i])
# }
# 
# #Plot F nets
# par(mfrow=c(2,4))
# for (i in 1: length(Clist_pcitF)){
#   Q(Clist_pcitF[[i]],title=names(Clist_pcitF)[i])
# }



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
################################
###### Convert trait values to 0-1 scale 
##
#MALES
rawmeansM<-aggregate(.~Pop,data=machos[,c("Pop",toi2)],mean)
#FEMALES
rawmeansF<-aggregate(.~Pop,data=hembras[,c("Pop",toi2)],mean)

#####
## convert to 0-1 scale
#MALES
Z_0to1M<-cbind(machos[,c("band","Sex","Pop","PopYr","Year")],apply(machos[,toi2],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops

#FEMALES
Z_0to1F<-cbind(hembras[,c("band","Sex","Pop","PopYr","Year")],apply(hembras[,toi2],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops

# #This CENTERS, and SCALES means
ZM<-cbind(machos[,c("band","Sex","Pop","PopYr","Year")],apply(machos[,toi2],2,function(x) scale(x,center=T,scale=T)))
ZF<-cbind(hembras[,c("band","Sex","Pop","PopYr","Year")],apply(hembras[,toi2],2,function(x) scale(x,center=T,scale=T)))

#MEANS by Country
ZMmeans<-aggregate(.~Pop,data=ZM[,c("Pop",toi2)],mean)
ZFmeans<-aggregate(.~Pop,data=ZF[,c("Pop",toi2)],mean)

#0-1 scaled means for each trait in each pop
#MALES
(Zmeans_0to1M<-aggregate(.~Pop,data=Z_0to1M[,c("Pop",toi2)],mean)) 
#FEMALEs
(Zmeans_0to1F<-aggregate(.~Pop,data=Z_0to1F[,c("Pop",toi2)],mean)) 

#Reverse brightness value...if you want
# Zmeans_0to1_revbri<-Zmeans_0to1
# Zmeans_0to1_revbri[,c(4,7,10,13)]<-apply(Zmeans_0to1_revbri[,c(4,7,10,13)],2,function(x) (1-x))

#calculate coefficients of variation (sd/mean) & SD
CVM<-aggregate(.~Pop,data=Z_0to1M[,c("Pop",toi2)],function(x) (sd(x,na.rm=T)/mean(x,na.rm=T)))
SDM<-aggregate(.~Pop,data=Z_0to1M[,c("Pop",toi2)],function(x) (sd(x,na.rm=T)))
CVF<-aggregate(.~Pop,data=Z_0to1F[,c("Pop",toi2)],function(x) (sd(x,na.rm=T)/mean(x,na.rm=T)))
SDF<-aggregate(.~Pop,data=Z_0to1F[,c("Pop",toi2)],function(x) (sd(x,na.rm=T)))

#####################################
#Plot Z-trait boxplots (* Note Z_0to1 are not normal Z-scales, they're bounded 0:1)
require(reshape2)
ZmeltM<-melt(Z_0to1M,id=c("band","Sex","Pop","Year","PopYr"))
ZmeltF<-melt(Z_0to1F,id=c("band","Sex","Pop","Year","PopYr"))

ggplot(ZmeltM,aes(x=Pop,y=value))+geom_boxplot()+facet_wrap(~variable,nrow=2)+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))
#ggsave("figs/Trait boxplots_BARS_males (PopYr).jpeg",width=12)
ggplot(ZmeltF,aes(x=PopYr,y=value))+geom_boxplot()+facet_wrap(~variable,nrow=2)+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))
#ggsave("figs/Trait boxplots_BARS_females (PopYr).jpeg",width=12)


#Test differences
for ( i in 1:length(toi2)){
  trait<-toi2[i]
  print(trait)
  print(teval(paste0("summary(aov(",trait,"~Pop,data=ZF))")))
  
} #All sig. different

#Traits of interest vectors
toi2
toi3<-c("TBri","THue","TChr","RBri","RHue","RChr","BBri","BHue","BChr","VBri","VHue","VChr")


#######
#-------------------------
# Plot Networks!!!
#MALES
pdf("figs/male networks.pdf",width=16,height=8)
#define vertex size
vertsize<-12#apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
popnames<-levels(D3$Pop)#sapply(levels(D3$PopYr),function(x) substr(x,1,2))#c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: length(Clist_pcitM)){
  mat<-Clist_pcitM[[i]]
  Q(mat,color="skyblue",border.color="gray20",labels=toi3,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="black",bg="white",label.color="black",label.cex=1.5)
  mtext(popnames[i],3,line=.5,at=-1.,,adj=.5,col="black")
 # mtext(paste0("n=",sum(!is.na(pcaListM[[i]]$scores[,1]))),side=3,line=-1.4,at=-1.,cex=.9,adj=.5,font=1,col="black")
}
dev.off()

#FEMALES
pdf("figs/female networks.pdf",width=16,height=8)
#define vertex size
vertsize<-12#apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
for (i in 1: length(Clist_pcitF)){
  mat<-Clist_pcitF[[i]]
  Q(mat,color="salmon",border.color="gray20",labels=toi3,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="black",bg="white",label.color="black",label.cex=1.5)
  mtext(popnames[i],3,line=.5,at=-1.,,adj=.5,col="black")
  #mtext(paste0("n=",sum(!is.na(pcaListF[[i]]$scores[,1]))),side=3,line=-1.4,at=-1.,cex=.9,adj=.5,font=1,col="black")

}
dev.off()


## NETWORK METRICS

## MALES
#Calculate |density| of whole network for each population
pop.netdensityM<-sapply(ClistM,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))

#calc degree
popwtdeg<-sapply(Clist_pcitM,function(x) sum(x[upper.tri(x)]!=0)/sum(upper.tri(x)))
popwtdeg

## FEMALES
#Calculate |density| of whole network for each population
pop.netdensityF<-sapply(ClistF,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))



###############################################################
###############################################################
### Calculate overall body brightness index (average chroma across 4 patches)

#MALES
darknessM<-sapply(DlistM,function(x) {df<-x[,c("T_Chrom","R_Chrom","B_Chrom","V_Chrom")]
            return(mean(colMeans(df,na.rm=T)))  })
integM<-data.frame(darknessM,pop.netdensityM,Pop=names(darknessM))
#FEMALES
darknessF<-sapply(DlistF,function(x) {df<-x[,c("T_Chrom","R_Chrom","B_Chrom","V_Chrom")]
            return(mean(colMeans(df,na.rm=T)))  })
integF<-data.frame(darknessF,pop.netdensityF,Pop=names(darknessF))

## Calculate Phenotypic Integration in the traditional way, based on var of eigenvalues of corr matrix correcting for # of traits & sample size
lapply(DlistM,function(x) pint(x[,toi2]))
(pintM<-unlist(lapply(DlistM,function(x) pint(x[,toi2])$PINT.c)))
(pintF<-unlist(lapply(DlistF,function(x) pint(x[,toi2])$PINT.c)))




tmp1<-integM
names(tmp1)<-c("darkness","netdensity","Pop")
tmp1$Sex="M"

tmp2<-integF
names(tmp2)<-c("darkness","netdensity","Pop")
tmp2$Sex="F"

integ<-rbind(tmp1,tmp2)
integ$pint<-c(pintM,pintF)

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# MAIN FIGURE

#MALES
(g1<-ggplot(data=integ,aes(x=darkness,y=pint,label=Pop))+geom_point(size=3,aes(col=darkness))+geom_point(size=3,shape=1,col="black")+geom_text_repel(size=5,col="black",aes(x=(darkness),y=(pint)),point.padding=.3,segment.color=NA)+xlab("Ventral Darkness")+ylab("Integration Index")+theme_bw()+scale_colour_gradient(limits=range(integ$darkness), low="#FFFFCC", high="#CC6600",guide=F)+facet_wrap(~Sex))#+annotate("text",x=12,y=.2,adj=0,label="rho= .809, p= 0.022",col="black",size=5) #paste0("rho == ~-.833~ p== 0.015"),parse=T,col="white")

ggsave("figs/Phenotypic integration by sex_pint metric.pdf",width=10,height =5 )

(g2<-ggplot(data=integ,aes(x=darkness,y=netdensity,label=Pop))+geom_point(size=3,aes(col=darkness))+geom_point(size=3,shape=1,col="black")+geom_text_repel(size=5,col="black",aes(x=(darkness),y=(netdensity)),point.padding=.3,segment.color=NA)+xlab("Ventral Darkness")+ylab("Network Density")+theme_bw()+scale_colour_gradient(limits=range(integ$darkness), low="#FFFFCC", high="#CC6600",guide=F)+facet_wrap(~Sex))#+annotate("text",x=12,y=.2,adj=0,label="rho= .809, p= 0.022",col="black",size=5) #paste0("rho == ~-.833~ p== 0.015"),parse=T,col="white")

ggsave("figs/Phenotypic integration by sex_network density.pdf",width=10,height =5 )

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

#Test the relationship
cor.test(darknessM,pop.netdensityM,method="k")
cor.test(darknessF,pop.netdensityF,method="k")
cor.test(darknessM,pintM,method="k")
cor.test(darknessF,pintF,method="k")

########
# Make funnel plots to test if either of these variables is affected by sample size
funnel<-as.data.frame(Ns)
funnel<-funnel[with(funnel,order(Sex)),]
funnel$darkness<-c(darknessF,darknessM)
funnel$netdensity<-c(pop.netdensityF,pop.netdensityM)
funnel

ggplot(melt(funnel,id=c("Sex","Pop","Freq")),aes(x=Freq,y=value,col=Pop,shape=Sex))+geom_point(size=2)+facet_grid(~variable) +xlab("Sample Size")+ylab("Value")+theme_bw()+theme(aspect.ratio = 3/4)

#No relationship with sample size

#############################################################################
###########

pintboth<-data.frame(pintM,pintF,darknessM,darknessF,Pop=names(pintM))#note this is out of order now...need to get some variables from code below
pintboth$Mts<-aggregate(maxTS~Pop,data=machos,FUN=mean)$maxTS
pintboth$Fts<-aggregate(maxTS~Pop,data=hembras,FUN=mean)$maxTS

## Test correlation between male & female traits
# Phen. Integ.
ggplot(pintboth,aes(x=pintM,y=pintF))+geom_point(size=3)+theme_bw()+geom_text(aes(x=pintM+0,y=pintF+.05,label=Pop))
cor.test(pintM,pintF) #M and F phenotypic integration for color not strongly correlated

ggplot(pintboth,aes(x=pop.netdensityM,y=pop.netdensityF))+geom_point(size=3)+theme_bw()+geom_text(aes(x=pop.netdensityM,y=pop.netdensityF+.005,label=Pop))
cor.test(pop.netdensityM,pop.netdensityF)


#Ventral Darkness
ggplot(pintboth,aes(x=darknessM,y=darknessF))+geom_point(size=3)+theme_bw()+geom_text(aes(x=darknessM+0,y=darknessF+.0025,label=Pop))
cor.test(darknessM,darknessF) #M & female *strongly* darkness correlated

#Tail streamer length
ggplot(pintboth,aes(x=Mts,y=Fts))+geom_point(size=3)+theme_bw()+geom_text(aes(x=Mts+0,y=Fts+1,label=Pop))+xlab("Male Tail Streamer Length (mm)")+ylab("Female Tail Streamer Length (mm)")
cor.test(pintboth$Mts,pintboth$Fts) #TS length almost perfectly correlated across sexes among pops

#get bootstrap CIs for each
# Time consuming
(pintMci<-lapply(DlistM,function(x) as.data.frame(pint.boot(x[,toi2],replicates=5000))) ) 
(pintFci<-lapply(DlistF,function(x) as.data.frame(pint.boot(x[,toi2],replicates=5000))) ) 
#test significance for each
(pintMp<-lapply(DlistM,function(x) pint.p(x[,toi2],n.replicates=5000)$Summary))
(pintFp<-lapply(DlistF,function(x) pint.p(x[,toi2],n.replicates=5000)$Summary))

#test calculation
names(princomp(cor(DlistM[[1]][,toi2])))
(wags<-unlist(lapply(DlistM,function(x) pint(x[,toi2])$PINT))) #Wagner pint measure (not controlling for anything..just variance of eigenvalues)
(eig<-unlist(lapply(DlistM,function(x) {Ex=eigen(cor(x[,toi2]))$values; return(sum((Ex-1)^2)/length(Ex))}))) #Note, pint calcs Ex (eigenvalue differences from 1 rather than mean)





#############
##### Plot networks ordered by modularity
#reorder factor levels west to east
machos$Pop2<-factor(machos$Pop)
levels(machos$Pop2)=rank(pop.netdensityM)


#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# OTHER MAIN FIGURES
# Plot networks ordered by integration
rawmeansF$Pop
countrynames<-c("US-Colorado","Czech Rep.","Israel","US-New York","Romania","Taiwan","Turkey","England")

pdf("figs/Male Networks_ordered.pdf",width=10.5,height=6)
#define vertex size
(NewOrderM<-order(rank(pop.netdensityM)))#order of pops by modularity

#vertsize<-12
par(mfrow=c(2,4),mar=rep(4,4),xpd=T,oma=rep(1,4),ps=18)
#specify palette
nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) #more realistic colors: #DED1BB","#7C472A"

#Calculate quantiles for each population's color values to color nodes
  scalar0<-sapply(names(rawmeansM)[-1],function(x) as.numeric(quantcut(rawmeansM[,x],q=50 ))) #make 50 quantiles for matching color scores
  rownames(scalar0)<-rawmeansM$Pop
  #reorder matrix
  scalar<-scalar0[NewOrderM,]
  scalar[,c(1:2,4:5,7:8,10:11)] <-51- scalar[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
#test that everything is as it should be
   cbind(normBri=scalar0[NewOrderM,1],revBri=scalar[,2],rawVal=rawmeansM[NewOrderM,2])

for (i in 1: length(Clist_pcitM)){
  mat<-Clist_pcitM[[NewOrderM[i] ]]
  nodecolor<-nodepal[scalar[i,]]
  Q(mat,color=nodecolor,border.color="gray20",labels=toi3,shape=shps,posCol="#181923",negCol=1,vsize=20,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,lay="spring")
  
  mtext(countrynames[NewOrderM[i]],3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)
}
dev.off()


#######
#fems
pdf("figs/Female Networks_ordered.pdf",width=10.5,height=6)
#define vertex size
(NewOrderF<-order(rank(pop.netdensityF)))#order of pops by modularity

#vertsize<-12
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
#Calculate quantiles for each population's color values to color nodes
  scalar0F<-sapply(names(rawmeansF)[-1],function(x) as.numeric(quantcut(rawmeansF[,x],q=50 ))) #make 50 quantiles for matching color scores
  rownames(scalar0F)<-rawmeansF$Pop
  #reorder matrix
  scalarF<-scalar0F[NewOrderF,]
  scalarF[,c(1:2,4:5,7:8,10:11)] <-51- scalarF[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
#test that everything is as it should be
   cbind(normBri=scalar0F[NewOrderF,1],revBri=scalarF[,2],rawVal=rawmeansF[NewOrderF,2])

for (i in 1: length(Clist_pcitF)){
  mat<-Clist_pcitF[[NewOrderF[i] ]]
  nodecolor<-nodepal[scalarF[i,]]
  Q(mat,color=nodecolor,border.color="gray20",labels=toi3,shape=shps,posCol="#181923",negCol=1,vsize=20,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,lay="spring")

  mtext(countrynames[NewOrderF[i]],3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)
}
dev.off()
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/



require(car)
leveneTest(R_Chrom~PopYr,machos) #variances aren't homogeneous
boxplot(R_Chrom~PopYr,machos)
tapply(machos$R_Chrom,machos$Pop,sd,na.rm=T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test for relationships between darkness, phenotypic integration and condition
# Predict that exaggeration comes from selection for condition-dependence, and this drives genetic integration
#*****!! Should probably control for capture date when calculating condition :/

use<-c("CO","CZ","IS","NY","TA","TU","UK")#names(goodpops)[which(goodpops>10)]

machos$condition=NA
dktraits<-c("T_Chrom","R_Chrom","B_Chrom","V_Chrom")
machos$darkness=rowMeans(machos[,dktraits])
  
integM$condep=NA# correlation between condition and darkness within populations
#calculate condition & condition dependence
for (i in 1: length(use))
{
  d<-subset(machos,Pop==use[i])
    mod<-lm(Mass~mRWL,data=d,na.action="na.exclude") #na.exclude pads NAs in residuals
  cond<-resid(mod)
  currows<-match(rownames(d),rownames(machos))
  machos$condition[currows]<-cond#assign condition scores for current set of individuals
  # machos$darkness[currows]<-d$darkness
    integM$condep[which(integM$Pop==use[i])]<-cor(d$darkness,cond,use="pairwise.complete.obs")
}


#add condition to integ df
#** Note, should probably aggregate condition by Pop, rather than year
#integM$condition<-aggregate(condition~Pop,data=machos,FUN=mean,na.rm=T)
integM$pint<-pintM
#now test condition ~ darkness
ggplot(machos,aes(x=darkness,y=condition,col=Pop))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Pop)

ggplot(machos,aes(x=maxTS,y=condition,col=Pop))+geom_point()+geom_smooth(method="lm")+facet_wrap(~Pop)

#Look at condition dependence and darkness/PI across pops
ggplot(integM,aes(x=darknessM,y=condep,col=darknessM,label=Pop))+geom_point(size=3,aes(col=darknessM))+geom_point(size=3,shape=1,col="black")+xlab("Ventral Darkness")+ylab("Condition Dependence")+theme_bw()+scale_colour_gradient(limits=range(integ$darkness), low="#FFFFCC", high="#CC6600",guide=F)+geom_text_repel(aes(text=Pop),col="black",point.padding = .2)

ggplot(integM,aes(x=pint,y=condep,col=darknessM,label=Pop))+geom_point(size=3,aes(col=darknessM))+geom_point(size=3,shape=1,col="black")+xlab("Phenotypic Integration")+ylab("Condition Dependence")+theme_bw()+scale_colour_gradient(limits=range(integ$darkness), low="#FFFFCC", high="#CC6600",guide=F)+geom_text_repel(aes(text=Pop),col="black",point.padding = .2)

#Nothing really there...



