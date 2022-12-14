#v1.1 add in larger CZ sample from 2013
#v1.2 add in fruchterman reingold layout, scale nodes by trait Z score
#v1.3 clean up, focus on figs for Evolution (color integration index)
#     no node-based analysis!
require(igraph);require(PCIT);require(qgraph);require(psych);require(plotrix);require(assortnet);require(RColorBrewer)
require(devtools) #necessary to import files from dropbox over https

#Define Custom theme
#Make figure for Presentations
#Make custom theme for black background presentation in PPT
PPTtheme<-theme(panel.background = element_rect(fill = "transparent", colour = "transparent"),
    plot.background = element_rect(fill = "black",colour = "transparent"),
    panel.grband.minor = element_blank(), 
    panel.grband.major = element_blank(),#element_line(colour="white", size = 0.25),
    axis.line=element_line(size=.25,colour="white"),
    axis.text=element_text(size=14, colour="white"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title=element_text(size=18,face="bold",colour="white"),
        strip.background = element_rect(fill = "transparent", color = "white", size = 1),
    strip.text=element_text(colour="white",size=14))

## Make custom network plot function
Q<-function(COR,...){
  G<-qgraph(COR,diag=F,fade=F,label.color="black",label.font=2,label.scale=T,label.norm="0000",negCol="black",layout="spring",...)
return(G)}


#Read in my custom function definitions
source_url("https://dl.dropboxusercontent.com/s/b32xvn9x18gahc6/bioworkflow.R") #bioworkflow.R script
D<-read.csv("data/8pop.MF.TCS.csv")
head(D)
#check that there are no issues with data entry
check_class(D)
NA_outliers(D,id="band")
#Get rband of this indiv that has outlier values
D2<-NA_outliers(D,id="band")$newdata
NA_outliers(D2,id="band") #outliers switched to NA

#Remove duplicates
D2$popband<-paste(D2$Pop,D2$band)
sum(duplicated(D2$popband)) #4 duplicates

D2<-subset(D2,!duplicated(D2$popband))
sum(duplicated(D2$popband)) #Removed, now


#Generate trait names so you don't have to type em out
cat(names(D)[c(7,8,11:22)],sep='","')

#define traits of interest
toi<-c("mRWL","maxTS","T_Avg.Brightness","T_Hue","T_Chrom","R_Avg.Brightness","R_Hue","R_Chrom","B_Avg.Brightness","B_Hue","B_Chrom","V_Avg.Brightness","V_Hue","V_Chrom")
toi2<-c("RWL","RTS","TBri","THue","TChr","RBri","RHue","RChr","BBri","BHue","BChr","VBri","VHue","VChr")
shps<-c("square","square",rep("triangle",12))

# #Combine this more limited dataset w/ bigger samples from 4 pop study
# pop4<-read.csv("/Users/mattwilkins/Dropbox/My Research/My Papers/bars phenotypic variation/data/4popdata.csv")
# pop4_culled<-pop4[,c("band","Sex","Pop","Year","Mass",toi[-2])]
# #Don't have mean tail streamers for Turkey; FOR NOW, treat maxTS as meanTS to merge with 8pop dataset
# pop4_culled$mRTS<-pop4$maxTS
# 
# NA_outliers(pop4_culled,band="band")
# #get rband of couple indivbanduals w/ extreme values (NA them)
# pop4b<-NA_outliers(pop4_culled,band="band")$newdata
# NA_outliers(pop4b,band="band") #no outliers now
# 
# names(pop4b)
# names(D2)
# names(D2)[5]<-"Year"
# names(pop4b)[which(is.na(match(names(pop4b),names(D2))))]
# D3<-D2[,names(pop4b)]
D3<-D2
tapply(D3$band,D3$Pop,length)#sample sizes for each pop

#######################################
##########################
##Just use males
machos<-subset(D3,Sex=="M")

#Make list w/ each country as a separate dataframe
Dlist<-lapply(levels(machos$Pop),function(x) (subset(machos,Pop==x)))
(names(Dlist)<-levels(machos$Pop))

# Extract PCA objects for all pops, 5PCs, varimax rotation
# Although pops had 4-6 eigenvalues greater than 1; but dband 5, as compromise
pcaList<-lapply(1:8,function(x) try(principal(Dlist[[x]][,toi],rot="varimax",nfactors=5)))
names(pcaList)<-names(Dlist)

#Make list of correlations for each country (for traits of interest)
Clist<-lapply(levels(machos$Pop), function(x) cor(as.matrix(Dlist[[x]][,toi]),method="s",use="pairwise.complete"))
(names(Clist)<-levels(machos$Pop))

##### Filter ze networks!!! (PCIT)
# Get list of PCIT-filtered correlation matrices
Clist_pcit<-lapply(levels(machos$Pop),function(x) {
  old<-Clist[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$bandx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcit)<-levels(machos$Pop)

par(mfrow=c(2,4))
for (i in 1: length(Clist_pcit)){
  Q(Clist_pcit[[i]],title=names(Clist_pcit)[i])
}


# convert to 0-1 scale
Z_0to1<-cbind(machos[,c("band","Sex","Pop","Year")],apply(machos[,toi],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops
#This CENTERS, and SCALES means
Z<-cbind(machos[,c("band","Sex","Pop","Year")],apply(machos[,toi],2,function(x) scale(x,center=T,scale=T)))

#MEANS by Country
(Zmeans_0to1<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],mean)) #Z-scaled means for each trait in each pop
(Zmeans<-aggregate(.~Pop,data=Z[,c("Pop",toi)],mean)) #Z-scaled means for each trait in each pop

#calculate coefficients of variation (sd/mean)
CV<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],function(x) (sd(x,na.rm=T)/mean(x,na.rm=T)))
SD<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],function(x) (sd(x,na.rm=T)))


#Plot Z-trait boxplots (* Note Z_0to1 are not normal Z-scales, they're bounded 0:1)
require(reshape)
Zmelt<-melt(Z_0to1,band=c("band","Sex","Pop","Year"))
ggplot(Zmelt,aes(x=Pop,y=value))+geom_boxplot(col="white")+facet_wrap(~variable,nrow=2)+PPTtheme+xlab("\nPopulation")+ylab("Scaled Trait Value\n")+theme(axis.title.y=element_text(vjust=.5))
#ggsave("figs/Evo Presentation_Trait boxplots_BARS.png",wbandth=12)

#######
#-------------------------
# Plot 8 population networks for presentations w throat and breast colored differently
#pdf("figs/Evo presentation_8pop-pcit.pdf",wbandth=16,height=8)
#define vertex size
vertsize<-12
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
popnames<-c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: 8){
  mat<-Clist_pcit[[i]]
  Q(mat,color="white",border.color="gray20",labels=toi2,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="white",bg="transparent",label.color="white")
  mtext(paste0("n=",sum(!is.na(pcaList[[i]]$scores[,1]))),sbande=1,line=-.75,at=1.1,cex=.9,font=1,col="white")
  mtext(popnames[i],3,line=-.6,at=-.8,,adj=.5,col="white")
}
#dev.off()

# Same plot, w throat and breast colored differently
#pdf("figs/Evo presentation_8pop-pcit_R&T colored diff.pdf",wbandth=16,height=8)
#define vertex size
vertsize<-12
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
nodecols<-c("white","white",rep("orangered",3),rep("orange",9))
popnames<-c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: 8){
  mat<-Clist_pcit[[i]]
  Q(mat,color=nodecols,border.color="gray20",labels=toi2,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="white",bg="transparent",label.color="white")
  mtext(paste0("n=",sum(!is.na(pcaList[[i]]$scores[,1]))),sbande=1,line=-.75,at=1.1,cex=.9,font=1,col="white")
  mtext(popnames[i],3,line=-.6,at=-.8,,adj=.5,col="white")
}
#dev.off()

#######
###
# Calculate Body color Integration Index
bodytraits<-c("R_AvgBri","R_Hue","R_Chrom","B_AvgBri","B_Hue","B_Chrom","V_AvgBri","V_Hue","V_Chrom")
body<-lapply(Clist,function(x) {tmp<-x[bodytraits,bodytraits] #Get corr matrix for body traits; set diag to NA
  diag(tmp)<-NA
  return(tmp)})
nonbody<-lapply(Clist,function(x) {x[bodytraits,bodytraits]<-0 #Get corr matrix for body traits; set diag to NA
  diag(x)<-NA
  return(x)})

#Calculate |density| of whole network for each population
pop.netdensity<-sapply(Clist,function(x) sum(abs(x),na.rm=T)/2/91)
#Calculate |density| of this cluster for each
pop.clusterdensity<-sapply(body,function(x) (sum(abs(x),na.rm=T)/2)/36)
#Calculate similar index, but for noncluster
pop.nonclusterdensity<-sapply(nonbody,function(x) (sum(abs(x),na.rm=T)/2)/55)
#Calculate ratio of cluster to network density (integration index)
II<-pop.clusterdensity/pop.netdensity

### Calculate overall body brightness index
meanbri<-sapply(Dlist,function(x) {df<-x[,c("R_Chrom","B_Chrom","V_Chrom")]
            return(mean(colMeans(df,na.rm=T)))  })

Ibandf<-data.frame(II,pop.netdensity,pop.clusterdensity,pop.nonclusterdensity,meanbri,Pop=names(II))

jitterx<-meanbri+.009
jitterx[5]<-jitterx[5]+.0012
jitterx[3]<-jitterx[3]-.02
jittery<-pop.clusterdensity+.03
jittery[6]<-jittery[6]-.05
  
ggplot(data=Ibandf,aes(x=meanbri,y=pop.clusterdensity,label=Pop))+geom_point(size=2,aes(col=meanbri))+geom_text(size=5,col="white",x=jitterx,y=jittery)+xlab("Body Brightness Index")+ylab("Integration Index")+stat_ellipse(col="white")+PPTtheme+ylim(c(0,.7))+scale_colour_gradient(limits=c(.3, .50), low="#FFFFCC", high="#CC6600",gubande=F)+annotate("text",x=.32,y=.6,adj=0,label="rho= .809, p= 0.022",col="white",size=5) #paste0("rho == ~-.833~ p== 0.015"),parse=T,col="white")
cor.test(meanbri,pop.clusterdensity,method="s")
#ggsave("figs/Evolution Presentation_Color Integ Index~Body Chroma.jpg")

# ggplot(Ibandf,aes(x=meansumbri,y=pop.nonclusterdensity,label=Pop))+geom_point(size=2,aes(col=meansumbri))+geom_text(size=5,col="white",aes(x=jitter(meansumbri,10),y=pop.nonclusterdensity+.03))+xlab("Body Brightness Index")+ylab("Integration Index")+stat_ellipse(col="white")+PPTtheme+ylim(c(0,.7))+scale_colour_gradient(limits=c(70, 120), high="#FFFFCC", low="#CC6600",gubande=F)+annotate("text",x=130,y=.05,label="rho= 0.429, p= 0.299",col="white",size=5) 
# cor.test(meansumbri,pop.nonclusterdensity,method="s")
# ggsave("figs/Evolution Presentation_Color Integ Index~NonBody Brightness.jpg")


#############
##### Plot networks ordered by modularity
#reorder factor levels west to east
df.comb$Pop2<-factor(df.comb$Pop)
levels(df.comb$Pop2)=rank(pop.clusterdensity)
df.comb$Pop2<-factor(df.comb$Pop2,levels=8:1)



#pdf("figs/Evo presentation_8pop-pcit_ordered by clustering.pdf",wbandth=16,height=8)
#define vertex size
NewOrder<-match(1:8,rank(pop.clusterdensity))#order of pops by modularity

#vertsize<-12
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
nodecols<-c("white","white",rep("white",3),rep("orange",9))
popnames<-c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
#Change Brightness to Darkness so node scaling is in consistent direction
Zmeans.dk<-Zmeans_0to1
Zmeans.dk[,c(4,7,10,13)]<-apply(Zmeans.dk[,c(4,7,10,13)],2,function(x) {1-x})
#new labls, changing Bri to Dk, bc of reversed order
toi3<-c("RWL","RTS","TDk","THue","TChr","RDk","RHue","RChr","BDk","BHue","BChr","VDk","VHue","VChr")

for (i in 1: 8){
  mat<-Clist_pcit[[NewOrder[i] ]]
  scalar<-as.numeric(10*Zmeans.dk[NewOrder[i],-1])
  Q(mat,color=nodecols,border.color="gray20",labels=toi3,shape=shps,
    vsize=scalar+2+scalar^.8,
    edge.color="white",bg="transparent",label.color="white")
  
  mtext(popnames[NewOrder[i]],3,line=-.6,at=-.8,,adj=.5,col="white")
}
#dev.off()

require(car)
leveneTest(R_Chrom~Pop,machos)
boxplot(R_Chrom~Pop,machos)
tapply(machos$R_Chrom,machos$Pop,sd,na.rm=T)
