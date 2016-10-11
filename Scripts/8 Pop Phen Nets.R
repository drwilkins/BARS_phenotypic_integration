#git config --global user.email "mrwilkins06@gmail.com"

require(igraph);require(PCIT);require(qgraph);require(psych);require(plotrix);require(assortnet)
require(devtools) #necessary to import files from dropbox over https
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }} #end multiplot definition

#Read in my custom function definitions
source_url("https://dl.dropboxusercontent.com/s/b32xvn9x18gahc6/bioworkflow.R") #bioworkflow.R script
D.raw<-read.csv("data/8pop.MF.TCS.csv")
quotenames(D.raw)
toi<-c('band','Pop','Year','bandyear','Sex','CI1','mRWL','maxTS','Mass','T_tcs.h.theta','T_tcs.h.phi','T_tcs.r.achieved','T_sum.B2','R_tcs.h.theta','R_tcs.h.phi','R_tcs.r.achieved','R_sum.B2','B_tcs.h.theta','B_tcs.h.phi','B_tcs.r.achieved','B_sum.B2','V_tcs.h.theta','V_tcs.h.phi','V_tcs.r.achieved','V_sum.B2')
D<-D.raw[,toi]
head(D)
#check that there are no issues with data entry
#class(D$Year)<-"factor"
check_class(D)
NA_outliers(D,c(.01,.99),id="band",ignore="CI1")
#Get rid of outlier values
D2<-NA_outliers(D,c(.01,.99),id="band",ignore="CI1")$newdata
NA_outliers(D2,id="band") #outliers switched to NA


#************************************************************************
###### D2 contains repeated individuals; D3 will just have first capture
#Remove duplicates
D2$popID<-paste(D2$Pop,D2$band)
sum(duplicated(D2$popID)) #620 duplicates
#Order D2 by Population, then Year
D2<-D2[order(D2$Pop,D2$Year),]
D3<-subset(D2,!duplicated(D2$popID))
sum(duplicated(D3$popID)) #Removed, now

###################################
#define traits of interest (toi2)
toi
(toi2<-toi[-c(1:6,9,11,15,19,23)])#get rid of descriptors & phi measures
shps<-c("square","square",rep("triangle",12))


tapply(D3$band,D3$Pop,length)#sample sizes for each pop

#reorder factor levels west to east
#D3$Pop<-factor(df.comb$Pop,levels=c("CO","IA","UK","CR","RM","TR","IL","TW"))

#######################################
##########################
##Sep sexes

machos<-subset(D3,Sex=="M")
hembras<-subset(D3,Sex=="F")
########################
####  Make list w/ each country as a separate dataframe
# Males
DlistM<-lapply(levels(machos$Pop),function(x) (subset(machos,Pop==x)))
(names(DlistM)<-levels(machos$Pop))
# Females
DlistF<-lapply(levels(hembras$Pop),function(x) (subset(hembras,Pop==x)))
(names(DlistF)<-levels(hembras$Pop))

########################
#### Make list of correlations for each country (for traits of interest)
# Males
ClistM<-lapply(levels(machos$Pop), function(x) cor(as.matrix(DlistM[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistM)<-levels(machos$Pop))
# Females
ClistF<-lapply(levels(hembras$Pop), function(x) cor(as.matrix(DlistF[[x]][,toi2]),method="s",use="pairwise.complete"))
(names(ClistF)<-levels(hembras$Pop))


### Make igraph understand triangle shapes
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
 
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                 plot=mytriangle)


#igraph test
gCO<-graph_from_adjacency_matrix(Clist$CO,"undirected",weighted=T,diag=F)
plot(gCO,vertex.size=40,layout=layout.circle,vertex.shape=shps)

######## Plot networks
#######
## Make custom plot function
Q<-function(COR,...){
  G<-qgraph(COR,diag=F,fade=F,label.color="black",label.font=2,label.scale=T,label.norm="0000",negCol="black",layout="spring",...)
return(G)}

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
## Males
  Clist_pcitM<-lapply(levels(machos$Pop),function(x) {
  old<-ClistM[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcitM)<-levels(machos$Pop)

## Females
 Clist_pcitF<-lapply(levels(hembras$Pop),function(x) {
  old<-ClistF[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcitF)<-levels(hembras$Pop)

#Plot M nets
par(mfrow=c(2,4))
for (i in 1: length(Clist_pcitM)){
  Q(Clist_pcitM[[i]],title=names(Clist_pcitM)[i])
}

#Plot F nets
par(mfrow=c(2,4))
for (i in 1: length(Clist_pcitF)){
  Q(Clist_pcitF[[i]],title=names(Clist_pcitF)[i])
}

### **** Some problems with the PCA for small populations (n approaches k vars to extract)
###Get PCA scores for all dataframes (Dlist)
# Count # factors w/ eigenvalues >=1
extractM<-vector(l=length(levels(D3$Pop)))
extractF<-vector(l=length(levels(D3$Pop)))
for( i in 1:length(levels(D3$Pop))){
pcaM<-principal(DlistM[[i]][,toi2],rot="none",nfactors=length(toi2))#Won't work for full set, calc 2 fewer eigenvectors
pcaF<-principal(DlistF[[i]][,toi2],rot="none",nfactors=length(toi2)-2)#WONT work, bc too few observations
print(names(DlistM)[i])
print(names(DlistF)[i])
extractM[i]<-sum(pcaM$values>=1)
extractF[i]<-sum(pcaF$values>=1)
print(pcaM)
print(pcaF)}
extractM
extractF#number of factors to extract for each population

# Extract PCA objects for all pops, 5PCs, varimax rotation
# Although pops had 4-6 eigenvalues greater than 1; but did 5, as compromise
pcaListM<-lapply(1:length(levels(D3$Pop)),function(x) principal(DlistM[[x]][,toi2],rot="varimax",nfactors=extractM[x]))
names(pcaListM)<-names(DlistM)
pcaListF<-lapply(1:length(levels(D3$Pop)),function(x) principal(DlistF[[x]][,toi2],rot="varimax",nfactors=5))
names(pcaListF)<-names(DlistF)

####
# Which correlations are in all networks? And what is the average value?
#MALES
Clist_pcitM_bin<-lapply(1:length(levels(D3$Pop)),function(x) apply(Clist_pcitM[[x]],1, function(xx) ifelse(xx==0,0,1)))#binary correlation matrices (as list)
names(Clist_pcitM_bin)<-names(Clist_pcitM)

#FEMALES
Clist_pcitF_bin<-lapply(1:length(levels(D3$Pop)),function(x) apply(Clist_pcitF[[x]],1, function(xx) ifelse(xx==0,0,1)))#binary correlation matrices (as list)
names(Clist_pcitF_bin)<-names(Clist_pcitF)

(commonedgesM<-teval(paste0("Clist_pcitM_bin[[",1:length(levels(D3$Pop)),"]]",collapse="+"))) #Adds all population binary edges together
(commonedgesF<-teval(paste0("Clist_pcitF_bin[[",1:length(levels(D3$Pop)),"]]",collapse="+"))) #Adds all population binary edges together

diag(commonedgesM)<-NA
colscaleM<-seq(1,length(levels(D3$Pop)),1) #color scale goes from 1 to 8; 8 if edge present in all pops (gray); 1 if only present once (black)
diag(commonedgesF)<-NA
colscaleF<-seq(1,length(levels(D3$Pop)),1) 

#####
## Which and how many edges don't exist?
#For MALES
idxM<-which(lower.tri(commonedgesM)&commonedgesM==0,arr.ind=T)
paste(rownames(commonedgesM)[idxM[,"row"]],colnames(commonedgesM)[idxM[,"col"]],sep="--") 
nrow(idxM)/length(which(upper.tri(Clist_pcitM[[1]])))
#16/91=17.6% not present in any pop

#For FEMALES
idxF<-which(lower.tri(commonedgesF)&commonedgesF==0,arr.ind=T)
paste(rownames(commonedgesF)[idxF[,"row"]],colnames(commonedgesF)[idxF[,"col"]],sep="--") 
nrow(idxF)/length(which(upper.tri(Clist_pcitF[[1]])))
#9/91=9.9% not present in any pop


#How many edges in all networks?
idx2<-which(lower.tri(commonedges)&commonedges==8,arr.ind=T)
paste(rownames(commonedges)[idx2[,"row"]],colnames(commonedges)[idx2[,"col"]],sep="--") 
#6/182=3.3% are present in all pops

# 
# col.indices<-apply(commonedges,1:2,function(x){max(which(x>=colscale))})
# j<-col.indices[upper.tri(col.indices)]
# k<-j[-which(is.infinite(j) )] #get rid of zero/undefined edges
# 
# pal<-colorRampPalette(c("black","gray90")) #palette ranging from (black, unique to 1 pop, to gray, in all pops)
# blacks<-pal(8)
# COLS<-blacks[k]
# #COLS<-COLS[which(!is.na(COLS))]
# 
# meanNet<-apply(simplify2array(Clist_pcit), 1:2, mean)
# diag(meanNet)<-NA
# length(Q(meanNet)$Edgelist[[1]])#79 edges
# 
# #test of colors
# plot(1:length(COLS),abs(as.vector(commonedges[which(upper.tri(commonedges)&commonedges!=0)])),col=COLS,pch=19,cex=3,xlab="edge index",ylab="# times in a network")
# 
# #This is the average
# dev.off()
# Q(meanNet,edge.color=COLS,title="Mean Phenotype Network",posCol=1,labels=toi2,shape=shps,vsize=8,vsize2=8)
# #Gray edges are super common; pure black edges are unique in 1 pop

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
################################
###### Get Z scores for traits
##
#MALES
rawmeansM<-aggregate(.~Pop,data=machos[,c("Pop",toi2)],mean)
#FEMALES
rawmeansF<-aggregate(.~Pop,data=hembras[,c("Pop",toi2)],mean)

# machos_sd_units<-apply(machos[,toi],2,function(x) (x/sd(x,na.rm=T)))#Put traits in sd units; not sure if necessary

#####
## convert to 0-1 scale
#MALES
Z_0to1M<-cbind(machos[,c("band","Sex","Pop","Year")],apply(machos[,toi2],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops

#FEMALES
Z_0to1F<-cbind(hembras[,c("band","Sex","Pop","Year")],apply(hembras[,toi2],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops

# #This CENTERS, and SCALES means
ZM<-cbind(machos[,c("band","Sex","Pop","Year")],apply(machos[,toi2],2,function(x) scale(x,center=T,scale=T)))
ZF<-cbind(hembras[,c("band","Sex","Pop","Year")],apply(hembras[,toi2],2,function(x) scale(x,center=T,scale=T)))

#MEANS by Country
#Z-scaled means for each trait in each pop
#MALES
(Zmeans_0to1M<-aggregate(.~Pop,data=Z_0to1M[,c("Pop",toi2)],mean)) 
#FEMALEs
(Zmeans_0to1F<-aggregate(.~Pop,data=Z_0to1F[,c("Pop",toi2)],mean)) 


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
ZmeltM<-melt(Z_0to1M,id=c("band","Sex","Pop","Year"))
ZmeltF<-melt(Z_0to1F,id=c("band","Sex","Pop","Year"))

ggplot(ZmeltM,aes(x=Pop,y=value))+geom_boxplot()+facet_wrap(~variable,nrow=2)+theme_bw()
ggsave("figs/Trait boxplots_BARS_males.jpeg",width=12)
ggplot(ZmeltF,aes(x=Pop,y=value))+geom_boxplot()+facet_wrap(~variable,nrow=2)+theme_bw()
ggsave("figs/Trait boxplots_BARS_females.jpeg",width=12)


#Make figure for Presentations
#Make custom theme for black background presentation in PPT
# PPTtheme<-theme(panel.background = element_rect(fill = "transparent", colour = "transparent"),
#     plot.background = element_rect(fill = "black",colour = "transparent"),
#     panel.grid.minor = element_blank(), 
#     panel.grid.major = element_blank(),#element_line(colour="white", size = 0.25),
#     axis.line=element_line(size=.25,colour="white"),
#     axis.text=element_text(size=14, colour="white"),
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     axis.title=element_text(size=18,face="bold",colour="white"),
#         strip.background = element_rect(fill = "transparent", color = "white", size = 1),
#     strip.text=element_text(colour="white",size=14))

#Test differences
for ( i in 1:14){
  trait<-toi2[i]
  print(trait)
  print(teval(paste0("summary(aov(",trait,"~Pop,data=ZF))")))
  
} #All sig. different


#new labls, changing Bri to Dk, bc of reversed order
toi2
toi3<-c("WL","TS","Tthet","TRach","TBri","Rthet","RRach","RBri","Bthet","BRach","BBri","Vthet","VRach","VBri")


#######
#-------------------------
# Plot Networks!!!
#MALES
pdf("figs/Evo presentation_8pop-pcit_males.pdf",width=16,height=8)
#define vertex size
vertsize<-12#apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
popnames<-sapply(levels(D3$Pop),function(x) substr(x,1,2))#c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: length(levels(D3$Pop))){
  mat<-Clist_pcitM[[i]]
  Q(mat,color="skyblue",border.color="gray20",labels=toi3,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="black",bg="transparent",label.color="black")
  mtext(popnames[i],3,line=.5,at=-1.,,adj=.5,col="black")
  mtext(paste0("n=",sum(!is.na(pcaListM[[i]]$scores[,1]))),side=3,line=-1.4,at=-1.,cex=.9,adj=.5,font=1,col="black")
}
dev.off()

#FEMALES
pdf("figs/Evo presentation_8pop-pcit_females.pdf",width=16,height=8)
#define vertex size
vertsize<-12#apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
popnames<-sapply(levels(D3$Pop),function(x) substr(x,1,2))#c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: length(levels(D3$Pop))){
  mat<-Clist_pcitF[[i]]
  Q(mat,color="salmon",border.color="gray20",labels=toi3,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="black",bg="transparent",label.color="black")
  mtext(popnames[i],3,line=.5,at=-1.,,adj=.5,col="black")
  mtext(paste0("n=",sum(!is.na(DlistF[[i]]$T_sum.B2)&!is.na(DlistF[[i]]$maxTS))),side=3,line=-1.4,at=-1.,cex=.9,adj=.5,font=1,col="black")

}
dev.off()

# Same plot, w throat and breast colored differently
pdf("figs/Evo presentation_8pop-pcit_R&T colored diff.pdf",width=16,height=8)
#define vertex size
vertsize<-12
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
nodecols<-c("white","white",rep("orangered",3),rep("orange",9))
popnames<-c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
for (i in 1: 8){
  mat<-Clist_pcit[[i]]
  Q(mat,color=nodecols,border.color="gray20",labels=toi2,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="white",bg="transparent",label.color="white")
  mtext(paste0("n=",sum(!is.na(pcaList[[i]]$scores[,1]))),side=1,line=-.75,at=1.1,cex=.9,font=1,col="white")
  mtext(popnames[i],3,line=-.6,at=-.8,,adj=.5,col="white")
}
dev.off()

#######
###
# Calculate Body color Integration Index
bodytraits<-toi2[-c(1:5)]
#MALES
bodyM<-lapply(ClistM,function(x) {tmp<-x[bodytraits,bodytraits] #Get corr matrix for body traits; set diag to NA
  diag(tmp)<-NA
  return(tmp)})
nonbodyM<-lapply(ClistM,function(x) {x[bodytraits,bodytraits]<-0 #Get corr matrix for body traits; set diag to NA
  diag(x)<-NA
  return(x)})
#FEMALES
bodyF<-lapply(ClistF,function(x) {tmp<-x[bodytraits,bodytraits] #Get corr matrix for body traits; set diag to NA
  diag(tmp)<-NA
  return(tmp)})
nonbodyF<-lapply(ClistF,function(x) {x[bodytraits,bodytraits]<-0 #Get corr matrix for body traits; set diag to NA
  diag(x)<-NA
  return(x)})

## MALES
#Calculate |density| of whole network for each population
pop.netdensityM<-sapply(ClistM,function(x) sum(abs(x),na.rm=T)/2/91)
#Calculate |density| of this cluster for each
pop.clusterdensityM<-sapply(bodyM,function(x) (sum(abs(x),na.rm=T)/2)/36)
#Calculate similar index, but for noncluster
pop.nonclusterdensityM<-sapply(nonbodyM,function(x) (sum(abs(x),na.rm=T)/2)/55)
#Calculate ratio of cluster to network density (integration index)
IIM<-pop.clusterdensityM/pop.netdensityM

## FEMALES
#Calculate |density| of whole network for each population
pop.netdensityF<-sapply(ClistF,function(x) sum(abs(x),na.rm=T)/2/91)
#Calculate |density| of this cluster for each
pop.clusterdensityF<-sapply(bodyF,function(x) (sum(abs(x),na.rm=T)/2)/36)
#Calculate similar index, but for noncluster
pop.nonclusterdensityF<-sapply(nonbodyF,function(x) (sum(abs(x),na.rm=T)/2)/55)
#Calculate ratio of cluster to network density (integration index)
IIF<-pop.clusterdensityF/pop.netdensityF

### Calculate overall body brightness index
#MALES
meanbriM<-sapply(DlistM,function(x) {df<-x[,c("R_sum.B2","B_sum.B2","V_sum.B2")]
            return(mean(colMeans(df,na.rm=T)))  })
IIdfM<-data.frame(IIM,pop.netdensityM,pop.clusterdensityM,pop.nonclusterdensityM,meanbriM,Pop=names(IIM))
#FEMALES
meanbriF<-sapply(DlistF,function(x) {df<-x[,c("R_sum.B2","B_sum.B2","V_sum.B2")]
            return(mean(colMeans(df,na.rm=T)))  })
IIdfF<-data.frame(IIF,pop.netdensityF,pop.clusterdensityF,pop.nonclusterdensityF,meanbriF,Pop=names(IIF))

jitterx<-meanbriM+.009
jitterx[5]<-jitterx[5]+.0012
jitterx[3]<-jitterx[3]-.02
jittery<-pop.clusterdensityM+.03
jittery[6]<-jittery[6]-.05

#MALES
g1<-ggplot(data=IIdfM,aes(x=meanbriM,y=pop.clusterdensityM,label=popnames))+geom_point(size=2,aes(col=meanbriM))+geom_text(size=5,col="black",x=jitterx,y=jittery)+xlab("")+ylab("Integration Index")+stat_ellipse(col="gray")+theme_bw()+scale_colour_gradient(limits=c(20, 40), low="#CC6600", high="#FFFFCC",guide=F)+ggtitle("Males")#+annotate("text",x=12,y=.2,adj=0,label="rho= .809, p= 0.022",col="black",size=5) #paste0("rho == ~-.833~ p== 0.015"),parse=T,col="white")
cor.test(meanbriM,pop.clusterdensityM,method="s")
ggsave("figs/Color Integ Index~Body Bri_males.jpg")

jitterxF<-meanbriF+.009
jitterxF[5]<-jitterxF[5]+.0012
jitterxF[3]<-jitterxF[3]-.02
jitteryF<-pop.clusterdensityF+.03
jitteryF[6]<-jitteryF[6]-.05
#FEMALES
g2<-ggplot(data=IIdfF,aes(x=meanbriF,y=pop.clusterdensityF,label=popnames))+geom_point(size=2,aes(col=meanbriM))+geom_text(size=5,col="black",x=jitterxF,y=jitteryF)+xlab("Body Brightness Index")+ylab("Integration Index")+stat_ellipse(col="gray")+theme_bw()+scale_colour_gradient(limits=c(20, 40), low="#CC6600", high="#FFFFCC",guide=F)+ggtitle("Females")#+annotate("text",x=12,y=.2,adj=0,label="rho= .809, p= 0.022",col="black",size=5) #paste0("rho == ~-.833~ p== 0.015"),parse=T,col="white")
cor.test(meanbriM,pop.clusterdensityM,method="s")
ggsave("figs/Color Integ Index~Body Bri_females.jpg")

require(grid)
jpeg("figs/Color Integ Index~Body Bri_both.jpg")
multiplot(g1,g2)
dev.off()



# ggplot(IIdf,aes(x=meansumbri,y=pop.nonclusterdensity,label=Pop))+geom_point(size=2,aes(col=meansumbri))+geom_text(size=5,col="white",aes(x=jitter(meansumbri,10),y=pop.nonclusterdensity+.03))+xlab("Body Brightness Index")+ylab("Integration Index")+stat_ellipse(col="white")+PPTtheme+ylim(c(0,.7))+scale_colour_gradient(limits=c(70, 120), high="#FFFFCC", low="#CC6600",guide=F)+annotate("text",x=130,y=.05,label="rho= 0.429, p= 0.299",col="white",size=5) 
# cor.test(meansumbri,pop.nonclusterdensity,method="s")
# ggsave("figs/Evolution Presentation_Color Integ Index~NonBody Brightness.jpg")


#############
##### Plot networks ordered by modularity
#reorder factor levels west to east
df.comb$Pop2<-factor(df.comb$Pop)
levels(df.comb$Pop2)=rank(pop.clusterdensity)
df.comb$Pop2<-factor(df.comb$Pop2,levels=8:1)



pdf("figs/Evo presentation_8pop-pcit_ordered by clustering.pdf",width=16,height=8)
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
dev.off()

require(car)
leveneTest(R_Chrom~Pop,machos)
boxplot(R_Chrom~Pop,machos)
tapply(machos$R_Chrom,machos$Pop,sd,na.rm=T)