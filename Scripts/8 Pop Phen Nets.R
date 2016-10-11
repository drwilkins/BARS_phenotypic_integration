#git config --global user.email "mrwilkins06@gmail.com"

require(igraph);require(PCIT);require(qgraph);require(psych);require(plotrix);require(assortnet)
require(devtools) #necessary to import files from dropbox over https

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
pcaM<-principal(DlistM[[i]][,toi2],rot="none",nfactors=length(toi2)-2)#Won't work for full set, calc 2 fewer eigenvectors
pcaF<-principal(DlistF[[i]][,toi2],rot="none",nfactors=length(toi2)-2)
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
pcaListF<-lapply(1:length(levels(D3$Pop)),function(x) principal(DlistF[[x]][,toi2],rot="varimax",nfactors=extractF[x]))
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
rawmeans<-aggregate(.~Pop,data=machos[,c("Pop",toi)],mean)

# machos_sd_units<-apply(machos[,toi],2,function(x) (x/sd(x,na.rm=T)))#Put traits in sd units; not sure if necessary
# convert to 0-1 scale
Z_0to1<-cbind(machos[,c("ID","Sex","Pop","Year")],apply(machos[,toi],2,function(x) ((x-min(x,na.rm=T))/diff(range(x,na.rm=T)))))#scales traits of interest across all pops
#This CENTERS, and SCALES means
Z<-cbind(machos[,c("ID","Sex","Pop","Year")],apply(machos[,toi],2,function(x) scale(x,center=T,scale=T)))

#MEANS by Country
(Zmeans_0to1<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],mean)) #Z-scaled means for each trait in each pop
Zmeans_0to1_revbri<-Zmeans_0to1
Zmeans_0to1_revbri[,c(4,7,10,13)]<-apply(Zmeans_0to1_revbri[,c(4,7,10,13)],2,function(x) (1-x))

#calculate coefficients of variation (sd/mean)
CV<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],function(x) (sd(x,na.rm=T)/mean(x,na.rm=T)))
SD<-aggregate(.~Pop,data=Z_0to1[,c("Pop",toi)],function(x) (sd(x,na.rm=T)))


#Plot Z-trait boxplots (* Note Z_0to1 are not normal Z-scales, they're bounded 0:1)
require(reshape)
Zmelt<-melt(Z_0to1,id=c("ID","Sex","Pop","Year"))
ggplot(Zmelt,aes(x=Pop,y=value))+geom_boxplot()+facet_wrap(~variable,nrow=2)+theme_bw()
ggsave("figs/Trait boxplots_BARS.jpeg",width=12)

#Make figure for Presentations
#Make custom theme for black background presentation in PPT
PPTtheme<-theme(panel.background = element_rect(fill = "transparent", colour = "transparent"),
    plot.background = element_rect(fill = "black",colour = "transparent"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),#element_line(colour="white", size = 0.25),
    axis.line=element_line(size=.25,colour="white"),
    axis.text=element_text(size=14, colour="white"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title=element_text(size=18,face="bold",colour="white"),
        strip.background = element_rect(fill = "transparent", color = "white", size = 1),
    strip.text=element_text(colour="white",size=14))

ggplot(Zmelt,aes(x=Pop,y=value))+geom_boxplot(col="white")+facet_wrap(~variable,nrow=2)+PPTtheme+xlab("\nPopulation")+ylab("Scaled Trait Value\n")+theme(axis.title.y=element_text(vjust=.5))
ggsave("figs/Evo Presentation_Trait boxplots_BARS.png",width=12)

#Test differences
for ( i in 1:14){
  trait<-toi[i]
  print(trait)
  print(teval(paste0("summary(aov(",trait,"~Pop,data=Z))")))
  
} #All sig. different


#new labls, changing Bri to Dk, bc of reversed order
toi3<-c("RWL","RTS","TDk","THue","TChr","RDk","RHue","RChr","BDk","BHue","BChr","VDk","VHue","VChr")

###########################################
#-------------------------------------------------------
#now plot the networks w/o arcs, relative to Colorado
pdf("figs/8pop-pcit filtered_gray_f-r.pdf",width=16,height=8)
#define vertex size
vertsize<-apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
for (i in 1: 8){
  mat<-Clist_pcit[[i]]
  j2<-col.indices*(apply(mat,1:2,function(x) ifelse(x==0,0,1))) #multiply col.index matrix by matrix at hand (to remove zero edges)
  j2<-j2[upper.tri(j2)] #get rid of bottom half
k2<-j2[-which(is.infinite(j) )] #get rid of zero/undefined edges
pal<-colorRampPalette(c("black","gray90")) #palette ranging from (black, unique to 1 pop, to gray, in all pops)
blacks<-pal(8)
icol<-blacks[k2]

    Q(mat,title=names(Clist_pcit)[i],posCol=1,labels=toi3,shape=shps,vsize=vertsize[i,],vsize2=vertsize[i,],edge.color=icol)
  mtext(paste0("n=",sum(!is.na(pcaList[[i]]$scores[,1]))),side=1,line=-.75,at=1.1,cex=.9,font=1)
}
dev.off()

#######
#-------------------------
# # Plot for presentations
# pdf("figs/Evo presentation_8pop-pcit.pdf",width=16,height=8)
# #define vertex size
# vertsize<-8#apply(Zmeans_0to1_revbri[,-1],2,function(x)as.numeric(15*x))+8
# par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
# popnames<-c("Colorado","New York","UK","Czech Rep","Romania","Turkey","Israel","Taiwan")
# for (i in 1: 8){
#   mat<-Clist_pcit[[i]]
#   Q(mat,color="white",border.color="gray20",labels=toi3,shape=shps,vsize=vertsize,vsize2=vertsize,edge.color="white",bg="transparent",label.color="white")
#   mtext(paste0("n=",sum(!is.na(pcaList[[i]]$scores[,1]))),side=1,line=-.75,at=1.1,cex=.9,font=1,col="white")
#   mtext(popnames[i],3,line=-.6,at=-.8,,adj=.5,col="white")
# }
# dev.off()


#####---------------------
## Test assortativity by size
assortbyZ<-list(0)
for(i in 1:length(Clist_pcit)) {
  mat<-Clist_pcit[[i]]
  diag(mat)<-0
  assortbyZ[[i]]<-assortment.continuous(abs(mat),Zmeans_0to1_revbri[i,-1],weighted=T)
}
assortbyZvec<-unlist(assortbyZ)
names(assortbyZvec)<-names(Clist_pcit)
assortbyZvec
assortment.continuous(mat,rnorm(20,1,1),weighted=T)
#--------------------------------------------
##  Make networks using igraph to calculate stats and whatnot

iCO<-graph_from_adjacency_matrix(abs(Clist_pcit$CO),"undirected",weighted=T,diag=F)
plot(iCO,vertex.size=40,layout=layout.circle,vertex.shape=shps)
iIA<-graph_from_adjacency_matrix(abs(Clist_pcit$IA),"undirected",weighted=T,diag=F)
iUK<-graph_from_adjacency_matrix(abs(Clist_pcit$UK),"undirected",weighted=T,diag=F)
iCR<-graph_from_adjacency_matrix(abs(Clist_pcit$CR),"undirected",weighted=T,diag=F)
iRM<-graph_from_adjacency_matrix(abs(Clist_pcit$RM),"undirected",weighted=T,diag=F)
iTR<-graph_from_adjacency_matrix(abs(Clist_pcit$TR),"undirected",weighted=T,diag=F)
iIL<-graph_from_adjacency_matrix(abs(Clist_pcit$IL),"undirected",weighted=T,diag=F)
iTW<-graph_from_adjacency_matrix(abs(Clist_pcit$TW),"undirected",weighted=T,diag=F)

items<-c("iCO","iIA","iUK","iCR","iRM","iTR","iIL","iTW") #vector of adjacency matrix names to iterate thru


#****************************************************************************
#****************************************************************************
####  Test relationship bw population pairwise diff in degree vs diff in mean Z-score
####   ...we predict a positive relationship between amount of divergence in means & divergence in corr structure

  # Extract degrees for each node for each population
degrees<-lapply(items,function(x) degree(teval(x)))
(degrees<-do.call(rbind,degrees))
viz(degrees[2:15])
deg.diff<-data.frame(Dyads=Z.meandiff$Dyads,apply(degrees,2,function(trait) apply(comparison,2,function(dyad) trait[dyad][1]-trait[dyad][2] )))
deg.diff



###########
## OK, now get clustering coeff diffs
  # Extract transitivity for each node for each population
trans<-lapply(items,function(x) transitivity(teval(x),"localundirected",isolates="zero"))
(trans<-do.call(rbind,trans))
viz
colnames(trans)<-names(V(teval(items[1])))


#output Z & degree diffs
degrees<-as.data.frame(degrees)
degrees$var<-"deg"
degrees$Pop<-Zmeans_0to1$Pop
degrees<-degrees[,c(16,1:15)]
#Z.meandiff$var<-"traitmean"
Zmeans_0to1$var<-"traitmean"
trans<-as.data.frame(trans)
trans$var<-"transitivity"
trans$Pop<-Zmeans_0to1$Pop
trans<-trans[,c(16,1:15)]
rawmeans$var<-"rawmean"
CV$var<-"CV"
SD$var<-"SD"

######### Just plot Degree & Transitivity as f(Z trait Mean)
require(reshape2)
df2<-rbind(Zmeans_0to1,degrees,trans,CV,SD,rawmeans)
write.csv(df2,"data/Trait means, C, & degree.csv")
df_melt<-melt(df2)
#dcast(melt(df2),Pop+var+variable~value,id.var=c("var","Pop","variable"))

## Graph Degree ~ Trait means for all traits
tmp<-subset(df_melt,var=="traitmean"|var=="deg")
df_deg<-dcast(tmp,variable+Pop~var)
degxZcor<-t(sapply(levels(df_deg$variable),function(x) {foo<-subset(df_deg,variable==x)
  results<-cor.test(foo$deg,foo$traitmean,method="s")
  return(c(results$estimate,results$p.value))})) #get correlation bw degree & trait mean for each trait
colnames(degxZcor)<-c("rho","p.val")
degxZcor<-data.frame(degxZcor,p.adj=p.adjust(degxZcor[,"p.val"],"fdr"))
degxZcor$sigcode<-ifelse(degxZcor$p.adj<=.05,2,1)
degxZcor

ggplot(df_deg,aes(x=traitmean,y=deg))+geom_point()+facet_wrap(~variable,nrow=2)+theme_bw()+annotate("text",x=1.2,y=0,label=paste0("rho=",round(degxZcor$rho,3)),adj=1,fontface=degxZcor$sigcode)+ylab("Node Degree")+xlab("Node Mean (Z Score)")

ggsave("figs/Degree~traitmeans.jpeg")

## Graph transitivity ~ Trait means for all traits
tmp<-subset(df_melt,var!="deg")
df_trans<-dcast(tmp,variable+Pop~var)
transxZcor<-t(sapply(levels(df_trans$variable),function(x) {foo<-subset(df_trans,variable==x)
  results2<-cor.test(foo$transitivity,foo$traitmean,method="s")
  return(c(results2$estimate,results2$p.value))})) #get correlation bw degree & trait mean for each trait
colnames(transxZcor)<-c("rho","p.val")
transxZcor<-data.frame(transxZcor,p.adj=p.adjust(transxZcor[,"p.val"],"fdr"))
transxZcor$sigcode<-ifelse(transxZcor$p.val<=.05,2,1) #not using p.adj value! Should, but need other populations
transxZcor

ggplot(df_trans,aes(x=traitmean,y=transitivity))+geom_point()+facet_wrap(~variable,nrow=2)+theme_bw()+annotate("text",x=1.2,y=0,label=paste0("rho= ",round(transxZcor$rho,3)),adj=1,fontface=transxZcor$sigcode)+ylab("Node Clustering Coefficient")+xlab("Node Mean (Z Score)")

ggsave("figs/C~traitmeans.jpeg")




#Plot transitivty~CV
xypanels<-function(X,Y,meltedData,color) #modified for presentation purposes
{  
  if(missing(color)){colorstring<-""}else{colorstring<-paste0(",col=",color)}
  xycor<-t(sapply(levels(meltedData$variable),function(i) 
    {
    SUB<-subset(meltedData,variable==i)   #pull out measures for each trait in turn (e.g. femur dark)
    x2<-eval(parse(text=paste0("SUB$",X)),envir=SUB)
    y2<-eval(parse(text=paste0("SUB$",Y)),envir=SUB)
    results<-cor.test(x2,y2,method="s")
    return(c(results$estimate,results$p.value)) 
    })) #get correlation bw Measures for all traits
  colnames(xycor)<-c("rho","p.val")
  xycor<-data.frame(xycor,p.adj=p.adjust(xycor[,"p.val"],"fdr"))
  xycor$sigcode<-ifelse(xycor$p.val<=.05,2,1) #not using p.adj value! Should, but need other populations
  maxx<-max(eval(parse(text=paste0("meltedData$",X))),na.rm=T)
  miny<-min(eval(parse(text=paste0("meltedData$",Y))),na.rm=T)
  graf<-
    ggplot(meltedData,eval(parse(text=paste0("aes(x=",X,",y=",Y,colorstring,")"))))+PPTtheme+stat_ellipse(alpha=.25,col="white")+
     geom_point(col="white")+facet_wrap(~variable,nrow=2)+annotate("text",x=maxx,y=miny-.4,col="white",label=paste("rho== ",round(xycor$rho,3)),parse=T,adj=1,fontface=xycor$sigcode)+ylab(Y)+xlab(X)
    resultado<-list(xycor,graf)
  plot(graf)
  names(resultado)<-c("xycor","graf")
  return(resultado)
  }#End xypanels

plot(xypanels("traitmean","transitivity",dcast(df_melt,variable+Pop~var))$graf+ylab("Node Clustering Coef")+xlab("Scaled Node Mean"))

ggsave("figs/Evolution Presentation_C~traitmeans.jpeg",width=12)

plot(xypanels("CV","transitivity",dcast(df_melt,variable+Pop~var))$graf+ylab("Node Clustering Coef")+xlab("Node Coefficient of Variation"))
ggsave("figs/Evolution presentation_trans~CV.jpeg",width=10,height=5)

plot(xypanels("SD","transitivity",dcast(df_melt,variable+Pop~var))$graf+ylab("Node Clustering Coef")+xlab("Node Standard Deviation"))
ggsave("figs/Evolution presentation_trans~SD.jpeg",width=10,height=5)

# Trait Means (scaled 1:0)~ Trait SD
gZmeanSD<-xypanels("SD","traitmean",dcast(df_melt,variable+Pop~var))
ggsave("figs/Evolution presentation_Mean~SD.jpeg",width=10,height=5)
#Test strongest relationship
cor.test(Zmeans_0to1$B_Chrom,SD$B_Chrom,method="s")
cor.test(Zmeans_0to1$B_AvgBri,SD$B_AvgBri,method="s")
ggplot(aes(x=B_AvgBri,y=B_Chrom),data=Z_0to1,groups=Pop)+geom_point(aes(groups="mean"))

gRawmeanSD<-xypanels("SD","rawmean",dcast(df_melt,variable+Pop~var))
plot(SD$B_Chrom,Zmeans_0to1$B_Chrom)
plot(SD$B_AvgBri,Zmeans_0to1$B_AvgBri) #Doesn't make any fucking sense! Should be opposite
plot(SD$B_AvgBri,SD$B_Chrom)
plot(Zmeans_0to1$B_AvgBri,Zmeans_0to1$B_Chrom)



