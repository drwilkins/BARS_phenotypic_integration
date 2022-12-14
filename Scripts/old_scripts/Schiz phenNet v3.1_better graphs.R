#Schiz phen network script v3

require(qgraph);require(PCIT);require(igraph); require(reshape2); require(ggplot2); require(grid); require(devtools)
source_url("https://dl.dropboxusercontent.com/s/b32xvn9x18gahc6/bioworkflow.R") #bioworkflow.R script
#>>----------------------------------------------------------->>
# Begin function definitions
#>>----------------------------------------------------------->>
### Make a function that will do the subsetting, correlation test, and everything, and plot the points for F & M on the same graph.
###### (1)
#Then Look at Node Strength ~ Node transitivity plots
xypanels<-function(X,Y,meltedData,color)
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
    ggplot(meltedData,eval(parse(text=paste0("aes(x=",X,",y=",Y,colorstring,")"))))+stat_ellipse(alpha=.25)+
     geom_point()+facet_wrap(~variable,nrow=2)+theme_bw()+annotate("text",x=maxx,y=miny,label=paste("rho== ",round(xycor$rho,3)),parse=T,adj=1,fontface=xycor$sigcode)+ylab(Y)+xlab(X)
    resultado<-list(xycor,graf)
  plot(graf)
  names(resultado)<-c("xycor","graf")
  return(resultado)
  }

######### (2)
## modified version of xypanels for plotting correlations for m & fem
xypanels.sexes<-function(X,Y,meltedData,color)
{  
  if(missing(color)){colorstring<-""}else{colorstring<-paste0(",col=",color)} #flexi code for incorporating color (sex) variable or not
 #function for getting trait correlations for fem and male datasets
   xycor.fun<-function(i,data,X,Y){
    SUB<-subset(data,variable==i)   #pull out measures for each trait in turn (e.g. femur dark)
    results<-cor.test(eval(parse(text=paste0("SUB$",X)),envir=SUB),eval(parse(text=paste0("SUB$",Y)),envir=SUB),method="s")
    return(c(results$estimate,results$p.value)) 
    } #get correlation bw Measures for all traits
 #Get fem correlations
    xycor.f<-t(sapply(levels(meltedData$variable),function(i) {xycor.fun(i,f.recast,X,Y)}))
  colnames(xycor.f)<-c("rho","p.val")
  xycor.f<-data.frame(xycor.f,p.adj=p.adjust(xycor.f[,"p.val"],"fdr"))
  xycor.f$sigcode<-ifelse(xycor.f$p.val<=.05,2,1) #not using p.adj value! Should, but need other populations
  maxx.f<-max(eval(parse(text=paste0("meltedData$",X))),na.rm=T)
  yvec.f<-eval(parse(text=paste0("f.recast$",Y)))
  maxy.f<-max(yvec.f,na.rm=T)+sd(yvec.f)/2

    #Get male correlations
    xycor.m<-t(sapply(levels(meltedData$variable),function(i) {xycor.fun(i,m.recast,X,Y)}))
  colnames(xycor.m)<-c("rho","p.val")
  xycor.m<-data.frame(xycor.m,p.adj=p.adjust(xycor.m[,"p.val"],"fdr"))
  xycor.m$sigcode<-ifelse(xycor.m$p.val<=.05,2,1) #not using p.adj value!
  maxx.m<-max(eval(parse(text=paste0("meltedData$",X))),na.rm=T)
  yvec.m<-eval(parse(text=paste0("m.recast$",Y)))
  miny.m<-min(yvec.m,na.rm=T)-sd(yvec.m)/2
  #extract colors from dummy graph
p<-ggplot(meltedData, eval(parse(text=paste0("aes(x=",X,",y=",Y,colorstring,")")))) + geom_point(shape=21, size=4)
ggcols<-unique(ggplot_build(p)$data[[1]]$colour)
  
  graf<-
    ggplot(meltedData,eval(parse(text=paste0("aes(x=",X,",y=",Y,colorstring,")"))))+stat_ellipse(aes_string(col=color,fill=color),alpha=.25)+
     geom_point()+facet_wrap(~variable,nrow=2)+theme_bw()+
  annotate("text",x=maxx.f,y=maxy.f,label=paste("rho== ",round(xycor.f$rho,3)),parse=T,adj=1,size=3+.75*xycor.f$sigcode,col=ggcols[1])+#Add Female value
   annotate("text",x=maxx.m,y=miny.m,label=paste("rho== ",round(xycor.m$rho,3)),parse=T,adj=1,size=3+.75*xycor.m$sigcode,col=ggcols[2])+#Add Male value
  ylab(Y)+xlab(X)
    resultado<-list(xycor.m,xycor.f,graf)
  plot(graf)
  names(resultado)<-c("xycor.m","xycor.f","graf")
  return(resultado)
  }
#<<-----------------------------------------------------------<<
# End of function definitions
#<<-----------------------------------------------------------<<


#read in data
schiz<-read.csv("/Users/mattwilkins/dropbox/My Research/My Papers/Alex Hansen Project/Data/schiz data (corrected).csv")
#Get rid of unnecessary columns in data (redefine variable)
names(schiz)#shows column names in data frame
schiz<-schiz[,-c(12:16,18)]
check_class(schiz)
schiz[,2:16]<-apply(schiz[,2:16],2,FUN=as.numeric)#make all measured columns numeric

#Get means for all traits
schiz.summary<-sapply(names(schiz)[-c(1,17,18)],function(x){aggregate(formula(paste(x,"~spp+sex")),FUN=function(x) {mean(x)},data=schiz)},simplify=F)

schiz.table<-as.data.frame(schiz.summary[-c(14,15)],check.names=T)
schiz.table<-schiz.table[,c(1:3,seq(6,39,3))] #removes extra spp columns
names(schiz.table)[2]<-"sex"
schiz.table


#Make list w/ each spp+sex as a separate dataframe
#females
f.schiz.DFs<-lapply(levels(schiz$spp),function(x) (subset(schiz,spp==x&sex=="F")))
names(f.schiz.DFs)<-levels(schiz$spp)
#males
m.schiz.DFs<-lapply(levels(schiz$spp),function(x) (subset(schiz,spp==x&sex=="M")))
names(m.schiz.DFs)<-levels(schiz$spp)

#calculate Z scores and means for traits of interest
TOI<-c('femur.mean.dark','femur.pixel.area','patella.mean.dark','patella.pixel.area','tibia.mean.dark','tibia.pixel.area','meta.mean.dark','meta.pixel.area','tarsus.mean.dark','tarsus.pixel.area','mean.ceph.width','femur.length','tibia.length')
# m.Z<-lapply(m.schiz.DFs,function(x) apply(x[,TOI],2,scale))
# #make sure the Z transformation worked properly
# plot(m.schiz.DFs[[1]][,"tibia.length"],m.Z[[1]][,"tibia.length"])

#get Z scores for all spp & sexes
Z<-as.data.frame(apply(schiz[,TOI],2,scale))
Z$spp<-schiz$spp
Z$sex<-schiz$sex
f.Zmeans<-aggregate(.~spp,data=subset(Z,sex=="F",select=c("spp",TOI)),FUN=mean,na.rm=T)
m.Zmeans<-aggregate(.~spp,data=subset(Z,sex=="M",select=c("spp",TOI)),FUN=mean,na.rm=T)
#check that these Z scores match untransformed means
plot(subset(schiz.table,sex=="M")$femur.mean.dark.femur.mean.dark,m.Zmeans$femur.mean.dark)#yep
plot(subset(schiz.table,sex=="F")$femur.mean.dark.femur.mean.dark,f.Zmeans$femur.mean.dark)  #yep

#Make list of correlations for each species (for traits of interest)
f.schiz.corrs<-lapply(levels(schiz$spp), function(x) cor(as.matrix(f.schiz.DFs[[x]][,TOI]),method="s",use="pairwise.complete"))
names(f.schiz.corrs)<-levels(schiz$spp)

m.schiz.corrs<-lapply(levels(schiz$spp), function(x) cor(as.matrix(m.schiz.DFs[[x]][,TOI]),method="s",use="pairwise.complete"))
names(m.schiz.corrs)<-levels(schiz$spp)

##### Filter ze networks!!! (PCIT)
# Get list of PCIT-filtered correlation matrices

#females
f.schiz_pcit<-lapply(levels(schiz$spp),function(x) {
  old<-f.schiz.corrs[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(f.schiz_pcit)<-levels(schiz$spp)
#males
m.schiz_pcit<-lapply(levels(schiz$spp),function(x) {
  old<-m.schiz.corrs[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(m.schiz_pcit)<-levels(schiz$spp)

## Make custom plot function
Q<-function(COR,...){
  G<-qgraph(COR,diag=F,fade=F,label.color="black",label.font=2,label.scale=T,label.norm="0000",negCol="black",layout="circle",...)
return(G)}

### Plot filtered outputs

#fems
par(mfrow=c(2,7))
for (i in 1: length(f.schiz_pcit)){
  Q(f.schiz_pcit[[i]],title=paste0(names(f.schiz_pcit)[i],"_Fems"))
}
#males
par(mfrow=c(2,7))
for (i in 1: length(f.schiz_pcit)){
  Q(m.schiz_pcit[[i]],title=paste0(names(m.schiz_pcit)[i],"_Males"))
}

####
# Which correlations are in all networks? And what is the average value?
#Make female binary network
f.schiz_pcit_bin<-lapply(1:length(f.schiz_pcit),function(x) apply(f.schiz_pcit[[x]],1, function(xx) ifelse(xx==0,0,1)))#binary correlation matrices (as list)
names(f.schiz_pcit_bin)<-names(f.schiz_pcit)
#males
m.schiz_pcit_bin<-lapply(1:length(m.schiz_pcit),function(x) apply(m.schiz_pcit[[x]],1, function(xx) ifelse(xx==0,0,1)))#binary correlation matrices (as list)
names(m.schiz_pcit_bin)<-names(m.schiz_pcit)

(f.commonedges<-teval(paste0("f.schiz_pcit_bin[[",1:length(f.schiz_pcit),"]]",collapse="+"))) #Adds all spp binary edges together
(m.commonedges<-teval(paste0("m.schiz_pcit_bin[[",1:length(f.schiz_pcit),"]]",collapse="+")))
diag(f.commonedges)<-NA
diag(m.commonedges)<-NA
colscale<-seq(1,length(f.schiz_pcit),1) #color scale goes from 1 to n; n if edge present in all pops (gray); 1 if only present once (black)

#Which and how many edges don't exist?
#fems
f.idx<-which(lower.tri(f.commonedges)&f.commonedges==0,arr.ind=T)
paste(rownames(f.commonedges)[f.idx[,"row"]],colnames(f.commonedges)[f.idx[,"col"]],sep="--") 
sum(!is.na(f.commonedges)) #how many possible edges
#1/156=0.6% not present in any pop

#males
m.idx<-which(lower.tri(m.commonedges)&m.commonedges==0,arr.ind=T) #none; all edges exist in some spp

#How many edges in all networks?
#fems
f.idx2<-which(lower.tri(f.commonedges)&f.commonedges==length(f.schiz_pcit),arr.ind=T)
paste(rownames(f.commonedges)[f.idx2[,"row"]],colnames(f.commonedges)[f.idx2[,"col"]],sep="--") 
#5/156=3.2%  present in all spp

#males
m.idx2<-which(lower.tri(m.commonedges)&m.commonedges==length(m.schiz_pcit),arr.ind=T)
paste(rownames(m.commonedges)[m.idx2[,"row"]],colnames(m.commonedges)[m.idx2[,"col"]],sep="--") 
#5/156=3.2%  present in all spp (but only tibia L--Femur Length shared across sexes)

pal<-colorRampPalette(c("black","gray90")) #palette ranging from (black, unique to 1 pop, to gray, in all pops)
blacks<-pal(length(m.schiz_pcit))

#now plot the networks 

#succinct node labels
Labels<-c("FemurMD","FemurPA","PatMD","PatPA","TibiaMD","TibiaPA","MetaMD","MetaPA","TarsMD","TarsPA","CephWid","FemurL","TibiaL")
# shapes
schizshapes<-rep(NA,length(Labels))
schizshapes[c(1,3,5,7,9)]<-"triangle"
schizshapes[c(2,4,6,8,10:13)]<-"rectangle"

#Make PDF of female networks
pdf("figs/schiz females filtered_gray.pdf",width=16,height=6)
par(mfrow=c(2,7),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
for (i in 1: length(f.schiz_pcit)){
  mat<-f.schiz_pcit[[i]]
  colindx.i<-f.commonedges*(apply(mat,1:2,function(x) ifelse(x==0,0,1))) #multiply col.index matrix by matrix at hand (to remove zero edges)
  colindx.i<-colindx.i[which(upper.tri(colindx.i)&!is.na(colindx.i)&colindx.i>0)] #get rid of bottom half & zero/undefined edges
  icol<-blacks[colindx.i]

 Q(mat,title=names(f.schiz_pcit)[i],posCol=1,labels=Labels,shape=schizshapes,vsize=10,vsize2=10,edge.color=icol)
 mtext(paste0("n=",nrow(f.schiz.DFs[[i]])),side=1,line=-.75,at=.9,cex=.8,font=1)
}
dev.off()

#Make PDF of Male networks
pdf("figs/schiz males filtered_gray.pdf",width=16,height=6)
par(mfrow=c(2,7),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
for (i in 1: length(m.schiz_pcit)){
  mat<-m.schiz_pcit[[i]]
  colindx.i<-m.commonedges*(apply(mat,1:2,function(x) ifelse(x==0,0,1))) #multiply col.index matrix by matrix at hand (to remove zero edges)
  colindx.i<-colindx.i[which(upper.tri(colindx.i)&!is.na(colindx.i)&colindx.i>0)] #get rid of bottom half & zero/undefined edges
  icol<-blacks[colindx.i]

 Q(mat,title=names(m.schiz_pcit)[i],posCol=1,labels=Labels,shape=schizshapes,vsize=10,vsize2=10,edge.color=icol)
 mtext(paste0("n=",nrow(m.schiz.DFs[[i]])),side=1,line=-.75,at=.9,cex=.8,font=1)
}
dev.off()


#--------------------------------------------
# Make networks using igraph to calculate stats and whatnot
levels(schiz$spp)
#FEMALES
f.igraf<-lapply(f.schiz_pcit, function(x) graph_from_adjacency_matrix(abs(x),"undirected",weighted=T,diag=F))
#MALES
m.igraf<-lapply(m.schiz_pcit, function(x) graph_from_adjacency_matrix(abs(x),"undirected",weighted=T,diag=F))

schiz.items<-c('avida','bilineata','crassipalpata','crassipes','duplex','floridana','maxima','mccooki','ocreata','retrorsa','rovneri','saltatrix','stridulans','uetzi') #vector of adjacency matrix names to iterate thru

  # Extract degrees for each node for each population
f.deg0<-lapply(schiz.items,function(x) degree(teval(paste0('f.igraf$',x))))
(f.deg<-as.data.frame(do.call(rbind,f.deg0)))
m.deg0<-lapply(schiz.items,function(x) degree(teval(paste0('m.igraf$',x))))
(m.deg<-as.data.frame(do.call(rbind,m.deg0)))

   # Extract transitivity for each node for each population
f.trans0<-lapply(schiz.items,function(x) transitivity(teval(paste0('f.igraf$',x)),"local",isolates="zero"))
(f.trans<-as.data.frame(do.call(rbind,f.trans0)))
m.trans0<-lapply(schiz.items,function(x) transitivity(teval(paste0('m.igraf$',x)),"local",isolates="zero"))
(m.trans<-as.data.frame(do.call(rbind,m.trans0)))
colnames(f.trans)<-colnames(m.trans)<-names(V(f.igraf$avida))
m.trans$spp<-f.trans$spp<-f.deg$spp<-m.deg$spp<-schiz.items
m.trans$var<-f.trans$var<-"transitivity"
m.deg$var<-f.deg$var<-"degree"
m.Zmeans$var<-f.Zmeans$var<-"traitmean"

 # Extract strength for each node for each population
f.stren0<-lapply(schiz.items,function(x) strength(teval(paste0('f.igraf$',x))))
(f.stren<-as.data.frame(do.call(rbind,f.stren0)))
m.stren0<-lapply(schiz.items,function(x) strength(teval(paste0('m.igraf$',x))))
(m.stren<-as.data.frame(do.call(rbind,m.stren0)))
colnames(f.stren)<-colnames(m.stren)<-names(V(f.igraf$avida))
f.stren$spp<-m.stren$spp<-schiz.items
f.stren$var<-m.stren$var<-"strength"

#########################
#Calculate dimorphism (in mean) for every trait
dimorph<-(m.Zmeans[,-c(1,15)]) - (f.Zmeans[,-c(1,15)])
dimorph$spp<-m.Zmeans$spp
dimorph$var<-"dimorphism" # negative means female more exaggerated

#########################
#Calculate dimorphism (in transitivity) for every trait
dimtrans<-m.trans[,1:13]-f.trans[,1:13]
dimtrans$spp<-m.trans$spp
dimtrans$var<-"dimtrans"

f.df<-rbind(f.trans,f.deg,f.Zmeans,f.stren,dimorph,dimtrans)
m.df<-rbind(m.trans,m.deg,m.Zmeans,m.stren,dimorph,dimtrans)
require(reshape2)
f.melt<-melt(f.df)
m.melt<-melt(m.df)
f.recast<-dcast(f.melt,variable+spp~var)
m.recast<-dcast(m.melt,variable+spp~var)
full.recast<-rbind(f.recast,m.recast)
full.recast$sex<-c(rep("F",nrow(f.recast)),rep("M",nrow(m.recast)))
######*******************************************************
## Graph transitivity ~ Trait means for all traits
xypanels.sexes("traitmean","transitivity",full.recast,"sex")
ggsave("figs/C~traitmeans_schiz.jpeg",width=10,height=5)


######*******************************************************
## Graph transitivity ~ Trait DIMORPHISM for all traits
xypanels.sexes("dimorphism","transitivity",full.recast,"sex")
ggsave("figs/C~dimorph_schiz.jpeg",width=10,height=5)

#################################
#plot dimorphism boxplots
ggplot(m.recast,aes(x=variable,y=dimorphism))+geom_boxplot()+theme(axis.text.x=element_text(angle=45,hjust=1))
#most traits are monomorphic in some spp. Exceptions metatarsus pixel area, tibia length, and ceph width

######*******************************************************
## Graph transitivity dimorphism ~ Trait DIMORPHISM for all traits
xypanels("dimorphism","dimtrans",full.recast)
ggsave("figs/Dimorphism of transitivity~dimorph of means_schiz.jpeg",width=10,height=5)

######*******************************************************
## Graph node strength ~ Trait Z-mean for all traits
xypanels.sexes("traitmean","strength",full.recast,"sex")
ggsave("figs/Node Strength~Node Mean_schiz.jpeg",width=10,height=5)


######*******************************************************
## Graph node strength ~ Dimorphism for all traits
xypanels.sexes("dimorphism","strength",full.recast,"sex")
ggsave("figs/Node strength~Node dimorphism_schiz.jpeg",width=10,height=5)

######*******************************************************
## Graph node degree ~ Trait Z-mean for all traits
xypanels.sexes("traitmean","degree",full.recast,"sex")
ggsave("figs/Node degree~Node Mean_schiz.jpeg",width=10,height=5)

######*******************************************************
## Graph node degree ~ Trait Z-mean for all traits
xypanels.sexes("dimorphism","degree",full.recast,"sex")
ggsave("figs/Node degree~dimorphism_schiz.jpeg",width=10,height=5)
