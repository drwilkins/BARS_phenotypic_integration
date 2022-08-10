require(igraph);require(PCIT);require(qgraph);require(psych);require(plotrix)
require(devtools) #necessary to import files from dropbox over https
#Read in my custom function definitions
source_url("https://dl.dropboxusercontent.com/s/b32xvn9x18gahc6/bioworkflow.R") #bioworkflow.R script
D<-read.csv("data/8pop.csv")
head(D)
names(D)[1]<-"Pop"
#check that there are no issues with data entry
check_class(D)
NA_outliers(D,id="ID")
#Get rid of this indiv that has outlier values
D2<-NA_outliers(D,id="ID")$newdata
NA_outliers(D2,id="ID") #outliers switched to NA

#Remove duplicates
D2$popID<-paste(D2$Pop,D2$ID)
sum(duplicated(D2$popID)) #4 duplicates

D2<-subset(D2,!duplicated(D2$popID))
sum(duplicated(D2$popID)) #Removed, now

#Make list w/ each country as a separate dataframe
Dlist<-lapply(levels(D2$Pop),function(x) (subset(D2,Pop==x)))
(names(Dlist)<-levels(D2$Pop))

#Generate trait names so you don't have to type em out
cat(names(D)[c(7,8,11:22)],sep='","')

#define traits of interest
toi<-c("mRWL","mRTS","T_AvgBri","T_Hue","T_Chrom","R_AvgBri","R_Hue","R_Chrom","B_AvgBri","B_Hue","B_Chrom","V_AvgBri","V_Hue","V_Chrom")
toi2<-c("RWL","RTS","TBri","THue","TChr","RBri","RHue","RChr","BBri","BHue","BChr","VBri","VHue","VChr")
shps<-c("square","square",rep("triangle",12))
#Make list of correlations for each country (for traits of interest)
Clist<-lapply(levels(D2$Pop), function(x) cor(as.matrix(Dlist[[x]][,toi]),method="s",use="pairwise.complete"))
(names(Clist)<-levels(D2$Pop))

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
  G<-qgraph(COR,diag=F,fade=F,label.color="black",label.font=2,label.scale=T,label.norm="0000",negCol="black",layout="circle",...)
return(G)}

#Unfiltered networks!
gCO<-Q(Clist$CO)
gCR<-Q(Clist$CR)
gIA<-Q(Clist$IA)
gIL<-Q(Clist$IL)
gRM<-Q(Clist$RM)
gTR<-Q(Clist$TR)
gTW<-Q(Clist$TW)
gUK<-Q(Clist$UK)

##### Filter ze networks!!! (PCIT)
# Get list of PCIT-filtered correlation matrices
Clist_pcit<-lapply(levels(D$Pop),function(x) {
  old<-Clist[[x]] #current dataset is ith element of list
  newdat<-old*0 #initialize 0 matrix w/ same dimensions as old data
  robustnodes<-pcit(old)$idx
  newdat[robustnodes]<-old[robustnodes]#copies only robust nodes over to new dataset (everything else 0)
  return (newdat)})
names(Clist_pcit)<-levels(D2$Pop)

par(mfrow=c(2,4))
for (i in 1: length(Clist_pcit)){
  Q(Clist_pcit[[i]],title=names(Clist_pcit)[i])
}

#############################################
#### Add ARCS to show PCA structure
PCarcs<-function(PCAobj,refPCAobj,...)
{
    if(missing(refPCAobj)){refPCAobj=PCAobj}
  #extract loadings (and order them, as they're sometimes not in numerical order, based on Eigenvalues)
  princ.load<-PCAobj$loadings[,sort(colnames(PCAobj$loadings))] 
  refprinc.load<-refPCAobj$loadings[,sort(colnames(refPCAobj$loadings))] 
  refnodegroups<-max.col(abs(refprinc.load))#Which PC does each trait load maximally on in reference?
  refGroups<-lapply(1:length(colnames(refPCAobj$loadings)),function(x) which(refnodegroups==x)) #List of row indexes for each PC's component traits (in reference)
  
    #object shows how the network should be ordered to match arcs (i.e., 1st element should be first index of trait loading highest on first ref. PC and so on)
  reorder.rows<-data.frame(Row=unlist(refGroups),row.names=rownames(princ.load)[unlist(refGroups)])
  # Reorder princ.load rows according to highest loadings of reference 
  princ.load<-princ.load[unlist(refGroups),]
  
  
  nodegroups<-max.col(abs(princ.load))  # figure out wich axis each node loads on maximally (in focal PCA object)
  Groups<-lapply(1:length(colnames(PCAobj$loadings)),function(x) which(nodegroups==x))#List of row indexes for each PC's component traits (in focal PCA object)
  contig<-lapply(1:length(Groups),function(x) if(length(Groups[[x]])==1){0}else{diff(Groups[[x]])}) #list for all PCs; step size for next arc starting point
  
  ###COMPLICATED for loop and nested if statements for figuring out contiguous segments from a series of numbers (row indexes)
  ## It's designed to figure out arc lengths for mapping PCs onto circular network figures, with nodes ordered by PC loadings 
  # It puts a zero if the trait is not contiguous (next row) to the next trait which loads highest onto the same axis 
  # (e.g. rows 1:4 load on PC1, so the connecting arc length is 3; the next PC1 is trait 6, followed by trait 8, which are disconnected. This block of code would produce 3,0,0 to represent 3 arc lengths)
  Segs=Indx=vector('list',length(Groups))
  for (i in 1: length(Groups))
  {
    sublist<-Groups[[i]]
    segs=exten=vector()
    if(length(sublist)==1){segs<-0;indx<-Groups[[i]][1]}else{
       for (ii in 1:(length(sublist)-1))
       {
              #First run behavior
            if(ii==1)
            {
              if(sublist[ii+1]-sublist[ii]==1){ #If next 2 indexes are 1 apart (contiguous)
                if(ii+1==length(sublist)){ #If end of string, 
                    segs<-1;indx<-Groups[[i]][1]}else{segs=vector();exten=1;indx<-Groups[[i]][1]} #and contiguous, segs=1; if not end of string, segs doesn't exist yet, extension=1
                }else{#If next 2 indexes are more than 1 apart
                   if(ii+1==length(sublist)){ #If end of string,
                    segs<-c(0,0);indx<-c(Groups[[i]][1],Groups[[i]][2])}else{exten<-0;segs=0; indx<-c(Groups[[i]][1],Groups[[i]][2])}#If 1st set nonconsecutive, put zero (put 2 zeros if end of series); 2 indexes should be stored, regardless bc you'll be making 2 arcs
                      }
            }else
              {
                 #Subsequent behavior
                if(sublist[ii+1]-sublist[ii]==1){exten<-exten+1 #If next 2 row indices are 1 apart
                  if(ii+1==length(sublist)){segs<-append(segs,exten)}#If end of string, append exten; no indx statements, bc no beginning of new arc
                }else{ #If next 2 indexes are more than 1 apart
                  if(ii+1==length(sublist)){segs<-append(segs,c(exten,0));indx<-append(indx,Groups[[i]][ii+1])#if end of string, add extra zero; and index for isolated node (i.e. row)
                  }else{segs<-append(segs,exten);exten=0;indx<-append(indx,Groups[[i]][ii+1])} #if not end of string, just put down extension value, and reset exten;put down indx of beginning new arc segment
              
                      }
              }
       
       }
    }
    Segs[[i]]<-segs
    Indx[[i]]<-indx
  }
  
  #Segs is now a list (as long as # of PCs), containing lengths of contiguous arc segments
  #Indx is now a list containing starting points (in terms of row number) for each of the arc segments
  
#   #The Indexes and Segments need to be reversed to accomodate counterclockwise drawing of arcs in draw.arc
#   Indx2<-lapply(Indx, rev)
#   Segs2<-lapply(Segs, rev)
  #Indexes need to be shifted 90 degrees clockwise, as well as being reversed to line up with the qgraph
    
  varnum<-sum(unlist(sapply(Groups,length)))
  Indx2<-lapply(Indx, function(x) sapply(x,function(y) (varnum+1)-y))# This subtracts each index value from the total+1, making them count down instead of up
  Indx3<-Indx2
  for(i in 1:length(Indx2)){for(ii in 1:length(Indx2[[i]])){Indx3[[i]][ii]<-Indx2[[i]][ii]-Segs[[i]][[ii]]  }}

  #Put circle on first for continuity
  draw.circle(0,0,1.25,lwd=1,border='gray')

  arcunits<-(2*pi)/(varnum) #how many radians per node?
  #startarc<-(pi/2)+arcunits#+.05*arcunits #Start 1 node/radians counterclockwise of 12 and draw arcs for PCs, in reverse order, down to PC1
  for (i in 1:length(Segs)){ #draw arc for the # of nodes in radians for each PC
    for (ii in 1:length(Segs[[i]]))
    {
      startarc<-(pi/2)+Indx3[[i]][ii]*arcunits   #Start drawing at pi/2 (12:00)+Node index in arcunits + 1 arcunit (i.e. radian units for 1 node)
      endarc<-startarc+arcunits*(Segs[[i]][ii])
      wiggle<-(arcunits/2.1)#how much to wiggle the arc to make it longer or shorter
      draw.arc(0,0,1.25,angle1=startarc-wiggle,angle2=endarc+wiggle,lwd=6,col=palette()[i+1]) 
    }
  }
  print("Make sure the networks are reordered to match loadings of reference dataset!")
  print("get row indices by doing x<-PCarcs(PCAobj,refPCAobj")
  return(reorder.rows) 
}#End PCarcs
#####################

###Get PCA scores for all dataframes (Dlist)
# Count # factors w/ eigenvalues >=1
extract<-vector(l=8)
for( i in 1:8){
pca<-principal(Dlist[[i]][,toi],rot="none",nfactors=14)
print(names(Dlist)[i]) 
extract[i]<-sum(pca$values>=1)
print(pca)}
extract#number of factors to extract for each population

# Extract PCA objects for all pops, 5PCs, varimax rotation
# Although pops had 4-6 eigenvalues greater than 1; but did 5, as compromise
pcaList<-lapply(1:8,function(x) principal(Dlist[[x]][,toi],rot="varimax",nfactors=extract[x]))
names(pcaList)<-names(Dlist)

#Make another list w/ more specific names for title
Clist_pcit2<-Clist_pcit
names(Clist_pcit2)<-c("Colorado (US)","Czech Rep","Ithaca (US)","Israel","Romania","Turkey","Taiwan","UK")
#Reorder,west to east
Clist_pcit2<-Clist_pcit2[c(1,3,8,2,5,6,4,7)]
names(Clist_pcit2) #Looks better!
#Also make reordered pca list
pcaList2<-pcaList[c(1,3,8,2,5,6,4,7)]

####
# Which correlations are in all networks? And what is the average value?
# * note that Clist_pcit_bin has been "reorganized here"
Clist_pcit_bin<-lapply(1:8,function(x) apply(Clist_pcit2[[x]][reorganize,reorganize],1, function(xx) ifelse(xx==0,0,1)))
names(Clist_pcit_bin)<-names(Clist_pcit2)
# bin8<-lapply(Clist_pcit_bin,function(x) apply(x,c(1,2),sum))
# bin8
commonedges<-teval(paste0("Clist_pcit_bin[[",1:8,"]]",collapse="+"))
diag(commonedges)<-NA
colscale<-seq(1,8,1)
col.indices<-sapply(commonedges,function(x){max(which(x>=colscale))})
pal<-colorRampPalette(c("black","gray")) #palette ranging from (black, unique to 1 pop, to gray, in all pops)
blacks<-pal(8)
COLS<-blacks[col.indices]
COLS<-COLS[which(!is.na(COLS))]

meanNet<-apply(simplify2array(Clist_pcit2), 1:2, mean)
diag(meanNet)<-NA

#test of colors
plot(1:length(COLS),abs(as.vector(commonedges[which(commonedges!=0)])),col=COLS,pch=19,cex=3,xlab="edge index",ylab="# times in a network")


Q(meanNet,edge.color=COLS,title="Mean Phenotype Network",posCol=1,labels=toi2[reorganize],shape=shps[reorganize],vsize=8,vsize2=8)






#run PCarcs function once to create reorganize object
plot(1,1)
reorganize<-PCarcs(pcaList2[[4]],pcaList2[[4]])$Row #ignore error

#now plot the networks w/ arcs, relative to Colorado
pdf("figs/8pop-pcit filtered_gray.pdf",width=16,height=8)
par(mfrow=c(2,4),mar=c(3,3,3,3),xpd=T,oma=rep(1,4),ps=18)
for (i in 1: 8){
  Q(Clist_pcit2[[i]][reorganize,reorganize],title=names(Clist_pcit2)[i],posCol=1,labels=toi2[reorganize],shape=shps[reorganize],vsize=10,vsize2=10,edge.color=COLS)
  PCarcs(pcaList2[[i]],pcaList2[[4]])
  mtext(paste0("n=",pcaList2[[i]]$n.obs),side=1,line=-.75,at=1.1,cex=.9,font=1)
}
dev.off()

pcaList[["CR"]]$loadings[reorganize,] #test it!

#--------------------------------------------
# Make networks using igraph to calculate stats and whatnot

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

#Calculate communities using 4 different algorithms
for(i in 1:8){
  sub<-teval(items[i]) #teval is custom functino in bioworkflow.R
  communities<-list(wt=walktrap.community(sub),fg=fastgreedy.community(sub),
    le=leading.eigenvector.community(sub),eb=edge.betweenness.community(sub),mat=sub)
  tmpname<-paste0(items[i],".comms")#name of object
  assign(tmpname,communities)
}

lay<-layout.fruchterman.reingold(iCO.comms$mat) #set default layout for plots
lay<-layout.circle(iCO.comms$mat)

#function for plotting different communities based on the different algorithms
plotcomm<-function(x,lay,population,...){
  par(mfrow=c(2,2),mar=c(1,1,1,1))
  if(missing(population)){population=""}
  plot(x$wt,x$mat, layout=lay, main="walktrap",edge.width=E(x$mat)$weight*4,vertex.shape=shps,vertex.label=toi2,vertex.label.cex=1.,vertex.label.color="blue",vertex.label.font=2)
  mtext(paste0("Q=",round(modularity(x$wt),4)),1,line=-1.6,at=1.05,cex=.8,font=2,adj=0)
  mtext(paste0("n.groups=",max(x$wt$membership)),1,line=-.9,at=1.05,cex=.8,font=2,adj=0)
  
  mtext(population,3,-.5,at=-1.5,adj=0,font=2,cex=1.05)
  
  plot(x$fg,x$mat, layout=lay, main="fast greedy",edge.width=E(x$mat)$weight*4,vertex.shape=shps,vertex.label=toi2,vertex.label.cex=1.,vertex.label.color="blue",vertex.label.font=2)
   mtext(paste0("Q=",round(modularity(x$fg),4)),1,line=-1.6,at=1.05,cex=.8,font=2,adj=0)
   mtext(paste0("n.groups=",max(x$fg$membership)),1,line=-.9,at=1.05,cex=.8,font=2,adj=0)
  plot(x$le,x$mat, layout=lay, main="leading eigenvector",edge.width=E(x$mat)$weight*4,vertex.shape=shps,vertex.label=toi2,vertex.label.cex=1.,vertex.label.color="blue",vertex.label.font=2)
   mtext(paste0("Q=",round(modularity(x$le),4)),1,line=-1.6,at=1.05,cex=.8,font=2,adj=0)
   mtext(paste0("n.groups=",max(x$le$membership)),1,line=-.9,at=1.05,cex=.8,font=2,adj=0)
  plot(x$eb,x$mat, layout=lay, main="edge betweenness",edge.width=E(x$mat)$weight*4,vertex.shape=shps,vertex.label=toi2,vertex.label.cex=1.,vertex.label.color="blue",vertex.label.font=2)
   mtext(paste0("Q=",round(modularity(x$eb),4)),1,line=-1.6,at=1.05,cex=.8,font=2,adj=0)
   mtext(paste0("n.groups=",max(x$eb$membership)),1,line=-.9,at=1.05,cex=.8,font=2,adj=0)
}

plotcomm(iCO.comms,lay,"Colorado")

par(mar=c(3,6,3,3))
#circular layouts
pdf("figs/8 pop modularity w diff algorithms_circles.pdf",onefile=T,height=8,width=10,colormodel="cmyk")
plotcomm(iCO.comms,lay,"Colorado")
plotcomm(iIA.comms,lay,"Ithaca")
plotcomm(iUK.comms,lay,"UK")
plotcomm(iCR.comms,lay,"Czech Rep")
plotcomm(iRM.comms,lay,"Romania")
plotcomm(iTR.comms,lay,"Turkey")
plotcomm(iIL.comms,lay,"Israel")
plotcomm(iTW.comms,lay,"Taiwan")
dev.off()

#Fruchterman-Reingold layouts
pdf("figs/8 pop modularity w diff algorithms_F-R lay.pdf",onefile=T,height=8,width=10,colormodel="cmyk")
plotcomm(iCO.comms,layout.fruchterman.reingold(iCO.comms$mat),"Colorado")
plotcomm(iIA.comms,layout.fruchterman.reingold(iIA.comms$mat),"Ithaca")
plotcomm(iUK.comms,layout.fruchterman.reingold(iUK.comms$mat),"UK")
plotcomm(iCR.comms,layout.fruchterman.reingold(iCR.comms$mat),"Czech Rep")
plotcomm(iRM.comms,layout.fruchterman.reingold(iRM.comms$mat),"Romania")
plotcomm(iTR.comms,layout.fruchterman.reingold(iTR.comms$mat),"Turkey")
plotcomm(iIL.comms,layout.fruchterman.reingold(iIL.comms$mat),"Israel")
plotcomm(iTW.comms,layout.fruchterman.reingold(iTW.comms$mat),"Taiwan")
dev.off()

