require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,ggiraph,glue,ggnetwork,gtools,colourvalues,PHENIX,dplyr)
# install_github("nathanhaigh/pcit@v1.6.0")#not on CRAN at the moment (also seems to be broken)

#Import all data
d0<-read_csv("Data/all_populations.csv")
nrow(d0)


# Setup -------------------------------------------------------------------
traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
traits_col <- traits[-c(1)]

#limit analysis to pops with at least N samples
(pop_summary<-d0 %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame())

#Let's say 20 is our minimium number of each sex
min_samples<-12
pops_w_min_samples<-pop_summary %>% filter(n_F>=min_samples & n_M>=min_samples)
nrow(pops_w_min_samples) #29 populations with at least 12 individuals
d<-d0 %>% filter(population %in% pops_w_min_samples$population)
nrow(d)
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)


#####################
#####################
# 1. Test `Network Density ~ Mean Darkness` across pops 
# Calculate population-level stats ----------------------------------------
#Male data subset by population
data_list_males<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="M")))
names(data_list_males)<-levels(d$population)

#Female data subset by population
data_list_females<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="F")))
names(data_list_females)<-levels(d$population)

#Male correlations by population
corr_list_males<-lapply(names(data_list_males), function(x) cor(as.matrix(data_list_males[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_males)<-levels(d$population)

#Female correlations by population
corr_list_females<-lapply(names(data_list_females), function(x) cor(as.matrix(data_list_females[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_females)<-levels(d$population)


pop_netdensity_females <- sapply(corr_list_females,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))

pop_netdensity_males <- sapply(corr_list_males,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))

pint_list_females=lapply(data_list_females, function(x){
  tb=na.omit(x[,traits_col])
  PINT=pint(tb)
})

pint_list_males=lapply(data_list_males, function(x){
  tb=na.omit(x[,traits_col])
  PINT=pint(tb)
})

pint_females=sapply(pint_list_females, function(x) x[[1]])
pint_males=sapply(pint_list_males, function(x) x[[1]])

#Make data frame for main figure (with throat and breast chroma and network density)
integ0<-d %>% group_by(population, sex) %>% summarise_at(c("t.chrom","r.chrom","t.avg.bright","r.avg.bright", "lat"),mean,na.rm=TRUE) %>% 
  arrange(sex,population) %>% 
  rename(mean.t.chrom=t.chrom,mean.r.chrom=r.chrom,mean.t.avg.bright=t.avg.bright,mean.r.chrom=r.chrom,mean.r.avg.bright=r.avg.bright, latitude=lat)

integ0$network_density <- c(pop_netdensity_females,pop_netdensity_males)
integ0$pint <- c(pint_females, pint_males)
integ <- integ0 %>% arrange(sex,desc(network_density))

#throat patch graph
G_t<-ggplot(integ,
       aes(x = mean.t.chrom, y = network_density, fill = mean.t.chrom)) + 
  stat_ellipse() +
  geom_point(size=3,pch=21,col="black") +
  scale_fill_gradient(
    limits = range(integ$mean.t.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
  ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
  xlab("Throat | Average Population Darkness (Chroma)")+
  ylab("Network Density")
#nonsignificant relationship with THROAT darkness & network density for both sexes
cor.test(subset(integ,sex=="F")$mean.t.chrom,
         subset(integ,sex=="F")$network_density,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.t.chrom,
         subset(integ,sex=="M")$network_density,method = "spearman")
                  
                

#breast patch graph
(G_r<-ggplot(integ,
       aes(x = mean.r.chrom, y = network_density, fill = mean.r.chrom)) + 
  stat_ellipse() +
  geom_point(size=3,pch=21,col="black") +
  scale_fill_gradient(
    limits = range(integ$mean.r.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
  ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
  xlab("Breast | Average Population Darkness (Chroma)")+
  ylab("Network Density")
  )

#throat patch graph with PINT (Wagner 1984 method for phenotypic integration)
(G_t_pint<-ggplot(integ,
            aes(x = mean.t.chrom, y = pint, fill = mean.t.chrom)) + 
  stat_ellipse() +
  geom_point(size=3,pch=21,col="black") +
  scale_fill_gradient(
    limits = range(integ$mean.t.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
  ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
  xlab("Throat | Average Population Darkness (Chroma)")+
  ylab("Phenotypic Integration (PINT)")
)


#breast patch graph with PINT
(G_r_pint<-ggplot(integ,
             aes(x = mean.r.chrom, y = pint, fill = mean.r.chrom)) + 
    stat_ellipse() +
    geom_point(size=3,pch=21,col="black") +
    scale_fill_gradient(
      limits = range(integ$mean.r.chrom),
      low = "#FFFFCC",
      high = "#CC6600",
      guide = "none"
    ) + 
    facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
    ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
    xlab("Breast | Average Population Darkness (Chroma)")+
    ylab("Phenotypic Integration (PINT)")
)

#nonsignificant relationship with THROAT darkness & network density for both sexes
cor.test(subset(integ,sex=="F")$mean.t.chrom,
         subset(integ,sex=="F")$network_density,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.t.chrom,
         subset(integ,sex=="M")$network_density,method = "spearman")


#Significant relationship with BREAST darkness & network density for both sexes
cor.test(subset(integ,sex=="F")$mean.r.chrom,
         subset(integ,sex=="F")$network_density,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.r.chrom,
         subset(integ,sex=="M")$network_density,method = "spearman")

# same with PINT -- stronger correlations
cor.test(subset(integ,sex=="F")$mean.r.chrom,
         subset(integ,sex=="F")$pint,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.r.chrom,
         subset(integ,sex=="M")$pint,method = "spearman")

## Quick Figures: Latitude by breast

#ggplot(integ, aes(x=latitude, y=mean.r.chrom)) +
  geom_point() +
  facet_wrap(~sex)

#ggplot(integ, aes(x=latitude, y=mean.t.chrom)) +
  geom_point()+
  facet_wrap(~sex)

#ggplot(integ, aes(x=latitude, y=pint)) +
  geom_point()+
  facet_wrap(~sex)

cor.test(subset(integ,sex=="F")$latitude,
         subset(integ,sex=="F")$pint,method = "spearman")

cor.test(subset(integ,sex=="M")$latitude,
         subset(integ,sex=="M")$pint,method = "spearman")
####

# Output Fig 1.  Darker birds have denser color networks (for R, but not T) --------

#patchwork syntax
(G_combined<-G_t/G_r)
ggsave("figs/Fig 1. network density ~ breast + throat chroma.png",dpi=300)

(G_combined<-G_t_pint/G_r_pint)
ggsave("figs/Fig 1.alternate PINT ~ breast + throat chroma.png",dpi=300,width=13,height=10,units="in")


#Pretty interesting that Egypt has such a low network density for its darkness. 



# Make SuppMat 1 (Phen Integ ~ Avg. Bright instead of Chrom) --------------
#throat patch graph with PINT (Wagner 1984 method for phenotypic integration)
(G_t_pint2<-ggplot(integ,
            aes(x = mean.t.avg.bright, y = pint, fill = mean.t.avg.bright)) + 
  stat_ellipse() +
  geom_point(size=3,pch=21,col="black") +
  scale_fill_gradient(
    limits = range(integ$mean.t.avg.bright),
    low = "#CC6600",
    high = "#FFFFCC",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
  ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
  xlab("Throat | Average Population Darkness (Avg. Brightness)")+
  ylab("Phenotypic Integration (PINT)")
)


#breast patch graph with PINT
(G_r_pint2<-ggplot(integ,
             aes(x = mean.r.avg.bright, y = pint, fill = mean.r.avg.bright)) + 
    stat_ellipse() +
    geom_point(size=3,pch=21,col="black") +
    scale_fill_gradient(
      limits = range(integ$mean.r.avg.bright),
      low = "#CC6600",
      high = "#FFFFCC",
      guide = "none"
    ) + 
    facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
    ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
    xlab("Breast | Average Population Darkness (Avg. Brightness)")+
    ylab("Phenotypic Integration (PINT)")
)

#nonsignificant relationship with THROAT darkness & network pint for both sexes
cor.test(subset(integ,sex=="F")$mean.t.avg.bright,
         subset(integ,sex=="F")$pint,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.t.avg.bright,
         subset(integ,sex=="M")$pint,method = "spearman")


#Significant relationship with BREAST darkness & network density for both sexes
cor.test(subset(integ,sex=="F")$mean.r.avg.bright,
         subset(integ,sex=="F")$network_density,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.r.avg.bright,
         subset(integ,sex=="M")$network_density,method = "spearman")

# same with PINT -- stronger correlations
cor.test(subset(integ,sex=="F")$mean.r.avg.bright,
         subset(integ,sex=="F")$pint,method = "spearman")

cor.test(subset(integ,sex=="M")$mean.r.avg.bright,
         subset(integ,sex=="M")$pint,method = "spearman")


# Output Fig S1.  Darker birds have denser color networks (for R, but not T) --------
(G_combined2<-G_t_pint2/G_r_pint2)
ggsave("figs/Fig S1. PINT ~ breast + throat avg brightness.png",dpi=300,width=13,height=10,units="in")





#make interactive version
G_r_inxn<-ggplot(integ,
       aes(x = mean.r.chrom, y = network_density, fill = mean.r.chrom)) + 
  stat_ellipse() +
  geom_point_interactive(size=3,pch=21,col="black",aes(tooltip=glue("Population: {population}\nBreast Chroma: {round(mean.r.chrom,2)}\nNetwork Density: {round(network_density,2)}"),data_id=population)) +
  scale_fill_gradient(
    limits = range(integ$mean.r.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + facet_wrap(~sex,labeller =as_labeller(c(M="Males",F="Females") ),ncol = 1)+
  ggrepel::geom_label_repel(aes(label =population),col="black",max.overlaps = 20,size=2)+
  xlab("Population Mean Breast Darkness (Chroma)")+
  ylab("Network Density")

htmlwidgets::saveWidget(ggiraph::girafe(ggobj=G_r_inxn),file = "figs/interactive_Fig1_network density~ throat chroma.html",selfcontained = TRUE)


#####################
#####################
# 2. Visualize phenotype networks for sampling of populations

#which populations have > 30 indivs sampled for both M & F?
min_samples_big<-15
pop_summary %>% filter(n_F>=min_samples_big & n_M>=min_samples_big)

#Define populations you want to include in phenonet figures
pops_of_interest<-c("CO","TA","TU","IS","UK","Morocco","Egypt","yekaterinburg","zakaltoose","zhangye")

#Make a handy function to order the vector of populations by network density
order_pops_by_net_density<-function(integ_df,pops,which_sex){
  new_df<-integ_df %>% filter(sex==which_sex,population%in% pops) %>% 
    arrange(network_density)
  new_df$population
}
male_pops<-order_pops_by_net_density(integ,pops_of_interest,"M")
female_pops<-order_pops_by_net_density(integ,pops_of_interest,"F")

#just a check cuz the spelling is all over the place on these names
if(sum(is.na(poi_ordered))>0){warning("Name mismatch. One of your pops_of_interest not matched to 'integ' df.")}


#Define f(x) for subsetting data & getting filtered correlation matrix
get_pop_cormat <- function(pop,which_sex,traits){
   d_cor<- d %>% 
  filter(population==pop & sex==which_sex) %>% 
  select(traits_col) %>% 
   cor(.,use="pairwise.complete",method = "spear")
  d_cor[diag(d_cor)]<-NA
  
  #Filter algorithm
  # Here, simply ≥|0.3|
  d_cor_bad_indx<-which(abs(d_cor)<0.3)
  d_cor[d_cor_bad_indx]<-0
  
  d_cor
}

# output Fig 2 & 3. -----------------------------------------------------------
###Setup
#Get means for traits in each population for each sex
rawmeansM<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="M") %>% summarise_at(traits_col,mean,na.rm=T)

rawmeansF<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)

# Function Definitions ----------------------------------------------------
####>>>>>>>>>>>>>>>>>>>>>
## Make custom plot function
Q<-function(COR,lab.col,lab.scale,lab.font,lay,...){
  if(missing(lab.col)){lab.col="black"}
  if(missing(lab.scale)){lab.scale=T}
  if(missing(lab.font)){lab.font=2}
  if(missing(lay)){lay="spring"}
  G<-qgraph(COR,diag=F,fade=F,label.color=lab.col,label.font=lab.font,label.scale=lab.scale,label.norm="0000",mar=c(4,7,7,4),...)
return(G)}
#<<<<<<<<<<<<<<<
#

### Generate male networks figure
png("figs/Fig 2. Male_10_Networks_ordered.png",width=13,height=6,units="in",res=300)
par(mfrow=c(2,5),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
  scalar<-sapply(names(rawmeansM)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansM[,x]),q=50 ))) 
  #make 50 quantiles for matching color scores
  rownames(scalar)<-rawmeansM$population
  scalar[,c(1:2,4:5,7:8,10:11)] <-51- scalar[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
  #define color ramp with 50 gradations
  nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 

for (i in 1: length(male_pops)){
  cur_pop<-male_pops[i]
  mat<-get_pop_cormat(cur_pop,"M",traits_col)
  nodecolor<-nodepal[scalar[as.character(cur_pop),]]
 # groupings<-list(throat=1:3,breast=4:6,belly=7:9,vent=10:12)
  Q(mat,color=nodecolor,border.color="gray20",labels=toi3,shape=shps,posCol="#181923",negCol=1,vsize=20,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout="circle",rescale=TRUE)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)

    #Add bounding rectangle for Egypt
  if(cur_pop=="Egypt"){
    box(which="figure",lwd=3)
    #rect(xleft = -1.6,ybottom = -1.25,xright = 1.25,ytop = 1.6,border="cyan",lwd=3)
  }
  

  
}
dev.off()

################
### Generate female networks figure
png("figs/Fig 3. Female_10_Networks_ordered.png",width=13,height=6,units="in",res=300)
par(mfrow=c(2,5),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
  scalar<-sapply(names(rawmeansF)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansF[,x]),q=50 ))) 
  #make 50 quantiles for matching color scores
  rownames(scalar)<-rawmeansF$population
  scalar[,c(1:2,4:5,7:8,10:11)] <-51- scalar[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
  #define color ramp with 50 gradations
  nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 

for (i in 1: length(female_pops)){
  cur_pop<-female_pops[i]
  print(i)
  mat<-get_pop_cormat(cur_pop,"F",traits_col)
  nodecolor<-nodepal[scalar[as.character(cur_pop),]]

  Q(mat,color=nodecolor,border.color="gray20",labels=toi3,shape=shps,posCol="#181923",negCol=1,vsize=20,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,lay="circle",rescale=TRUE)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)

    #Add bounding rectangle for Egypt
  if(cur_pop=="Egypt"){
    box(which="figure",lwd=3)
    #rect(xleft = -1.6,ybottom = -1.25,xright = 1.25,ytop = 1.6,border="cyan",lwd=3)
  }
  

  
}
dev.off()


# Looking at selection  ---------------------------------------------------

#define function to scale,plot, and test linear relationship for a trait
selection<-function(df,pop,year, which_sex,rawtrait,rawfitmetric){ 
  require(ggplot2)
  #remove incomplete rows
  dataset<-df %>% dplyr::filter(population==pop,year==year,sex %in% which_sex,complete.cases(rawtrait),complete.cases(rawfitmetric))
  trait.Z<-dataset[,rawtrait] %>% scale() %>% as.vector()
  fit<-dataset[,rawfitmetric]%>% unlist()
  relfit<-fit/mean(fit,na.rm=T) 
  
  newdata<-tibble(X=trait.Z,Y=relfit)
  model<-lm(relfit~trait.Z+I(trait.Z^2)) #Is this right? Should differential include squared term?
  summ<-summary(model)

  fitplot <-
    ggplot(newdata, aes(x = X, y = Y))+
    geom_smooth(method = "loess", col ="black") + 
    geom_point() + 
    theme_bw() + 
    ggtitle(pop) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 18, face = "plain"),
      axis.title = element_text(size = 18)
    ) + xlab(rawtrait) +
    ylab("Rel CI")
  fitplot
  s=coef(summ)[2,1]
  se=coef(summ)[2,2]
  g=2*coef(summ)[3,1]#These are supposed to be doubled (Stinchcombe Evolution 2008)
  se2=2*coef(summ)[3,2]#These are supposed to be doubled
  t.val=coef(summ)[2,3]
  t.val2=coef(summ)[3,3]
  p.val=coef(summ)[2,4]
  p.val2=coef(summ)[3,4]
  stats<-data.frame(
    population = pop,
    sex=paste(which_sex,collapse="+"),
    Trait = rawtrait,
    FitMetric = rawfitmetric,
    s,
    se,
    t.val,
    p.val,
    g,
    se2,
    t.val2,
    p.val2
  )
  output<-list(fitplot,stats)
  names(output)<-c("fitplot","stats")
  return(output)
}

selection(d,"colorado",2013,c("M"),"r.chrom", "rs")
selection(d,"colorado",2013,c("M"),"r.chrom", "ci1")

selection(d,"colorado",2013,c("M"),"t.chrom", "rs")
selection(d,"colorado",2013,c("M"),"t.chrom", "ci1")
