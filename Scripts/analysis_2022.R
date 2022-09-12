require(pacman)
p_load(tidyverse,qgraph,igraph,devtools)
install_github("nathanhaigh/pcit@v1.6.0")#not on CRAN at the moment (also seems to be broken)

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
pops_w_min_samples_F<-pop_summary %>% filter(n_F>=min_samples & n_M>=min_samples)
nrow(pops_w_min_samples) #22 populations with at least this many individuals
d<-d0 %>% filter(population %in% pops_w_min_samples$population)
nrow(d)
d$population<-as.factor(d$population)


# Function Definitions ----------------------------------------------------
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
#


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

#Feale correlations by population
corr_list_females<-lapply(names(data_list_females), function(x) cor(as.matrix(data_list_females[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_females)<-levels(d$population)


pop_netdensity_females <- sapply(corr_list_females,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))

pop_netdensity_males <- sapply(corr_list_males,function(x) sum(abs(x[upper.tri(x)]))/sum(upper.tri(x)))

pop_darkness_males <- sapply(data_list_males, function(x) {
  df <- x[, c("r.chrom")]#, "r.chrom", "b.chrom", "v.chrom")]
  return(mean(colMeans(df, na.rm = T)))
})

pop_darkness_females <- sapply(data_list_females, function(x) {
  df <- x[, c("r.chrom")]#, "r.chrom", "b.chrom", "v.chrom")]
  return(mean(colMeans(df, na.rm = T)))
})

integ<-tibble(
  population = rep(names(data_list_females), 2),
  pop_darkness = c(pop_darkness_males, pop_darkness_females),
  network_density = c(pop_netdensity_males, pop_netdensity_females),
  sex = c(rep("M", length(pop_darkness_males)), rep("F", length(
    pop_darkness_females
  )))
)
#something wrong with Colorado data
ggplot(integ,aes(x=pop_darkness,y=network_density,col=pop_darkness))+geom_point()+stat_ellipse()+facet_wrap(~sex)+xlim(0.2,.55)

