##breaking the correlation into 2 modules 

require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,ggiraph,glue,ggnetwork,gtools,colourvalues,PHENIX,dplyr,rsample,pbapply,RColorBrewer)
require(cowplot)
require(colorBlindness)
require(viridis)
#Import all data
d0<-read_csv("Data/all_populations.csv")
nrow(d0)


# Setup -------------------------------------------------------------------
traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
traits_col <- traits[-c(1)]

#limit analysis to pops with at least N samples
(pop_summary<-d0 %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame())

#Let's say 20 is our minimum number of each sex
min_samples<-12
pops_w_min_samples<-pop_summary %>% filter(n_F>=min_samples & n_M>=min_samples)
nrow(pops_w_min_samples) #28 populations with at least 12 individuals
d<-d0 %>% filter(population %in% pops_w_min_samples$population)
nrow(d)
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)

#Let's only work with Colorado samples from 2008 (before many experiments)
d<-d %>% filter(population!="colorado"|population=="colorado"&year==2008)
#Now CO has a more comparable N to other pops
d %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame()

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

##DS Trying to calculate modularity, ala Mel & Marroig
#define a matrix of traits in same modules
# patches=c(rep(1, 3), rep(2, 3), rep(3,3), rep(4,3))
# same.patch=outer(patches, patches, FUN="==")
# same.patch
# patch.names=c("Throat", "Breast", "Belly", "Vent")
# modules=matrix(nrow=length(patches), ncol=length(patches))
# for(i in 1:4){
# modules[which(patches==i), which(patches==i)] = i
# }
# modules

#throat vs. others
patches=c(rep(1, 3), rep(2, 9))
same.patch=outer(patches, patches, FUN="==")
same.patch
patch.names=c("Throat", "Breast-Belly-Vent")
modules=matrix(nrow=length(patches), ncol=length(patches))
for(i in 1:2){
  modules[which(patches==i), which(patches==i)] = i
}
modules
#now take the list of correlation matrices for males across populations and calculate average correlation coefficients within modules (1-4) and between modules (originating from module 1-4)
# mods.list.male=lapply(corr_list_males, function(x) {
#   diag(x)=NA
#   wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
#   wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
#   wi_mod3=mean(abs(x[which(modules==3)]), na.rm=T)
#   wi_mod4=mean(abs(x[which(modules==4)]), na.rm=T)
#   btw_mod=mean(abs(x[which(is.na(modules))]))
#   btw_12=mean(abs(x[1:3, 4:6]))
#   data.frame(wi_mod1, wi_mod2, wi_mod3, wi_mod4, btw_mod, btw_12)
# })

#throat vs. not throat
mods.list.male=lapply(corr_list_males, function(x) {
  diag(x)=NA
  wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
  wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
  wi_mod_both=mean(abs(x[which(is.na(modules)==FALSE)]), na.rm=T)
  btw_mod=mean(abs(x[which(is.na(modules))]))
  data.frame(wi_mod1, wi_mod2, wi_mod_both, btw_mod)
})
#organize results into dataframe
mods.dat.male=tibble(bind_rows(mods.list.male))
mods.dat.male$sex="M"
mods.dat.male$population=names(corr_list_males)

#now do the same for female
# mods.list.female=lapply(corr_list_females, function(x) {
#   diag(x)=NA
#   wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
#   wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
#   wi_mod3=mean(abs(x[which(modules==3)]), na.rm=T)
#   wi_mod4=mean(abs(x[which(modules==4)]), na.rm=T)
#   btw_mod=mean(abs(x[which(is.na(modules))]))
#   btw_12=mean(abs(x[1:3, 4:6]))
#   data.frame(wi_mod1, wi_mod2, wi_mod3, wi_mod4, btw_mod, btw_12)
# })

#throat vs. not throat
mods.list.female=lapply(corr_list_females, function(x) {
  diag(x)=NA
  wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
  wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
  wi_mod_both=mean(abs(x[which(modules==1|modules==2)]), na.rm=T)
  btw_mod=mean(abs(x[which(is.na(modules))]))
  data.frame(wi_mod1, wi_mod2, wi_mod_both, btw_mod)
})

mods.dat.female=tibble(bind_rows(mods.list.female))
mods.dat.female$sex="F"
mods.dat.female$population=names(corr_list_females)

#merge the male and female data by population
mods.dat=mods.dat.female %>% full_join(., mods.dat.male) %>% mutate(avgratio_1=wi_mod1/btw_mod, avgratio_2=wi_mod2/btw_mod, avgratio_both=wi_mod_both/btw_mod)



#Make data frame for main figure (with throat and breast chroma and network density)
integ0<-d %>% group_by(population, sex) %>% summarise_at(c("t.chrom","r.chrom","t.avg.bright","r.avg.bright", "b.chrom", "v.chrom", "lat"),mean,na.rm=TRUE) %>% 
  arrange(sex,population) %>% 
  rename(mean.t.chrom=t.chrom,mean.r.chrom=r.chrom,mean.t.avg.bright=t.avg.bright,mean.r.chrom=r.chrom,mean.r.avg.bright=r.avg.bright, mean.b.chrom=b.chrom, mean.v.chrom=v.chrom, latitude=lat)

integ0$network_density <- c(pop_netdensity_females,pop_netdensity_males)
integ0$pint <- c(pint_females, pint_males)
integ0 <- integ0 %>% arrange(sex,desc(network_density)) 

#now combine the population-level color data with modularity data
integ=mods.dat %>% left_join(., integ0) 
###

##modularity figures

head(integ)

dat2=integ %>% select(wi_mod1, wi_mod2, btw_mod, ends_with("chrom"), sex, population) %>% 
  rename(wi_t=wi_mod1, wi_r=wi_mod2, btw=btw_mod) %>%
  pivot_longer(-c(starts_with("mean"),sex, population), names_to="edge.type", values_to="edge.weight") %>%
  mutate(wi_btw=str_sub(edge.type, start=1, end=2)) %>%
  mutate(wi_btw = replace(wi_btw, wi_btw=="wi", 1)) %>%
  mutate(wi_btw = replace(wi_btw, wi_btw=="bt", 2)) %>%
  mutate(patch=str_sub(edge.type, start=4, end=6)) %>%
  mutate(patch = replace(patch, patch=="t", "within throat")) %>%
  mutate(patch = replace(patch, patch=="r", "within other patches")) %>%
  mutate(patch = replace(patch, patch=="", "between modules")) 

## plot avg ratio
# ggplot(integ %>% filter(mean.t.chrom > 0.45, sex=="M"), aes(x=mean.t.chrom, y=avgratio_1)) +
#   geom_smooth(method="lm")
# 
# summary(lm(avgratio_1~mean.t.chrom, data=integ %>% filter(mean.t.chrom > 0.45, sex=="M")))
# summary(lm(avgratio_1~mean.t.chrom, data=integ %>% filter(mean.t.chrom > 0.45, sex=="F")))
# 
# summary(lm(avgratio_1~mean.r.chrom, data=integ %>% filter(sex=="M")))
###

modplot1m=ggplot(dat2 %>% filter(mean.t.chrom > 0.45, sex=="M"), aes(x=mean.t.chrom, y=edge.weight, fill=patch, color=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  xlim(0.469, 0.581) +
  ylim(0,0.7) +
  theme_cowplot() +
  ylab("Average edge weight") +
  xlab("Average throat chroma of population") +
  ggtitle("Male") +
  guides(linetype=FALSE)

modplot1m_nolegend=modplot1m + theme(legend.position="none")

modplot1f=ggplot(dat2 %>% filter(mean.t.chrom > 0.45, sex=="F"), aes(x=mean.t.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  xlim(0.469, 0.581) +
  ylim(0,0.7) +
  theme_cowplot() +
  ylab("Average edge weight") +
  xlab("Average throat chroma of population") +
  theme(legend.position="none") +
  ggtitle("Female") +
  guides(linetype=FALSE)

#ggsave("modularity_t.r_throat.jpg", width=10, height=4, bg="white")

##unfiltered
# ggplot(dat2, aes(x=mean.t.chrom, y=edge.weight, color=patch)) +
#   geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
#   geom_point()+
#   facet_wrap(~sex) +
#   theme_cowplot() +
#   ylab("Average edge eeight within/across modules") +
#   xlab("Average throat chroma of population") 

modplot2m=ggplot(dat2 %>% filter(sex=="M"), aes(x=mean.r.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average breast chroma of population") +
  theme(legend.position="none") +
  ggtitle("Male") +
  guides(linetype=FALSE) 
  
modplot2f=ggplot(dat2 %>% filter(sex=="F"), aes(x=mean.r.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average breast chroma of population") +
  theme(legend.position="none") +
  ggtitle("Female") +
  guides(linetype=FALSE) 

legend_plot=get_legend(modplot1m)

plot_grid(modplot1m_nolegend, modplot2m, NULL, modplot1f, modplot2f, legend_plot, nrow=2, rel_widths=c(2,2,1))

ggsave("figs/modularityplot.pdf", width=11, height=8.5, units="in")

#regressions
t_wit_m=lm(edge.weight~mean.t.chrom, data=dat2 %>% filter(mean.t.chrom > 0.45, sex=="M", edge.type=="wi_t"))
summary(t_wit_m)

t_wit_f=lm(edge.weight~mean.t.chrom, data=dat2 %>% filter(mean.t.chrom> 0.45, sex=="F", edge.type=="wi_t"))
summary(t_wit_f)

t_btw_m=lm(edge.weight~mean.t.chrom, data=dat2 %>% filter(mean.t.chrom> 0.45, sex=="M", edge.type=="btw"))
summary(t_btw_m)

t_btw_f=lm(edge.weight~mean.t.chrom, data=dat2 %>% filter(mean.t.chrom> 0.45, sex=="F", edge.type=="btw"))
summary(t_btw_f)

#interaction between within-module and between module
t_int_m=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter(mean.t.chrom> 0.45, sex=="M"))
summary(t_int_m)

t_int_f=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter(mean.t.chrom> 0.45, sex=="F"))
summary(t_int_f)

r_int_m=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter(mean.r.chrom> 0.45, sex=="M"))
summary(r_int_m)

r_int_f=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter(mean.r.chrom> 0.45, sex=="F"))
summary(r_int_f)

#ggsave("modularity_t.r_breast.jpg", width=10, height=4, bg="white")


# 
# ggplot(dat2, aes(x=mean.b.chrom, y=edge.weight, color=patch)) +
#   geom_smooth( method="lm", se=F) +
#   geom_point()+
#   facet_wrap(~sex) +
#   theme_classic() +
#   ggtitle("by belly chroma")
# 
# ggplot(dat2, aes(x=mean.v.chrom, y=edge.weight, color=patch)) +
#   geom_smooth( method="lm", se=F) +
#   geom_point()+
#   facet_wrap(~sex) +
#   theme_classic() +
#   ggtitle("by vent chroma")

### plot avg ratio
# avgratios=dat2 %>% group_by(population, sex, patch) %>% 
#   summarise(means=mean(edge.weight)) %>% 
#   pivot_wider(id_cols=c(population, sex), names_from=patch, values_from=means) %>% 
#   mutate(avgratio_1=throat-between, avgratio_2=`breast/belly/vent`-between) %>% 
#   inner_join(., dat2, multiple="first") %>%
#   select(population, sex, avgratio_1, avgratio_2, starts_with("mean")) %>%
#   pivot_longer(cols=c(starts_with("avg")), names_to="type", values_to="ratio")
# 
# ggplot(avgratios%>% filter(mean.t.chrom > 0.45), aes(x=mean.t.chrom, y=ratio, color=type)) +
#   geom_smooth( method="lm", se=F) +
#   geom_point()+
#   facet_wrap(~sex) +
#   theme_classic() +
#   ggtitle("by throat chroma")
# 
# ggplot(avgratios, aes(x=mean.r.chrom, y=ratio, color=type)) +
#   geom_smooth( method="lm", se=F) +
#   geom_point()+
#   facet_wrap(~sex) +
#   theme_classic() +
#   ggtitle("by breast chroma")

###run community detection on each matrix to quantitatively determine good "modules"

nets_male=lapply(corr_list_males, function(x) {
  absmat=abs(x)
  diag(absmat)=0
  graph_from_adjacency_matrix(absmat, "undirected", weighted=T)
})

nets_female=lapply(corr_list_females, function(x) {
  absmat=abs(x)
  diag(absmat)=0
  graph_from_adjacency_matrix(absmat, "undirected", weighted=T)
})

## just shorthand for now, removing lower 20% of correlations. Need to figure out a package to use for filtering now that PCIT is defunct.
clusters_male=lapply(nets_male, function(x) {
g=delete.edges(x, which(E(x)$weight<quantile(E(x)$weight, probs=0.2)))
cluster_fast_greedy(g, weights=E(g)$weight)
})

clusters_female=lapply(nets_female, function(x) {
  g=delete.edges(x, which(E(x)$weight<quantile(E(x)$weight, probs=0.2)))
  cluster_fast_greedy(g, weights=E(g)$weight)
})

memberships_male=lapply(clusters_male, membership)
comembers_male=lapply(memberships_male, function(x) outer(x, x, "==")+0)

memberships_female=lapply(clusters_female, membership)
comembers_female=lapply(memberships_female, function(x) outer(x, x, "==")+0)

library(abind)
comembers_male_array=abind(comembers_male, along=3)
comembers_female_array=abind(comembers_female, along=3)

sum_mat_male=apply(comembers_male_array, c(1,2), sum)
sum_mat_male

sum_mat_female=apply(comembers_female_array, c(1,2), sum)


map.data_male=data.frame(expand.grid(rownames(sum_mat_male), colnames(sum_mat_male)), expand.grid(sum_mat_male))
names(map.data_male)=c("Rows", "Columns", "Values")

matrixplot_male=ggplot(map.data_male, aes(x=Rows, y=Columns, fill=Values)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank()) 

map.data_female=data.frame(expand.grid(rownames(sum_mat_female), colnames(sum_mat_female)), expand.grid(sum_mat_female))
names(map.data_female)=c("Rows", "Columns", "Values")

matrixplot_female=ggplot(map.data_female, aes(x=Rows, y=Columns, fill=Values)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme(legend.position="none", axis.text.x=element_blank(), axis.title=element_blank())

plot_grid(matrixplot_male, matrixplot_female, nrow=1)

#ggsave("figs/modularity_heatmap.pdf", width=11.5, height=5)

net=graph_from_adjacency_matrix(sum_mat_male, "undirected", weighted=T, diag=FALSE)
V(net)$membership=membership(cluster_optimal(net, weights=E(net)$weight))

plot(net, edge.width=E(net)$weight/2, vertex.color=V(net)$membership, vertex.label.dist=2, vertex.label.color="black")

####see if throat comes out more as module when average throat chroma is higher
integ.male=integ %>% filter(sex=="M")
qs.male=quantile(integ.male$mean.t.chrom, probs=c(0.25, 0.5, 0.75))
q1.male=which(integ.male$mean.t.chrom<qs.male[1])
q2.male=which(integ.male$mean.t.chrom<qs.male[2]&integ.male$mean.t.chrom>=qs.male[1])
q3.male=which(integ.male$mean.t.chrom<qs.male[3]&integ.male$mean.t.chrom>=qs.male[2])
q4.male=which(integ.male$mean.t.chrom>=qs.male[3])

get_clusters=function(z,y){
  q.nets=lapply(z[y], function(x) {
    absmat=abs(x)
    diag(absmat)=0
    graph_from_adjacency_matrix(absmat, "undirected", weighted=T)
  })
  
  q.clusters=lapply(q.nets, function(x) {
    g=delete.edges(x, which(E(x)$weight<0.3))
    cluster_fast_greedy(g, weights=E(g)$weight)
  })
  q.clusters
  }

q1.clusters_male=get_clusters(corr_list_males, q1.male)
q2.clusters_male=get_clusters(corr_list_males, q2.male)
q3.clusters_male=get_clusters(corr_list_males, q3.male)
q4.clusters_male=get_clusters(corr_list_males, q4.male)

get_comember_net=function(x){
  require(abind)
  memberships=lapply(x, membership)
  comembers=lapply(memberships, function(x) outer(x, x, "==")+0)
  comembers_array=abind(comembers, along=3)
  sum_mat=apply(comembers_array, c(1,2), sum)
  net=graph_from_adjacency_matrix(sum_mat, "undirected", weighted=T, diag=FALSE)
  V(net)$membership=membership(cluster_optimal(net, weights=E(net)$weight))
  net
}

q1.net=get_comember_net(q1.clusters_male)
q2.net=get_comember_net(q2.clusters_male)
q3.net=get_comember_net(q3.clusters_male)
q4.net=get_comember_net(q4.clusters_male)

net=q4.net
plot(net, edge.width=E(net)$weight/2, vertex.color=V(net)$membership, vertex.label.dist=2, vertex.label.color="black")

par(mfrow=c(2,2))
for(i in 1:4){
  net=list(q1.net, q2.net, q3.net, q4.net)[[i]]
  plot(net, edge.width=E(net)$weight/2, vertex.color=V(net)$membership, vertex.label.dist=2, vertex.label.color="black", main=i)
}

###


summary(lm(wi_mod1~mean.t.chrom, data=integ %>% filter(sex=="F")))
summary(lm(wi_mod2~mean.r.chrom, data=integ %>% filter(sex=="F")))
#summary(lm(wi_mod3~mean.b.chrom, data=integ %>% filter(sex=="F")))
#summary(lm(wi_mod4~mean.v.chrom, data=integ %>% filter(sex=="F")))

summary(lm(avgratio_1~mean.t.chrom, data=integ %>% filter(sex=="F")))
summary(lm(avgratio_2~mean.r.chrom, data=integ %>% filter(sex=="F")))
#summary(lm(avgratio_3~mean.b.chrom, data=integ %>% filter(sex=="F")))
#summary(lm(avgratio_4~mean.v.chrom, data=integ %>% filter(sex=="F")))

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

# ggplot(integ, aes(x=latitude, y=mean.r.chrom)) +
#   geom_point() +
#   facet_wrap(~sex)

# ggplot(integ, aes(x=latitude, y=mean.t.chrom)) +
#   geom_point()+
#   facet_wrap(~sex)

# ggplot(integ, aes(x=latitude, y=pint)) +
#    geom_point()+
#    facet_wrap(~sex)

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
#male_pops<-c("baotu", "morocco", "taiwan")
#female_pops<-c("baotu", "morocco", "taiwan")

male_pops<-c("baotu", "czech rep", "morocco")
female_pops<-c("baotu", "czech rep", "morocco")


#just a check cuz the spelling is all over the place on these names
#if(sum(is.na(poi_ordered))>0){warning("Name mismatch. One of your pops_of_interest not matched to 'integ' df.")}


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
# rawmeansM<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="M") %>% summarise_at(traits_col,mean,na.rm=T)
# 
# rawmeansF<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)
# 

rawmeansM<-d %>% group_by(population) %>% filter(population %in% male_pops, sex=="M") %>% summarise_at(traits_col,mean,na.rm=T) 

rawmeansF<-d %>% group_by(population) %>% filter(population %in% female_pops, sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)
# Function Definitions ----------------------------------------------------
####>>>>>>>>>>>>>>>>>>>>>
## Make custom plot function
Q<-function(COR,lab.col="black",lab.scale=T,lab.font=2,lay="spring",...){
  G<-qgraph(COR,diag=F,fade=F,label.color=lab.col,label.font=lab.font,label.scale=lab.scale,label.norm="0000",mar=c(4,7,7,4),...)
return(G)}

net_layout=layout=matrix(c(1,3,
                           0,2,
                           1,2,
                           
                           2,3,
                           3,2,
                           2,2,
                           
                           2,0,
                           3,1,
                           2,1,
                           
                           1,0,
                           0,1,
                           1,1
                           ), byrow = T, ncol=2)


#<<<<<<<<<<<<<<<
#

traits_col
t.lab=c("TBri", "THue", "TChr", "RBri", "RHue", "RChr", "BBri", "BHue", "BChr", "VBri", "VHue", "VChr")
shps=c("triangle", "triangle", "triangle", "circle", "circle", "circle", "square", "square", "square", "diamond", "diamond", "diamond")

pops=unique(d$population)
### Generate male networks figure
#png("figs/Fig 2. Male_10_Networks_ordered.png",width=13,height=6,units="in",res=300)
pdf("figs/NewFig 2. Male_Networks_modules_all.pdf",width=10,height=14)
par(mfrow=c(7,4),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
  scalarM<-sapply(names(rawmeansM)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansM[,x]),q=50 ))) 
  #make 50 quantiles for matching color scores
  rownames(scalarM)<-rawmeansM$population
  scalarM[,c(1:2,4:5,7:8,10:11)] <-51- scalarM[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
  #define color ramp with 50 gradations
  nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 
  nodepal=brewer.pal(5,"YlGnBu")[c(1,5,2,4,3)]

for (i in 1: length(pops)){
  cur_pop<-pops[i]
  mat<-get_pop_cormat(cur_pop,"M",traits_col)
  g=graph_from_adjacency_matrix(abs(mat), diag=FALSE, weighted=T, mode="undirected")
  nodecolor=nodepal[membership(cluster_fast_greedy(g))]
  #nodecolor<-nodepal[scalarM[as.character(cur_pop),]]
 # groupings<-list(throat=1:3,breast=4:6,belly=7:9,vent=10:12)
  Q(abs(mat),color=nodecolor,border.color="gray20",labels=t.lab,shape=shps,posCol="#181923",negCol=1,vsize=15,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout=net_layout,rescale=TRUE, maximum=1)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)

    #Add bounding rectangle for Egypt
  # if(cur_pop=="Egypt"){
  #   box(which="figure",lwd=3)
    #rect(xleft = -1.6,ybottom = -1.25,xright = 1.25,ytop = 1.6,border="cyan",lwd=3)
  #}
}
dev.off()

################
### Generate female networks figure
#png("figs/Fig 3. Female_10_Networks_ordered.png",width=13,height=6,units="in",res=300)
#png("figs/NewFig 1c. Female_Networks_modules_circle.png",width=13,height=5,units="in",res=300)
pdf("figs/NewFig 2. Female_Networks_modules.pdf",width=13,height=5)
par(mfrow=c(1,3),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
  scalarF<-sapply(names(rawmeansF)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansF[,x]),q=50 ))) 
  #make 50 quantiles for matching color scores
  rownames(scalarF)<-rawmeansF$population
  scalarF[,c(1:2,4:5,7:8,10:11)] <-51- scalarF[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
  #define color ramp with 50 gradations
  nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 

for (i in 1: length(female_pops)){
  cur_pop<-female_pops[i]
  mat<-get_pop_cormat(cur_pop,"F",traits_col)
  nodecolor<-nodepal[scalarF[as.character(cur_pop),]]

  Q(mat,color=nodecolor,border.color="gray20",labels=t.lab,shape=shps,posCol="#181923",negCol=1,vsize=15,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout=net_layout,rescale=TRUE, maximum=1)
  
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
selection<-function(df,pop,year=NULL, which_sex,rawtrait,rawfitmetric){ 
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

ds<-d %>% filter(complete.cases(ci1))
unique(ds$population)

selection(d,"ithaca",2009,c("M"),"r.chrom", "rs")
selection(d,"colorado",2009,c("M"),"r.chrom", "ci1")

selection(d,"colorado",2013,c("M"),"t.chrom", "rs")
selection(d,"colorado",2013,c("M"),"t.chrom", "ci1")


# Looking at Dimorphism ---------------------------------------------------
# Code from Wilkins et al. "Analysis of female song provides insight..." AnBeh 2020

##### Bootstrap correlation estimates to get confidence intervals
dimorphCal <- function(split, columns) {
  require(tidyr)
  require(dplyr)
  #split= split list from rsample::samples(); 
  #columns= columns you want to calculate differences from
  
  if (is.data.frame(split)) {
    d <- split
  } else{
    d <- analysis(split)
  }
  if (length(columns) == 1) {
    sd.f = sd(unlist(d[which(d$sex == "F" |
                               d$sex == "f"), columns]), na.rm =
                T)
    sd.m = sd(unlist(d[which(d$sex == "M" |
                               d$sex == "m"), columns]), na.rm = T)
    pooledSD <- sqrt((sd.f ^ 2 + sd.m ^ 2) / 2)
  } else{
     pooledSD <-apply(d[, columns], 2, function(x)
    {
      sd.f = sd(x[which(d$sex == "F" | d$sex == "f")], na.rm = T)
      sd.m = sd(x[which(d$sex == "M" | d$sex == "m")], na.rm = T)
      sqrt((sd.f ^ 2 + sd.m ^ 2) / 2)
    })
  }
  
  dfmeans <-
    as_tibble(d[, c("sex", columns)]) %>%  group_by(sex) %>% summarize_if(is.numeric, mean, na.rm =
                                                                            T)
  (dfmeans[which(tolower(dfmeans$sex) == "f"),-1] - dfmeans[which(tolower(dfmeans$sex) == "m"),-1]) / pooledSD
  
}

#SETUP
set.seed(99)
#populations to iterate over
d

populations<-unique(d$population)
#bootstraps
bootn<-1000
#Bootstrap dichromatism
#Time Consuming
message("Bootstrapping populations ",bootn," times & calculating dichromatism across populations.")
bootDC0<-pbapply::pblapply(1:length(populations),function(i){
                      #cols to carry to output
                      keep=c("population","country","lat","long")
                      
                      d_i<-subset(d,population==populations[i])
                      years<-paste(unique(d_i$year),collapse=", ")
                      boots <- bootstraps(d_i,bootn,strata=c("sex"))
                      boot_diffs <- lapply(boots$splits,dimorphCal,columns=traits_col) %>% 
                        dplyr::bind_rows()
                      boot_means<-boot_diffs %>% 
                        dplyr::summarise_all(mean,na.rm=TRUE)%>% 
                        dplyr::rename_with(~paste0("DC_",.))
                      out<-dplyr::tibble(d_i[1, keep], 
                                         years = years, 
                                         n_boot = nrow(boots), 
                                         m_samples_per_boot=min_samples,
                                         boot_means)
}) %>% dplyr::bind_rows() 
#Bootstrap trait dimorphism (F-M)/SDpooled
bootDC0

pint_info<-tidyr::pivot_wider(integ %>% dplyr::select(population, sex,pint,mean.r.chrom) ,
                   names_from= sex,
                   names_glue="{sex}_{.value}",
                   values_from= c(pint,mean.r.chrom),
                   id_cols = population)
#Add in PI for M & F
bootDC<-bootDC0 %>% left_join(.,pint_info)

#Want to plot Dichromatism against Phenotypic Integration in M & F
# Expect high dichromatism to correlate with phenotypic integration
# 
# PLOT Phenotypic Integration against dichromatism for all populations
g_m<-bootDC %>% ggplot(aes(x=DC_r.chrom,y=M_pint)) +
  geom_point()+ 
  stat_ellipse()+
  ggrepel::geom_label_repel(aes(label=population))+
  labs(x=expression(atop(bold(Dichromatism~"in"~Breast~Chroma),"<--Darker Males")),
       y=expression(bold("Phenotypic Integration of All Ventral Color Traits")),title="Males")

g_f <- bootDC %>% ggplot(aes(x=DC_r.chrom,y=F_pint)) +
  geom_point()+ 
  stat_ellipse()+
  ggrepel::geom_label_repel(aes(label=population))+
    labs(x=expression(atop(bold(Dichromatism~"in"~Breast~Chroma),"<--Darker Males")),
       y="",title="Females")

#################################
#Output Combined Plot
g_m+g_f  
ggsave("figs/Fig 4. Phenotypic Integration ~ Dichromatism in Breast Chroma.png")
# This figure shows average Pint for both sexes (not bootstrapped currently...) plotted against Dichromatism of Breast Chroma. 
#################################


#Testing whether Lower breast dichromatism (>breast chroma in males) correlates with higher phenotypic integration...
cor.test(bootDC$DC_r.chrom,bootDC$M_pint,method = "s") #nonsignificant nonparametric test
cor.test(bootDC$DC_r.chrom,bootDC$M_pint,method = "p") #significant parametric correlation
#Very nonsignificant for females (both types of test)
cor.test(bootDC$DC_r.chrom,bootDC$F_pint,method = "s")
cor.test(bootDC$DC_r.chrom,bootDC$F_pint,method = "p")

# #What about looking at average dichromatism across all patches?
#  --this is flawed b/c directionality of Avg.Brightness vs Chroma should be accounted for
#  
# bootDC2 <- bootDC %>% dplyr::group_by(population) %>%  
#   dplyr::mutate(mean_DC=rowMeans(dplyr::across(dplyr::starts_with("DC_"))))
# 
# bootDC2 %>% ggplot(aes(x=mean_DC,y=M_pint)) +
#   geom_point()+ 
#   stat_ellipse()+
#   ggrepel::geom_label_repel(aes(label=population))

###map
library(ggmap)
library(maps)
pop_locations=d %>% group_by(population) %>% 
  summarize(lat=mean(lat), long=mean(long)) %>%
  drop_na() %>%
  mutate(is_moscow=population=="moscow")

world_coordinates <- map_data("world") 
ggplot()+
  geom_map( data = world_coordinates, map = world_coordinates, aes(long, lat, map_id = region), color = "black", fill = "white", size = 0.2) +
  geom_point(aes(x=long, y=lat, color=is_moscow), data=pop_locations)  +
  scale_color_manual(values=c("black", "red"))
