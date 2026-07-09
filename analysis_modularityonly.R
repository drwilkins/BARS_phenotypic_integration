## separate out the modularity analysis to make the script simpler

require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,glue,gtools,PHENIX,dplyr,rsample,pbapply,parallel,lme4,broom, cowplot, RColorBrewer, pheatmap, abind, infotheo)


#import data for 28 populations for which we have at least 12 individuals of each sex
#for Colorado, include only samples from 2008 (before many experiments)
d=read.csv("Data/data_for_submission.csv")
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)

# 
# # Run Subanalysis to account for genetics ---------------------------------
# fst<-read.csv("Data/pairwise_Fst/pairwise_Fst_table.csv")
# 
# pops_w_fst<-fst$population
# #test if any mismatches with population names (next line should be character(0) )
# pops_w_fst[which(is.na(match(pops_w_fst,d$population)))]
# 
# # integ0<-d %>% group_by(population, sex) %>% summarise_at(c("t.chrom","r.chrom", "b.chrom", "v.chrom", "lat"),mean,na.rm=TRUE) %>% 
# #   arrange(sex,population) %>% 
# #   rename(avg_t.chrom=t.chrom,avg_r.chrom=r.chrom, avg_b.chrom=b.chrom, avg_v.chrom=v.chrom,latitude=lat)
# # 
# # integ0$PINT.c <- c(pint.c_females, pint.c_males)
# # 
# # integ_gen=integ0 %>% mutate(use=population %in% rownames(grm)) %>% filter(use==TRUE)
# # 

### Modularity
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)
### Figure 4A: run community detection on each matrix to quantitatively determine good "modules"
#Male data subset by population
data_list_males<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="M")))
names(data_list_males)<-levels(d$population)

#Female data subset by population
data_list_females<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="F")))
names(data_list_females)<-levels(d$population)

traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
traits_col <- traits[-c(1)]

#Male correlations by population
corr_list_males<-lapply(names(data_list_males), function(x) cor(as.matrix(data_list_males[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_males)<-levels(d$population)

#Female correlations by population
corr_list_females<-lapply(names(data_list_females), function(x) cor(as.matrix(data_list_females[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_females)<-levels(d$population)

#filter by signficance
library(Hmisc)
corr_test_males<-lapply(names(data_list_males), function(x) rcorr(as.matrix(data_list_males[[x]][,traits_col]),type="pearson"))
names(corr_test_males)<-levels(d$population)

corr_sigp_males=lapply(corr_test_males, function(x){
  p_mat=x$P
  apply(p_mat, c(1,2), function(x) (x<0.05)+0)
})

corr_test_females<-lapply(names(data_list_females), function(x) rcorr(as.matrix(data_list_females[[x]][,traits_col]),type="pearson"))
names(corr_test_females)<-levels(d$population)

corr_sigp_females=lapply(corr_test_females, function(x){
  p_mat=x$P
  apply(p_mat, c(1,2), function(x) (x<0.05)+0)
})

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



#function to create matrix plots with hierarchical clustering
make_module_mat=function(net.list, threshold=0.2) {
  clusters=lapply(net.list, function(x) {
    g=delete_edges(x, which(E(x)$weight<threshold))
    #g=x
    cluster_fast_greedy(g, weights=E(g)$weight)
  })
  memberships=lapply(clusters, membership)
  comembers=lapply(memberships, function(x) outer(x, x, "==")+0)
  comembers_array=abind(comembers, along=3)
  sum_mat=apply(comembers_array, c(1,2), sum)
  map.data=data.frame(expand.grid(rownames(sum_mat), colnames(sum_mat)), expand.grid(sum_mat))
  names(map.data)=c("Rows", "Columns", "Values")
  map.data$Values=map.data$Values/length(net.list)
  mat=base::as.matrix(pivot_wider(map.data%>%filter(), names_from = Rows, values_from = Values)[,-1])
  rownames(mat)=colnames(mat)
  mat
}
dist_cor <- function(x) as.dist(1-x)
hclust_complete <- function(x) hclust(x, method = "complete")
colors=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)

trait_names=c("Throat Brightness", "Throat Hue", "Throat Chroma", "Breast Brightness", "Breast Hue", "Breast Chroma", "Belly Brightness", "Belly Hue", "Belly Chroma", "Vent Brightness", "Vent Hue", "Vent Chroma")
mat1_male=make_module_mat(nets_male, threshold=0.1)
rownames(mat1_male) = colnames(mat1_male) = trait_names
heatmap(mat1_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,4))

mat2_male=make_module_mat(nets_male, threshold=0.2)
rownames(mat2_male) = colnames(mat2_male) = trait_names
heatmap(mat2_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete)

mat3_male=make_module_mat(nets_male, threshold=0.3)
rownames(mat3_male) = colnames(mat3_male) = trait_names
heatmap(mat3_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete)

pdf("heatmap_male_threshold0.1.pdf")
heatmap(mat1_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

pdf("heatmap_male_threshold0.2.pdf")
heatmap(mat2_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

pdf("heatmap_male_threshold0.3.pdf")
heatmap(mat3_male, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

mat1_female=make_module_mat(nets_female, threshold=0.1)
rownames(mat1_female) = colnames(mat1_female) = trait_names
mat2_female=make_module_mat(nets_female, threshold=0.2)
rownames(mat2_female) = colnames(mat2_female) = trait_names
mat3_female=make_module_mat(nets_female, threshold=0.3)
rownames(mat3_female) = colnames(mat3_female) = trait_names

pdf("heatmap_female_threshold0.1.pdf")
heatmap(mat1_female, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

pdf("heatmap_female_threshold0.2.pdf")
heatmap(mat2_female, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

pdf("heatmap_female_threshold0.3.pdf")
heatmap(mat3_female, scale="none", col=colors, distfun=dist_cor, hclustfun=hclust_complete, margins=c(8,8))
dev.off()

library(vegan)
mantel(mat1_male, mat2_male)
mantel(mat1_male, mat3_male)
mantel(mat2_male, mat3_male)

###make null model using permutation of module assignment for each population, calculate upper 95% confidence threshold to determine module assignments

calc_upper95=function(net.list, threshold=0.2){
upper95_value=vector(length=1000)
for(i in 1:1000){
clusters=lapply(net.list, function(x) {
  g=delete_edges(x, which(E(x)$weight<threshold))
  #g=x
  cluster_fast_greedy(g, weights=E(g)$weight)
})
memberships=lapply(clusters, membership)
memberships_permuted=lapply(memberships, function(x) sample(x, length(x), replace = F))
comembers_permuted=lapply(memberships_permuted, function(x) outer(x, x, "==")+0)
comembers_array_permuted=abind(comembers_permuted, along=3)
sum_mat_permuted=apply(comembers_array_permuted, c(1,2), sum)
map.data_permuted=data.frame(expand.grid(rownames(sum_mat_permuted), colnames(sum_mat_permuted)), expand.grid(sum_mat_permuted))
names(map.data_permuted)=c("Rows", "Columns", "Values")
map.data_permuted$Values=map.data_permuted$Values/length(net.list)
mat_permuted=base::as.matrix(pivot_wider(map.data_permuted%>%filter(), names_from = Rows, values_from = Values)[,-1])
rownames(mat_permuted)=colnames(mat_permuted)
diag(mat_permuted)=NA
upper95_value[i]=quantile(mat_permuted, na.rm=T, probs=c(0.97))
}
upper95_value
}

upper95_male_1=calc_upper95(nets_male, threshold=0.1)
sig_threshold=1-quantile(upper95_male_1, probs=c(0.975))
plot(hclust(as.dist(1-mat1_male)), main="Males, threshold=0.1")
abline(h=sig_threshold, lty=2)

upper95_female_1=calc_upper95(nets_female, threshold=0.1)
sig_threshold=1-quantile(upper95_female_1, probs=c(0.975))
plot(hclust(as.dist(1-mat1_female)), main="Females, threshold=0.1")
abline(h=sig_threshold, lty=2)

upper95_male_2=calc_upper95(nets_male, threshold=0.2)
sig_threshold=1-quantile(upper95_male_2, probs=c(0.975))
plot(hclust(as.dist(1-mat2_male)), main="Males, threshold=0.2")
abline(h=sig_threshold, lty=2)

upper95_female_2=calc_upper95(nets_female, threshold=0.2)
sig_threshold=1-quantile(upper95_female_2, probs=c(0.975))
plot(hclust(as.dist(1-mat2_female)), main="Females, threshold=0.2")
abline(h=sig_threshold, lty=2)

upper95_male_3=calc_upper95(nets_male, threshold=0.3)
sig_threshold=1-quantile(upper95_male_3, probs=c(0.975))
plot(hclust(as.dist(1-mat3_male)), main="Males, threshold=0.3")
abline(h=sig_threshold, lty=2)

upper95_female_3=calc_upper95(nets_female, threshold=0.3)
sig_threshold=1-quantile(upper95_female_3, probs=c(0.975))
plot(hclust(as.dist(1-mat3_female)), main="Females, threshold=0.3")
abline(h=sig_threshold, lty=2)

hc1=hclust(as.dist(1-(mat2_female+mat2_male)/2))
str(hc1)
hc1$order


hc1$labels[hc1$order]
hc1$merge
hc_df=data.frame(hc1$merge, hc1$height)
hc_df
write.csv(hc_df, "hierarchicalclustering_result.csv")

###


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

####
#Make data frame for main figure 
integ0<-d %>% group_by(population, sex) %>% summarise_at(c("t.chrom","r.chrom","t.avg.bright","r.avg.bright", "b.chrom", "v.chrom", "lat"),mean,na.rm=TRUE) %>% 
  arrange(sex,population) %>% 
  rename(mean.t.chrom=t.chrom,mean.r.chrom=r.chrom,mean.t.avg.bright=t.avg.bright,mean.r.chrom=r.chrom,mean.r.avg.bright=r.avg.bright, mean.b.chrom=b.chrom, mean.v.chrom=v.chrom, latitude=lat)

# integ0$network_density <- c(pop_netdensity_females,pop_netdensity_males)
# integ0$pint <- c(pint_females, pint_males)
# integ0 <- integ0 %>% arrange(sex,desc(network_density)) 

mods.dat=bind_rows(list(mods.dat.female, mods.dat.male))
#now combine the population-level color data with modularity data
integ=mods.dat%>% left_join(., integ0) 


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

dat2=dat2 %>% left_join(., fst)
dat2
# ## avg ratio analyses
# summary(lm(avgratio_1~mean.t.chrom, data=integ %>% filter(mean.t.chrom > 0.45, sex=="M")))
# summary(lm(avgratio_1~mean.t.chrom, data=integ %>% filter(mean.t.chrom > 0.45, sex=="F")))
# 
# summary(lm(avgratio_1~mean.r.chrom, data=integ %>% filter(sex=="M")))

### Build Figure 4C,D
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

##
#interaction between within-module and between module
## revision--realizing that we should actually only look at within the focal patch vs. between patch

t_int_m=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter(sex=="M"))
summary(t_int_m)

t_int_f=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter( sex=="F"))
summary(t_int_f)

r_int_m=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter( sex=="M"))
summary(r_int_m)

r_int_f=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter( sex=="F"))
summary(r_int_f)

