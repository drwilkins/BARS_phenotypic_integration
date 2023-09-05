## Generate Map of sampling sites

library(tidyverse)
library(ggmap)

d00<-read_csv("Data/all_populations.csv")
nrow(d00)

traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
#Just the color traits
traits_col <- traits[-c(1)]

#limit analysis to pops with at least N samples
(pop_summary<-d00 %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame())

#Let's say 20 is our minimum number of each sex
min_samples<-12
pops_w_min_samples<-pop_summary %>% filter(n_F>=min_samples & n_M>=min_samples)
pops_w_20_samples<-pop_summary %>% filter(n_F>=20 & n_M>=20)
nrow(pops_w_min_samples) #28 populations with at least 12 individuals
d0<-d00 %>% filter(population %in% pops_w_min_samples$population)
nrow(d0)

pop_location=d0 %>% group_by(population) %>% summarise(pop.lat=mean(lat), pop.long=mean(long))


ggplot(pop_location) +borders("world", colour="black", fill="gray90") +
  theme_bw() +
  ylim(0,90) +
  geom_point(aes(x=pop.long, y=pop.lat), pch=21, fill="red", alpha=0.8, size=4) +
  labs(y="Latitude", x="Longitude") +
  geom_label(data=pop_location %>% filter(population=="baotu"|population=="morocco"|population=="taiwan"), aes(label=population, x=pop.long, y=pop.lat), nudge_y=-5) 

ggsave("figs/Map_samplesites.png", width=10, height=4)
