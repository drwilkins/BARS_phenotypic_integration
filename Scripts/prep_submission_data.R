## clean up population names for submission data


#import data for 28 populations for which we have at least 12 individuals of each sex
#for Colorado, include only samples from 2008 (before many experiments)
d=read.csv("Data/data_for_submission.csv")
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)

#make table of sample sizes
pop_info_n=d %>% group_by(population, sex) %>% summarise(n=n()) %>% pivot_wider(id_cols=population, names_from=sex, values_from=n)

pop_info_latlong= d %>% group_by(population) %>% summarise(lat=first(lat), long=first(long))

pop_info_all=pop_info_latlong %>% left_join(pop_info_n)
pop_info_all

pop_lookup=read_excel("Data/location<->population_key.xlsx")

d_new = d %>% left_join(pop_lookup, by=join_by(population==population)) %>% select(band, Municipality, year, sex, t.avg.bright, t.hue, t.chrom, r.avg.bright, r.hue, r.chrom, b.avg.bright, b.hue, b.chrom, v.avg.bright, v.hue, v.chrom, lat, long) %>%
  dplyr::rename(population=Municipality)

d_new
#write.csv(pop_info_all, "data/population_info_for_supplement.population#write.csv(pop_info_all, "data/population_info_for_supplement.csv")

write.csv(d_new, "Data/data_resubmission.csv")

###same with Fst data

fst<-read.csv("Data/pairwise_Fst/pairwise_Fst_table.csv")

fst_new=fst %>% left_join(pop_lookup, by=join_by(population==population)) %>%
  dplyr::select(Municipality, locality2, weighted_Fst) %>%
  dplyr::rename(population=Municipality)

write.csv(fst_new, "Data/pairwise_Fst/pairwise_Fst_correctpopulation.csv")

## same with pairwise distance data
gr_pair=read.csv("Data/pairwise_distance_pop.csv") 
gr_pair_new= gr_pair %>%
  left_join(pop_lookup, by=join_by(pop1==other_pop_name)) %>%
  dplyr::select(Municipality, pop2, mean.gen.dist) %>%
  dplyr::rename(pop1=Municipality) %>%
  left_join(pop_lookup, by=join_by(pop2==other_pop_name)) %>%
  dplyr::select(pop1, Municipality, mean.gen.dist) %>%
  dplyr::rename(pop2=Municipality)

write.csv(gr_pair_new, "Data/pairwise_distance_correctpopulation.csv", row.names=F)
