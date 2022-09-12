require(pacman)
p_load(tidyverse,skimr)


# Import all data files ---------------------------------------------------
# original 8pop data
# CO=Colorado, CZ=Czech Rep., IS=Israel, NY=New York, RO=Romania, TA=Taiwan, TU=Turkey, UK
pops_8<-read_csv("Data/8pop.MF.TCS.csv") %>% 
        mutate(population=recode(Pop,Colorado="CO")) #some vals are Colorado; should be CO
#just make all the names lowercase as a first step to merging.
names(pops_8)<-tolower(names(pops_8))
#change _ to . for conformity with other pops
names(pops_8) <- gsub("_","\\.",names(pops_8))

#Asian pops (Russia, China, Mongolia, Japan)
pops_asia<-read_csv("Data/all individual phenotype data/asia_pheno_all.csv")
names(pops_asia)<-tolower(names(pops_asia))

#Morocco
pops_mor<-read_csv("Data/all individual phenotype data/morocco_pheno_all.csv")
names(pops_mor)<-tolower(names(pops_mor))

#Egypt
pops_egy<-read_csv("Data/all individual phenotype data/egypt_pheno_all.csv")
names(pops_egy)<-tolower(names(pops_egy))


# Define traits of interest -----------------------------------------------
# And ensure all pops have same names

toi<-c('band','population','year','sex','ci1','rs','tail.mean','mass','date','t.avg.brightness','t.hue','t.chrom','r.avg.brightness','r.hue','r.chrom','b.avg.brightness','b.hue','b.chrom','v.avg.brightness','v.hue','v.chrom')

toi_names_match<-function(df){
  df_names<-names(df)
  sapply(1:length(toi),function(i){
    matched <- toi[i]%in%df_names
    ifelse(!matched,"MISSING","---")
  })
}

#What names are missing/mismatched?
data.frame(TOI=toi,
       not_in_pop8=toi_names_match(df=pops_8),
       not_in_asia=toi_names_match(df=pops_asia),
       not_in_mor=toi_names_match(df=pops_mor),
       not_in_egypt=toi_names_match(df=pops_egy)
       )


# Fix missing variables across all data sets ------------------------------
pops_8 <-
  pops_8 %>% mutate(tail.mean = mean(c(mlts, mrts), na.rm = TRUE)) #add tail.mean
pops_asia <-
  pops_asia %>%
  rename(
    mass = weight,
    t.avg.brightness = t.avg.bright,
    b.avg.brightness = b.avg.bright,
    r.avg.brightness = r.avg.bright,
    v.avg.brightness = v.avg.bright) %>%
  mutate(
    population= pop3, #I believe this is the most relevant column
    date = NA,
    year = NA,
    ci1 = NA,
    rs = NA
  )


# Remove outliers in each pop ---------------------------------------------


