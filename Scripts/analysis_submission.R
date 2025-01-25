require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,glue,gtools,PHENIX,dplyr,rsample,pbapply,parallel,lme4,broom, cowplot, RColorBrewer)
#remotes::install_github("galacticpolymath/galacticEdTools")
#require(galacticEdTools)


# Function Definitions -----------------------------------------------------

# BOOTSTRAP PROCEDURE
# Does the following:
#   1. Runs analysis over the traits specified by 'columns', with IDs in 'id'
#   2. Iterates this whole procedure over iterate_over column (default="population")
#   3. Bootstraps data 'boot_n' times (default=1e4), balancing sample sizes over 'strata_var' (default="sex")
#       *'strata_levels'= case-insensitive 2 levels of strata, 
#                         ordered Female first (default=c("f","m"))      
#     ** Bootstraps are parallelized, change cores with the 'cor_num' parameter. 
#        default: cor_num=parallel::detectCores()-1
#   4. Calculates averages for traits specified in 'columns' within each bootstrap sample
#   5. Calculates the corrected phenotypic integration (PINT.c) value using PHENIX::pint()
#   6. Calculates dichromatism (or dimorphism) for traits specified in 'columns' and sex defined by 'strata' using dimorphCal()
#   7. Calculates correlations of 'cor_type' (default = "p"; passed to cor()) between:
#       a. (PI) and trait mean for "toi"
#       b. (PI) and dichromatism across trait(s) in 'toi'
#   NOTES:
#   *   'keep_cols' refers to columns you want to keep in your outputs
#   **  'seed' is set to 99 by default, allowing for repeatable bootstrap results.
#   
#   
boot_analy <- function(df = NULL,
                       columns = NULL,
                       toi = "r.chrom",
                       boot_n = 1000,
                       strata_var = "sex",
                       strata_levels = c("f", "m"),
                       id_col = "band",
                       keep_cols=NULL,
                       iterate_over = "population",
                       cor_type = "p",
                       seed = 99,
                       cor_num=parallel::detectCores()-1) {
  
  set.seed(seed) #to make this repeatable
  if (is.null(df)) {
    message("must input your data as the 'df' parameter")
  }
  if (is.null(columns)) {
    message("specify what columns you want to analyze (i.e. traits to measure dichromatism for)")
  }
  #make case-insensitive levels for filtering
  stratum1 <- c(toupper(strata_levels[1]), tolower(strata_levels[1]))
  stratum2 <- c(toupper(strata_levels[2]), tolower(strata_levels[2]))
 
  #And test just as an extra precaution
  sum_levels_match_expected<-sum(unlist(unique(df[,substitute(strata_var)])) %in% strata_levels)
  
  if(sum_levels_match_expected!=2){
    stop("Levels don't match. Make sure strata_levels=c(val1,val2), found in your data frame.")
  }
  
  #Get a vector of unique entries to iterate over
  uniques <-
    unique(df[, iterate_over]) %>% na.omit() %>% unlist() %>% as.vector()
  
  
  # Iterate analysis over all populations/subgroups -------------------------
  all_subgroups <- pbapply::pblapply(1:length(uniques), function(i) {
    subgroup <- uniques[i]
    message("\n\n# Analyzing ",
            iterate_over,
            " ",
            i,
            " of ",
            length(uniques),
            ": '",
            subgroup,
            "'")
    # Get the subgroup (i.e. subset with population==x)
    d_i <-
      df %>% dplyr::filter(eval(as.symbol(iterate_over)) == subgroup)
    
    #Now do bootstrapping with each subgroup_i, stratified by 'strata'
    boots <- rsample::bootstraps(d_i, boot_n, strata = strata_var)
    
    #######
    # ANALYZE each bootstrap, within population i -----------------------------
    message("   Bootstrap Progress")
    boot_results <- pbapply::pblapply(1:length(boots$splits),cl = parallel::detectCores()-1,FUN= function(ii) {
      
      boot_ii <- boots$splits[[ii]] %>%
        rsample::analysis() %>%
        select(all_of(c(id_col,strata_var,iterate_over,columns,keep_cols))) %>%  #we want the value of these vars
        group_by(!!sym(strata_var))#man this escaping w/ !! vs {{ }} distinction is annoying!
      
      # A) Trait Means
      trait_means_ii <- boot_ii  %>% 
        dplyr::summarise(across(where(is.double), mean, na.rm = TRUE)) %>%
        dplyr::rename_with(~ paste0("avg_", .), where(is.double))
      #preserve keep_cols
      if(!is.null(keep_cols)){
        trait_means_ii <- bind_cols(boot_ii[1:2,keep_cols],trait_means_ii)
      }
      
      # B) Sex Differences in Traits
      
      sex_diffs_ii <-
        dimorphCal(boot_ii,
                   columns,
                   strata_var = strata_var,
                   strata_levels = strata_levels,
                   keep_cols=keep_cols)%>% 
        rename_with(~paste0("SD_",.),all_of(columns))
      
      
      # C) Phenotypic Integration across all traits
      PI_ii <- lapply(strata_levels, function(stratum) {
        
        boot_ii_stratum <- boot_ii %>%
          dplyr::filter(eval(as.symbol(strata_var)) == stratum) %>% ungroup %>%  select(all_of(columns))
        
        PI_ii <- PHENIX::pint(na.omit(boot_ii_stratum))
        dplyr::tibble({
          {
            strata_var
          }
        } := stratum,
        N_PINT = PI_ii$N,
        PINT = PI_ii$PINT,
        PINT.c = PI_ii$PINT.c)
      }) %>% dplyr::bind_rows()
      
      
      # D) Aggregate summary of PI and trait avgs
      
      summary_ii <-PI_ii %>% 
        mutate({{iterate_over}} := subgroup, boot_i = ii) %>%
        relocate({{iterate_over}}, boot_i) %>%
        left_join(x = ., y = trait_means_ii, by = {{strata_var}})
      # E) Return summary data
      list(
        summary = summary_ii,
        sex_diffs = as_tibble(sex_diffs_ii) %>% mutate({
          {
            iterate_over
          }
        } := subgroup, boot_i = ii) %>%
          relocate({
            {
              iterate_over
            }
          }, boot_i)
      )
      
    })  # End bootstraps within population_i
    
    #Combine sex diffs and summary table across all
    boot_summ <-
      lapply(boot_results, function(x)
        x$summary) %>% bind_rows()
    boot_sex_diffs <-
      lapply(boot_results, function(x)
        x$sex_diffs) %>% bind_rows()
    
    
    subgroup_output <- list(boot_summ = boot_summ ,
                            boot_sex_diffs = boot_sex_diffs)
    
    subgroup_output
    
  })# End lapply over all populations/subgroups
  
  #Organize final output
  boot_summ_tibble <-
    lapply(1:length(all_subgroups), function(i) {
      all_subgroups[[i]]$boot_summ
    }) %>% bind_rows
  boot_sex_diffs_tibble <-
    lapply(1:length(all_subgroups), function(i) {
      all_subgroups[[i]]$boot_sex_diffs
    }) %>% bind_rows
  
  mean_traits_tibble <- boot_summ_tibble %>% select(-boot_i) %>%
    group_by(!!sym(iterate_over), !!sym(strata_var),!!!syms(keep_cols)) %>% 
    summarize_if(is.double, mean, na.rm =TRUE)
  
  mean_sex_diffs_tibble <-
    boot_sex_diffs_tibble %>% select(-boot_i) %>%
    group_by(!!sym(iterate_over),!!!syms(keep_cols)) %>% 
    summarize(across(all_of(paste0("SD_",columns)), mean, na.rm=TRUE))
  
  #calculate correlations bw phenotypic integration and Trait Means for Trait(s) of Interest
  cor_trait_avgs_df <- lapply(1:length(toi), function(i) {
    lapply(1:boot_n, function(j) {
      #stratum1 is usually going to be female (the first entry in 'strata_levels')
      df_j1 <-
        boot_summ_tibble %>% dplyr::filter(boot_i == j &
                                             eval(as.symbol(strata_var)) %in% stratum1)
      df_j2 <-
        boot_summ_tibble %>% dplyr::filter(boot_i == j &
                                             eval(as.symbol(strata_var)) %in% stratum2)
      #usually the female correlation
      #Note, we're using the PINT.c value, which corrects for number of individuals and traits 
      #for each population
      cor_1 <-
        cor(df_j1$PINT.c ,
            df_j1[, paste0("avg_", toi[i])] %>% unlist(),
            method = cor_type,
            use = "complete.obs")
      #usually the male correlation
      cor_2 <-
        cor(df_j2$PINT.c ,
            df_j2[, paste0("avg_", toi[i])] %>% unlist(),
            method = cor_type,
            use = "complete.obs")
      tibble({
        {
          strata_var
        }
      } := strata_levels, trait_x1 = "PINT.c", trait_x2 = paste0("avg_", toi[i]), cor =
        c(cor_1, cor_2))
    }) %>% bind_rows
  }) %>% bind_rows
  
  #Collapse correlations across bootstrap correlations for trait means
  #Just temporarily disabling this b/c of a persistent message
  options(dplyr.summarise.inform = F)
  cor_trait_avgs_tibble <- cor_trait_avgs_df %>%
    group_by(!!sym(strata_var), trait_x1, trait_x2) %>% summarize(
      mean_cor = mean(cor, na.rm = T),
      ci_low = quantile(
        cor,
        probs = .025,
        na.rm = T,
        names = F,
        type = 7
      ),
      ci_high = quantile(
        cor,
        probs = .975,
        na.rm = T,
        names = F,
        type = 7
      )
    ) %>% arrange(trait_x2)
  
  
  #Calculate correlations between phenotypic integration and sex difference in toi(s)
  cor_sex_diff_df <- lapply(1:length(toi), function(i) {
    lapply(1:boot_n, function(j) {
      #stratum1 is usually going to be female (the first entry in 'strata_levels')
      sex_diff_j <-
        boot_sex_diffs_tibble %>% dplyr::filter(boot_i == j)
      df_j1 <-
        boot_summ_tibble %>% dplyr::filter(boot_i == j &
                                             eval(as.symbol(strata_var)) %in% stratum1)
      df_j2 <-
        boot_summ_tibble %>% dplyr::filter(boot_i == j &
                                             eval(as.symbol(strata_var)) %in% stratum2)
      #usually the female correlation
      cor_1 <- cor(sex_diff_j[, paste0("SD_",toi[i])] %>% unlist() ,
                   df_j1$PINT.c,
                   method = cor_type,
                   use = "complete.obs")
      #usually the male correlation
      cor_2 <- cor(sex_diff_j[,paste0("SD_",toi[i])] %>% unlist() ,
                   df_j2$PINT.c,
                   method = cor_type,
                   use = "complete.obs")
      tibble({
        {
          strata_var
        }
      } := strata_levels, trait_x1 = "PINT.c", trait_x2 = paste0("sex_diff_", toi[i]), cor =
        c(cor_1, cor_2))
    }) %>% bind_rows
  }) %>% bind_rows
  
  #Collapse correlations across bootstrap correlations for trait means
  cor_trait_sex_diff_tibble <- cor_sex_diff_df %>%
    group_by(!!sym(strata_var), trait_x1, trait_x2) %>%
    summarize(
      mean_cor = mean(cor, na.rm = T),
      ci_low = quantile(
        cor,
        probs = .025,
        na.rm = T,
        names = F,
        type = 7
      ),
      ci_high = quantile(
        cor,
        probs = .975,
        na.rm = T,
        names = F,
        type = 7
      )
    ) %>%
    arrange(trait_x2)
  options(dplyr.summarise.inform = T)#turn messaging back on
  
  message(":) ALL DONE!\n")
  SUMMARY <- list(
    boot_n= boot_n,
    boot_summ = boot_summ_tibble,
    boot_sex_diffs = boot_sex_diffs_tibble,
    mean_traits = mean_traits_tibble,
    mean_sex_diffs = mean_sex_diffs_tibble,
    cor_pint_x_avg_trait = cor_trait_avgs_tibble,
    cor_pint_x_sex_diff = cor_trait_sex_diff_tibble
  )
  
  
}
    
  

# 
# 1. Setup -------------------------------------------------------------------
# #Import all data
# d00<-read_csv("Data/all_populations.csv")
# nrow(d00)
# 
# traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
# #Just the color traits
# traits_col <- traits[-c(1)]
# 
# #limit analysis to pops with at least N samples
# (pop_summary<-d00 %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame())
# 
# #Let's say 12 is our minimum number of each sex
# min_samples<-12
# pops_w_min_samples<-pop_summary %>% filter(n_F>=min_samples & n_M>=min_samples)
# pops_w_20_samples<-pop_summary %>% filter(n_F>=20 & n_M>=20)
# nrow(pops_w_min_samples) #28 populations with at least 12 individuals
# d0<-d00 %>% filter(population %in% pops_w_min_samples$population)
# nrow(d0)
# d0$population<-as.factor(d0$population)
# d0$sex<-as.factor(d0$sex)
# 
# #Let's only work with Colorado samples from 2008 (before many experiments)
# d<-d0 %>% filter(population!="colorado"|population=="colorado"&year==2008)
# #Now CO has a more comparable N to other pops
# d %>% group_by(population,sex) %>%  summarise(n=n()) %>%  pivot_wider(names_from=sex,values_from=n,names_prefix = "n_") %>% mutate(n_TOT=n_F+n_M) %>% as.data.frame()
# 
# #Export info on just the populations we're using &rename Taiwan to Taipei, Taiwan for consistencey
# d %>% select(population,location,year,lat,long,hybrid_zone,zone) %>%
#   mutate(location=case_when(location=="Taiwan"~"Taipei, Taiwan",.default=location)) %>%
#   distinct(population,.keep_all = T) %>% write_csv(.,file="Data/populations_analyzed.csv")
#d=d%>%dplyr::select(band, population, year, sex, tidyselect::starts_with("t."), starts_with("r."), starts_with("b."), starts_with("v"), lat, long)

write.csv(d, "Data/data_for_submission.csv")

# RUNNING THE BOOTSTRAP PROCEDURE

#import data for 28 populations for which we have at least 12 individuals of each sex
#for Colorado, include only samples from 2008 (before many experiments)
d=read.csv("Data/data_for_submission.csv")
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)
# VERY Time consuming!
# Uncomment if you want to run the full 10e4 bootstraps

# res <- boot_analy(
#   df = d,
#   columns = traits_col,
#   toi = c("r.avg.bright","r.chrom","b.chrom","b.avg.bright", "t.chrom"),
#   boot_n = 1e4,
#   strata_var="sex",
#   strata_levels = c("F", "M"),
#   keep_cols = c("location", "lat", "long")
# )

#If you have access to our google drive, you can read in the large data file. Not on github.
#results


res<-readRDS("/Users/dshizuka2/Dropbox/Dai_Research/Main Projects/BARS_phenotypicintegration/results_10k_bootstraps.RDS")

#res<-readRDS("/Users/dshizuka2/Dropbox/Dai_Research/manuscripts/BARS_PhenotypicIntegration/results_10k_bootstraps.RDS")

##Dai's link
# res<-readRDS("/Users/daishizuka/Dropbox/Dai_Research/Main Projects/BARS_phenotypicintegration/results_10k_bootstraps.RDS")



# Run Subanalysis to account for genetics ---------------------------------
fst<-read.csv("Data/pairwise_Fst/pairwise_Fst_table.csv")

pops_w_fst<-fst$population
#test if any mismatches with population names (next line should be character(0) )
pops_w_fst[which(is.na(match(pops_w_fst,d$population)))]

# do main analysis as a glmm of PINT.c~avg_r.chrom + fst from Egypt
d_gen<-res$boot_sum %>% filter(population %in% pops_w_fst)
#how many pops in d_gen?
d_gen$population %>% unique() %>% length


#Add genetic info
d_gen<-left_join(d_gen,fst[,c("population","weighted_Fst")],by="population") %>% select(population,location,sex,boot_i,PINT,PINT.c,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,weighted_Fst)
d_gen_pops<-d_gen %>% distinct(population) %>% unlist
d_gen_m<-d_gen %>% filter(sex=="M")
d_gen_f<-d_gen %>% filter(sex=="F")

#run the linear models for each color patch using the bootstrap results and calculate average and CI of estimates. 
# #Run all of the models and save the results so that we don't have to run them each time.
# 
# lm_results_r.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
#   pops_boot_i <- d_gen %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_r.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_r.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
# }) %>% bind_rows
# 
# lm_results_t.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
#   pops_boot_i <- d_gen %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_t.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_t.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
# }) %>% bind_rows
# 
# #Running model with sex as a covariate and belly chroma
# lm_results_b.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
#   pops_boot_i <- d_gen %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_b.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_b.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
# }) %>% bind_rows
# 
# #Running model with sex as a covariate and vent chroma
# lm_results_v.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
#   pops_boot_i <- d_gen %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_v.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_v.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
# }) %>% bind_rows
# 
# save(lm_results_r.chrom_fst, lm_results_t.chrom_fst, lm_results_b.chrom_fst, lm_results_v.chrom_fst, file="Data/PINTanalysis_results_w_Fst.rds")

load("Data/PINTanalysis_results_w_Fst.rds")
# This is the result reported in paper
#95CIs for analysis with sex as covariate, rather than splitting up analyses
#Significant effect of breast darkness on phenotypic integration (PCIT)
quantile(lm_results_r.chrom_fst$est_avg_r.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_r.chrom_fst$est_avg_r.chrom,na.rm=T,names=F)
#Significant effect of population genetics (Fst) on phenotypic integration (PCIT)
quantile(lm_results_r.chrom_fst$est_Fst,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_r.chrom_fst$est_Fst,na.rm=T,names=F)
#Nonsignificant effect of sex on phenotypic integration (PCIT)
quantile(lm_results_r.chrom_fst$est_sexM,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_r.chrom_fst$est_sexM,na.rm=T)


# This is the result reported in paper
#95CIs for analysis with sex as covariate, rather than splitting up analyses
#Significant effect of breast darkness on phenotypic integration (PCIT)
quantile(lm_results_t.chrom_fst$est_avg_t.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_t.chrom_fst$est_avg_t.chrom,na.rm=T,names=F)
#Significant effect of population genetics (Fst) on phenotypic integration (PCIT)
quantile(lm_results_t.chrom_fst$est_Fst,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_t.chrom_fst$est_Fst,na.rm=T,names=F)
#Nonsignificant effect of sex on phenotypic integration (PCIT)
quantile(lm_results_t.chrom_fst$est_sexM,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_t.chrom_fst$est_sexM,na.rm=T)




# This is the result reported in paper
#95CIs for analysis with sex as covariate, rather than splitting up analyses
#Significant effect of breast darkness on phenotypic integration (PCIT)
quantile(lm_results_b.chrom_fst$est_avg_b.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_b.chrom_fst$est_avg_b.chrom,na.rm=T,names=F)
#Significant effect of population genetics (Fst) on phenotypic integration (PCIT)
quantile(lm_results_b.chrom_fst$est_Fst,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_b.chrom_fst$est_Fst,na.rm=T,names=F)
#Nonsignificant effect of sex on phenotypic integration (PCIT)
quantile(lm_results_b.chrom_fst$est_sexM,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_b.chrom_fst$est_sexM,na.rm=T)



# This is the result reported in paper
#95CIs for analysis with sex as covariate, rather than splitting up analyses
#Significant effect of breast darkness on phenotypic integration (PCIT)
quantile(lm_results_v.chrom_fst$est_avg_v.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_v.chrom_fst$est_avg_v.chrom,na.rm=T,names=F)
#Significant effect of population genetics (Fst) on phenotypic integration (PCIT)
quantile(lm_results_v.chrom_fst$est_Fst,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_v.chrom_fst$est_Fst,na.rm=T,names=F)
#Nonsignificant effect of sex on phenotypic integration (PCIT)
quantile(lm_results_v.chrom_fst$est_sexM,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_results_v.chrom_fst$est_sexM,na.rm=T)

###supplemental analysis: do this for all 28 populations by excluding genomic data

boot_f<-res$boot_sum %>% filter(sex=="F") %>% select(population, location, boot_i, sex,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,PINT.c)
boot_m<-res$boot_sum %>% filter(sex=="M") %>% select(population, location, boot_i, sex,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,PINT.c)
boot_lm<-rbind(boot_f,boot_m)

# #Run models with sex as a factor
# lm_res_R<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
#   pops_boot_i <- boot_both %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_r.chrom+sex,data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_r.chrom=mod$estimate[2], est_sex=mod$estimate[3])
# }) %>% bind_rows
# 
# #Throat Chroma
# lm_res_T<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
#   pops_boot_i <- boot_both %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_t.chrom+sex,data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_t.chrom=mod$estimate[2], est_sex=mod$estimate[3])
# }) %>% bind_rows
# 
# #Belly Chroma
# lm_res_B<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
#   pops_boot_i <- boot_both %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_b.chrom+sex,data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_b.chrom=mod$estimate[2], est_sex=mod$estimate[3])
# }) %>% bind_rows
# 
# #Vent Chroma
# lm_res_V<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
#   pops_boot_i <- boot_both %>% filter(boot_i == i)
#   mod<-lm(PINT.c~avg_v.chrom+sex,data=pops_boot_i) %>% tidy()
#   out<-tibble(boot=i,est_avg_v.chrom=mod$estimate[2], est_sex=mod$estimate[3])
# }) %>% bind_rows
# 
# 
# save(lm_res_R, lm_res_T, lm_res_B, lm_res_V, file="Data/PINTanalysis_results_no_Fst.rds")

load("Data/PINTanalysis_results_no_Fst.rds")
#Significant effect of chroma
quantile(lm_res_R$est_avg_r.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_R$est_avg_r.chrom,na.rm=T)

#Nonsignificant effect of sex
quantile(lm_res_R$est_sex,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_R$est_sex,na.rm=T)


# Look at negligible effect of throat on PINT -----------------------------
#Run both in a model with sex as a factor


#Significant effect of chroma
quantile(lm_res_T$est_avg_t.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_T$est_avg_t.chrom,na.rm=T)

#Nonsignificant effect of sex
quantile(lm_res_T$est_sex,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_T$est_sex,na.rm=T)


# Look at effect of belly on PINT -----------------------------
#Run both in a model with sex as a factor


#Significant effect of chroma
quantile(lm_res_B$est_avg_b.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_B$est_avg_b.chrom,na.rm=T)

#Nonsignificant effect of sex
quantile(lm_res_B$est_sex,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_B$est_sex,na.rm=T)

# Look at effect of vent on PINT -----------------------------
#Run both in a model with sex as a factor


#Significant effect of chroma
quantile(lm_res_V$est_avg_v.chrom,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_V$est_avg_v.chrom,na.rm=T)

#Nonsignificant effect of sex
quantile(lm_res_V$est_sex,probs=c(.025,.975),na.rm=T,names=F,type=7)
mean(lm_res_V$est_sex,na.rm=T)

####Figure

#Make males show up on the left to be consistent with other figures.
G_r_simple<-res$mean_traits %>%  mutate(sex=factor(sex, levels=c("M", "F"))) %>%
    ggplot(aes(x = avg_r.chrom, y = PINT.c)) +
    geom_point(size=3,pch=21,col="black", aes(fill = avg_r.chrom)) +
    scale_fill_gradient(
      limits = range(res$mean_traits$avg_r.chrom),
      low = "#FFFFCC",
      high = "#CC6600",
      guide = "none"
    ) + 
    geom_smooth(method="lm", color="black") +
    facet_wrap( ~ sex) + 
    xlab("Breast | Average Population Color (Chroma)")+
    ylab("Phenotypic Integration") +
    theme_bw() +
  ylim(0.75, 3.0) +
    theme(strip.text=element_blank(),axis.text=element_text(size=14), axis.title=element_text(size=16), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) 

G_t_simple<-res$mean_traits %>%  mutate(sex=factor(sex, levels=c("M", "F"))) %>%
  ggplot(aes(x = avg_t.chrom, y = PINT.c)) +
  geom_point(size=3,pch=21,col="black", aes(fill = avg_t.chrom)) +
  scale_fill_gradient(
    limits = range(res$mean_traits$avg_t.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  geom_smooth(method="lm", color="black") +
  facet_wrap( ~ sex) + 
  xlab("Throat | Average Population Color (Chroma)")+
  ylab("Phenotypic Integration") +
  theme_bw() +
  ylim(0.75, 3.0) +
  theme(strip.text=element_blank(),axis.text=element_text(size=14), axis.title=element_text(size=16), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) 

plot_grid(G_r_simple, G_t_simple, nrow=2)
ggsave("figs/Fig3_maleonleft.pdf")



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
t_int_m=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter(sex=="M"))
summary(t_int_m)

t_int_f=lm(edge.weight~mean.t.chrom*edge.type, data=dat2 %>% filter( sex=="F"))
summary(t_int_f)

r_int_m=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter( sex=="M"))
summary(r_int_m)

r_int_f=lm(edge.weight~mean.r.chrom*edge.type, data=dat2 %>% filter( sex=="F"))
summary(r_int_f)

## supplemental figure 2
modplot3m=ggplot(dat2 %>% filter(sex=="M"), aes(x=mean.b.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average belly chroma of population") +
  theme(legend.position="none") +
  ggtitle("Male") +
  guides(linetype=FALSE) 

modplot3f=ggplot(dat2 %>% filter(sex=="F"), aes(x=mean.b.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average belly chroma of population") +
  theme(legend.position="none") +
  ggtitle("Female") +
  guides(linetype=FALSE) 

modplot4m=ggplot(dat2 %>% filter(sex=="M"), aes(x=mean.v.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average vent chroma of population") +
  theme(legend.position="none") +
  ggtitle("Male") +
  guides(linetype=FALSE) 

modplot4f=ggplot(dat2 %>% filter(sex=="F"), aes(x=mean.v.chrom, y=edge.weight, color=patch, fill=patch)) +
  geom_smooth( method="lm", se=F, mapping=aes(linetype=wi_btw)) +
  geom_point(alpha=0.8, pch=21, color="black")+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)+
  theme_cowplot() +
  ylim(0,0.7) +
  ylab("Average edge weight") +
  xlab("Average vent chroma of population") +
  theme(legend.position="none") +
  ggtitle("Female") +
  guides(linetype=FALSE) 

plot_grid(modplot2m, modplot2f, modplot3m, modplot3f, modplot4m, modplot4f,nrow=3, rel_widths=c(2,2,1))

#ggsave("supplemental_figure_1.pdf", width=10, height=12)
##########
#########

## phenotype network plots


#Define f(x) for subsetting data & getting filtered correlation matrix
get_pop_cormat <- function(pop,which_sex,traits){
  d_cor<- d %>% 
    filter(population==pop & sex==which_sex) %>% 
    select(all_of(traits_col)) %>% 
    cor(.,use="pairwise.complete",method = "spear")
  d_cor[diag(d_cor)]<-NA
  
  #Filter algorithm
  # Here, simply ≥|0.3|
  d_cor_bad_indx<-which(abs(d_cor)<0.3)
  d_cor[d_cor_bad_indx]<-0
  
  d_cor
}

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



rawmeansM<-d %>% group_by(population) %>% filter(sex=="M") %>% summarise_at(traits_col,mean,na.rm=T) 

rawmeansF<-d %>% group_by(population) %>% filter(sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)

traits_col
t.lab=c("TBri", "THue", "TChr", "RBri", "RHue", "RChr", "BBri", "BHue", "BChr", "VBri", "VHue", "VChr")
shps=c("triangle", "triangle", "triangle", "circle", "circle", "circle", "square", "square", "square", "diamond", "diamond", "diamond")

pops=unique(d$population)
### Generate male networks figure
#png("figs/Fig 2. Male_10_Networks_ordered.png",width=13,height=6,units="in",res=300)
#pdf("figs/NewFig 2. Male_Networks_modules_all_trial.pdf",width=10,height=14)
par(mfrow=c(7,4),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
scalarM<-sapply(names(rawmeansM)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansM[,x]),q=50 ))) 


#make 50 quantiles for matching color scores
rownames(scalarM)<-rawmeansM$population
scalarM[,c(1:2,4:5,7:8,10:11)] <-51- scalarM[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
#define color ramp with 50 gradations
#nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 
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



par(mfrow=c(7,4),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#Calculate quantiles for each population's color values to color nodes
scalarF<-sapply(names(rawmeansF)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansF[,x]),q=50 ))) 
#make 50 quantiles for matching color scores
rownames(scalarF)<-rawmeansF$population
scalarF[,c(1:2,4:5,7:8,10:11)] <-51- scalarF[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
#define color ramp with 50 gradations
nodepal=brewer.pal(5,"YlGnBu")[c(1,5,2,4,3)]

for (i in 1: length(pops)){
  cur_pop<-pops[i]
  mat<-get_pop_cormat(cur_pop,"F",traits_col)
  g=graph_from_adjacency_matrix(abs(mat), diag=FALSE, weighted=T, mode="undirected")
  nodecolor=nodepal[membership(cluster_fast_greedy(g))]
  
  Q(mat,color=nodecolor,border.color="gray20",labels=t.lab,shape=shps,posCol="#181923",negCol=1,vsize=15,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout=net_layout,rescale=TRUE, maximum=1)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)
  
  # #Add bounding rectangle for Egypt
  # if(cur_pop=="Egypt"){
  #   box(which="figure",lwd=3)
  #   #rect(xleft = -1.6,ybottom = -1.25,xright = 1.25,ytop = 1.6,border="cyan",lwd=3)
  # }
}


