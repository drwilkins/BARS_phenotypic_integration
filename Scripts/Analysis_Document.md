---
title: "Supplemental R Script for Shizuka et al., Global analysis in a widespread songbird reveals that phenotypic integration and modularity evolve along with plumage coloration"
output:
  html_document:
    toc: true
    toc_depth: 3
    keep_md: true
date: "updated 01/25/25 "
---

***

This document contains the R codes used for analyses in Shizuka et al., **Global analysis in a widespread songbird reveals that phenotypic integration and modularity evolve along with plumage coloration.** The analysis codes were first developed by Matthew Wilkins and subsequently edited by Dai Shizuka. 

# Load Libraries

``` r
require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,glue,gtools,PHENIX,dplyr,rsample,pbapply,parallel,lme4,broom, cowplot, RColorBrewer, abind)
library(xtable)
```

# Source Data

## Plumage data


``` r
d=read.csv("../Data/data_for_submission.csv")
d=d%>%dplyr::select(band, population, year, sex, tidyselect::starts_with("t."), starts_with("r."), starts_with("b."), starts_with("v"), lat, long)
```

## Genomic data

This data is organized as pairwise $F_{ST}$ from reference population (Egypt).


``` r
fst<-read.csv("../Data/pairwise_Fst/pairwise_Fst_table.csv")
```

# Bootstrap procedure

The phenotype sample sizes vary between populations. Therefore, we implemented a bootstrap procedure to 

BOOTSTRAP PROCEDURE

Does the following:
  
  1. Runs analysis over the traits specified by 'columns', with IDs in 'id'
  
  2. Iterates this whole procedure over iterate_over column (default="population")
  
  3. Bootstraps data 'boot_n' times (default=1e4), balancing sample sizes over 'strata_var' (default="sex")
      *'strata_levels'= case-insensitive 2 levels of strata,
                        ordered Female first (default=c("f","m"))
    ** Bootstraps are parallelized, change cores with the 'cor_num' parameter.
       default: cor_num=parallel::detectCores()-1
  
  4. Calculates averages for traits specified in 'columns' within each bootstrap sample
  
  5. Calculates the corrected phenotypic integration (PINT.c) value using PHENIX::pint()
  
  6. Calculates dichromatism (or dimorphism) for traits specified in 'columns' and sex defined by 'strata' using dimorphCal() (**calculated in earlier iterations of the analyses but taken out in final version**)
  
  7. Calculates correlations of 'cor_type' (default = "p"; passed to cor()) between:
      a. (PI) and trait mean for "toi"
      b. (PI) and dichromatism across trait(s) in 'toi'
  
  NOTES:
  *   'keep_cols' refers to columns you want to keep in your outputs
  **  'seed' is set to 99 by default, allowing for repeatable bootstrap results.


``` r
# Function Definitions -----------------------------------------------------
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
      
      # B) Sex Differences in Traits (used in earlier interations of analyses but taken out in final draft because it does not explain any patterns)
      
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
```


## Run the bootstrap procedure

Here is how we ran the bootstrapping procedure using the function above (not run here).

``` r
d=read.csv("Data/data_for_submission.csv")

# VERY Time consuming!
# Uncomment if you want to run the full 10e4 bootstraps

res <- boot_analy(
  df = d,
  columns = traits_col,
  toi = c("r.avg.bright","r.chrom","b.chrom","b.avg.bright", "t.chrom"),
  boot_n = 1e4,
  strata_var="sex",
  strata_levels = c("F", "M"),
  keep_cols = c("location", "lat", "long")
  )
```


## Read in the bootstrap results

``` r
res<-readRDS("/Users/dshizuka2/Dropbox/Dai_Research/Main Projects/BARS_phenotypicintegration/results_10k_bootstraps.RDS")
```

***

# Phenotypic Integration analysis

These are the analyses that Table 1 is built from.

## Dataset for populations with phenotype & genomic data.

From the imported datasets, create a dataset that only includes populations with Fst data. Then, separate out male and female data

``` r
#identify the populations with Fst data
pops_w_fst<-fst$population
#trim the bootstrapped phenotype data to these populations
d_gen<-res$boot_sum %>% filter(population %in% pops_w_fst)
#combine the datasets and reduce columns to those we need.
d_gen<-left_join(d_gen,fst[,c("population","weighted_Fst")],by="population") %>% select(population,location,sex,boot_i,PINT,PINT.c,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,weighted_Fst)
#clean up and separate males and females.
d_gen_pops<-d_gen %>% distinct(population) %>% unlist
d_gen_m<-d_gen %>% filter(sex=="M")
d_gen_f<-d_gen %>% filter(sex=="F")
```

### linear models with phenotypic integration as response variable and color (chroma), Fst, and sex as predictor variables. Do this for each patch: throat (t), breast (r), belly (b), vent (v).


``` r
#Run all of the models and save the results so that we don't have to run them each time.
#throat
lm_results_t.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
  pops_boot_i <- d_gen %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_t.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_t.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
}) %>% bind_rows

#breast
lm_results_r.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
  pops_boot_i <- d_gen %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_r.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_r.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
}) %>% bind_rows

#belly
lm_results_b.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
  pops_boot_i <- d_gen %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_b.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_b.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
}) %>% bind_rows

#vent
lm_results_v.chrom_fst <- pbapply::pblapply(1:max(d_gen$boot_i), function(i) {
  pops_boot_i <- d_gen %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_v.chrom+ weighted_Fst+as.factor(sex),data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_v.chrom=mod$estimate[2],est_Fst=mod$estimate[3],est_sexM=mod$estimate[4])
}) %>% bind_rows

save(lm_results_r.chrom_fst, lm_results_t.chrom_fst, lm_results_b.chrom_fst, lm_results_v.chrom_fst, file="Data/PINTanalysis_results_w_Fst.rds")
```

Load in the saved results. 

``` r
load("../Data/PINTanalysis_results_w_Fst.rds")
```

Quick function to extract average coefficient and confidence intervals from bootstrap results.

``` r
get_mean_ci95=function(x) {
  avg=round(mean(x,na.rm=T,names=F), digits=2)
  quant=round(quantile(x,probs=c(.025,.975),na.rm=T,names=F,type=7), digits=2)
  lower.ci=quant[1]
  upper.ci=quant[2]
  paste(avg, " (", lower.ci, " ± ", upper.ci, ")", sep="")
}
```

Extract mean and 95% confidence intervals of coefficient for the linear models and organize into a table. 

``` r
t_mod_out=apply(lm_results_t.chrom_fst[,2:4], 2, get_mean_ci95)
r_mod_out=apply(lm_results_r.chrom_fst[,2:4], 2, get_mean_ci95)
b_mod_out=apply(lm_results_b.chrom_fst[,2:4], 2, get_mean_ci95)
v_mod_out=apply(lm_results_v.chrom_fst[,2:4], 2, get_mean_ci95)

mod.results.df=as.data.frame(matrix(c(t_mod_out, r_mod_out, b_mod_out, v_mod_out), nrow=4, byrow=T))
names(mod.results.df)=c("Chroma", "Fst", "Sex")
rownames(mod.results.df)=c("Throat", "Breast", "Belly", "Vent")

knitr::kable(mod.results.df, caption="Table 1 (note: the column orders are different than in main text)")
```



Table: Table 1 (note: the column orders are different than in main text)

|       |Chroma              |Fst                 |Sex                 |
|:------|:-------------------|:-------------------|:-------------------|
|Throat |0.35 (-8.55 ± 9.35) |5.5 (0.74 ± 10.09)  |0.15 (-0.13 ± 0.43) |
|Breast |3.98 (1.36 ± 6.64)  |6.77 (2.12 ± 11.21) |0.1 (-0.17 ± 0.38)  |
|Belly  |4.04 (1.5 ± 6.66)   |7.42 (2.88 ± 11.81) |0.09 (-0.18 ± 0.37) |
|Vent   |4.58 (1.86 ± 7.26)  |6.5 (1.83 ± 10.91)  |0.08 (-0.19 ± 0.36) |


### Supplemental analysis: do this for all 28 populations by excluding genomic data

The main analysis uses 20 populations for which we have both phenotype AND genomic data. Here, we exclude the genomic data and repeat the analysis with all 28 populations for which we have phenotype data to see if the results are robust (they are).

Build the dataset:


Run all analyses and save results:

``` r
boot_f<-res$boot_sum %>% filter(sex=="F") %>% select(population, location, boot_i, sex,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,PINT.c)

boot_m<-res$boot_sum %>% filter(sex=="M") %>% select(population, location, boot_i, sex,avg_r.chrom,avg_t.chrom,avg_b.chrom,avg_v.chrom,PINT.c)

boot_both<-rbind(boot_f,boot_m)
```


``` r
#Run models with color and sex as covariates
lm_res_R<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
  pops_boot_i <- boot_both %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_r.chrom+sex,data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_r.chrom=mod$estimate[2], est_sex=mod$estimate[3])
}) %>% bind_rows

#Throat Chroma
lm_res_T<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
  pops_boot_i <- boot_both %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_t.chrom+sex,data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_t.chrom=mod$estimate[2], est_sex=mod$estimate[3])
}) %>% bind_rows

#Belly Chroma
lm_res_B<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
  pops_boot_i <- boot_both %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_b.chrom+sex,data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_b.chrom=mod$estimate[2], est_sex=mod$estimate[3])
}) %>% bind_rows

#Vent Chroma
lm_res_V<- pbapply::pblapply(1:max(boot_both$boot_i), function(i) {
  pops_boot_i <- boot_both %>% filter(boot_i == i)
  mod<-lm(PINT.c~avg_v.chrom+sex,data=pops_boot_i) %>% tidy()
  out<-tibble(boot=i,est_avg_v.chrom=mod$estimate[2], est_sex=mod$estimate[3])
}) %>% bind_rows


save(lm_res_R, lm_res_T, lm_res_B, lm_res_V, file="Data/PINTanalysis_results_no_Fst.rds")
```

Re-load the dataset (if run previously)

``` r
load("../Data/PINTanalysis_results_no_Fst.rds")
```


``` r
head(lm_res_T)
```

```
## # A tibble: 6 × 3
##    boot est_avg_t.chrom est_sex
##   <int>           <dbl>   <dbl>
## 1     1          0.464   0.119 
## 2     2          1.27    0.251 
## 3     3         -1.51    0.163 
## 4     4          1.48   -0.0146
## 5     5          0.0744  0.254 
## 6     6          0.995   0.0391
```
Supplemental Table 1: 


``` r
t_mod_nogenome=apply(lm_res_T[,2:3], 2, get_mean_ci95)
r_mod_nogenome=apply(lm_res_R[,2:3], 2, get_mean_ci95)
b_mod_nogenome=apply(lm_res_B[,2:3], 2, get_mean_ci95)
v_mod_nogenome=apply(lm_res_V[,2:3], 2, get_mean_ci95)

mod_nogenome_results.df=as.data.frame(matrix(c(t_mod_nogenome, r_mod_nogenome, b_mod_nogenome, v_mod_nogenome), nrow=4, byrow=T))
names(mod_nogenome_results.df)=c("Chroma", "Sex")
rownames(mod_nogenome_results.df)=c("Throat", "Breast", "Belly", "Vent")

knitr::kable(mod_nogenome_results.df, caption="Table 1 (note: the column orders are different than in main text)")
```



Table: Table 1 (note: the column orders are different than in main text)

|       |Chroma              |Sex                 |
|:------|:-------------------|:-------------------|
|Throat |0.36 (-2.17 ± 2.95) |0.15 (-0.08 ± 0.37) |
|Breast |3.77 (1.36 ± 6.28)  |0.1 (-0.12 ± 0.33)  |
|Belly  |1.95 (0.29 ± 3.73)  |0.12 (-0.1 ± 0.34)  |
|Vent   |3.73 (1.37 ± 6.08)  |0.09 (-0.14 ± 0.31) |



#### Figure 3

``` r
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
```

![](Analysis_Document_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


# Modularity analyses

Analysis for this part of the paper has two parts: 

1. Determine if there are particular sets of traits that act as a separate module. 

    + First, for each population, build a phenotype network, which is essentially the correlation matrix of traits with weakest correlations dropped. 

    + Then, run a 'community detection' algorithm  to determine which nodes form modules in each population. 

    + Add up how many times each pair of traits end up being assigned to the same module.

    + Visualize results as a matrix with heat map colors. 
  
(2) Having determined that throat tends to act as a separate module, we then ask whether the plumage coloration becomes more modular as predicted by Melo & Marroig's model (i.e., stronger correlation within module under selection and weaker correlation between modules).

    + calculate the average correlation (edge weight) within and across modules, with throat patch traits assigned as module 1 and traits in breast/belly/vent patches assigned as module 2. Do this for each population for each sex.
    
    + linear models to estimate interaction effect within vs. between modules. 

## Make male and female datasets 

For these analyses, we will use the raw data (not the bootstrapped data).


``` r
d$population<-as.factor(d$population)
d$sex<-as.factor(d$sex)

#Male data subset by population
data_list_males<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="M")))
names(data_list_males)<-levels(d$population)

#Female data subset by population
data_list_females<-lapply(levels(d$population),function(x) (subset(d,population==x&sex=="F")))
names(data_list_females)<-levels(d$population)
```

## Run pairwise correlations for all traits, by sex.

``` r
traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
traits_col <- traits[-c(1)]

#Male correlations by population
corr_list_males<-lapply(names(data_list_males), function(x) cor(as.matrix(data_list_males[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_males)<-levels(d$population)

#Female correlations by population
corr_list_females<-lapply(names(data_list_females), function(x) cor(as.matrix(data_list_females[[x]][,traits_col]),method="s",use="pairwise.complete"))
names(corr_list_females)<-levels(d$population)
```

## From these correlation matrices, make networks

``` r
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
```

Here, we filter out the lowest 20% of edges and then run fast-greedy community detection on that filtered network. Save which cluster each trait was assigned to, and for each pair of traits, save whether they were assigned to the same module or not into a 'co-membership matrix'. 


``` r
## removing lower 20% of correlations. 
clusters_male=lapply(nets_male, function(x) {
  g=delete.edges(x, which(E(x)$weight<quantile(E(x)$weight, probs=0.2)))
  cluster_fast_greedy(g, weights=E(g)$weight)
})
```

```
## Warning: `delete.edges()` was deprecated in igraph 2.0.0.
## ℹ Please use `delete_edges()` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
clusters_female=lapply(nets_female, function(x) {
  g=delete.edges(x, which(E(x)$weight<quantile(E(x)$weight, probs=0.2)))
  cluster_fast_greedy(g, weights=E(g)$weight)
})

memberships_male=lapply(clusters_male, membership)
comembers_male=lapply(memberships_male, function(x) outer(x, x, "==")+0)

memberships_female=lapply(clusters_female, membership)
comembers_female=lapply(memberships_female, function(x) outer(x, x, "==")+0)
```

## Now, for each sex, stack the co-membership matrix across all populations to form a 3-dimensional array.


``` r
comembers_male_array=abind(comembers_male, along=3)
comembers_female_array=abind(comembers_female, along=3)
```

## Sum the co-membership matrices along populations. The result is a matrix showing how many times each pair of traits are assigned to the same cluster.

``` r
sum_mat_male=apply(comembers_male_array, c(1,2), sum)

sum_mat_female=apply(comembers_female_array, c(1,2), sum)
```

## Plot Figure 4A. 

``` r
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
```

![](Analysis_Document_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

## Avg correlations within and between modules

First, set up a dummy matrix that assigns traits into module 1 (throat color traits) vs. module 2 (breast/belly/vent color traits).

``` r
#throat vs. others
patches=c(rep(1, 3), rep(2, 9))
same.patch=outer(patches, patches, FUN="==")

patch.names=c("Throat", "Breast-Belly-Vent")
modules=matrix(nrow=length(patches), ncol=length(patches))
for(i in 1:2){
  modules[which(patches==i), which(patches==i)] = i
}
modules
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
##  [1,]    1    1    1   NA   NA   NA   NA   NA   NA    NA    NA    NA
##  [2,]    1    1    1   NA   NA   NA   NA   NA   NA    NA    NA    NA
##  [3,]    1    1    1   NA   NA   NA   NA   NA   NA    NA    NA    NA
##  [4,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
##  [5,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
##  [6,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
##  [7,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
##  [8,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
##  [9,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
## [10,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
## [11,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
## [12,]   NA   NA   NA    2    2    2    2    2    2     2     2     2
```

We can now use the above matrix as reference to say which pairwise correlations between traits should fall into one of the three categories (1) within module 1, (2) within module 2, or (3) between modules.

``` r
#throat vs. not throat
mods.list.male=lapply(corr_list_males, function(x) {
  diag(x)=NA
  wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
  wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
  btw_mod=mean(abs(x[which(is.na(modules))]))
  data.frame(wi_mod1, wi_mod2, btw_mod)
})
#organize results into dataframe
mods.dat.male=tibble(bind_rows(mods.list.male))
mods.dat.male$sex="M"
mods.dat.male$population=names(corr_list_males)


#now do the same for females
mods.list.female=lapply(corr_list_females, function(x) {
  diag(x)=NA
  wi_mod1=mean(abs(x[which(modules==1)]), na.rm=T)
  wi_mod2=mean(abs(x[which(modules==2)]), na.rm=T)
  btw_mod=mean(abs(x[which(is.na(modules))]))
  data.frame(wi_mod1, wi_mod2, btw_mod)
})

mods.dat.female=tibble(bind_rows(mods.list.female))
mods.dat.female$sex="F"
mods.dat.female$population=names(corr_list_females)
```

Now we have estimates of within- and between-module correlations. Let's integrate that within the population data.


``` r
#extract population means of color by sex
integ0<-d %>% group_by(population, sex) %>% summarise_at(c("t.chrom","r.chrom","t.avg.bright","r.avg.bright", "b.chrom", "v.chrom", "lat"),mean,na.rm=TRUE) %>% 
  arrange(sex,population) %>% 
  rename(mean.t.chrom=t.chrom,mean.r.chrom=r.chrom,mean.t.avg.bright=t.avg.bright,mean.r.chrom=r.chrom,mean.r.avg.bright=r.avg.bright, mean.b.chrom=b.chrom, mean.v.chrom=v.chrom, latitude=lat)

#now combine the population-level color data with modularity data
mods.dat=bind_rows(mods.dat.female, mods.dat.male)

integ=mods.dat %>% left_join(., integ0) 
```

```
## Joining with `by = join_by(sex, population)`
```

We will now simplify this dataset to columns that we need, then reshape this data to 'longer' format so that we can plot it using ggplot. 

``` r
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
```

### Build Figure 4C,D

``` r
#male plot for throat coloration
modplot1m=ggplot(dat2 %>% filter(sex=="M"), aes(x=mean.t.chrom, y=edge.weight, fill=patch, color=patch)) +
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

#make a version of this that doesn't have the legend
modplot1m_nolegend=modplot1m + theme(legend.position="none")

#female plot for throat color (no legend)
modplot1f_nolegend=ggplot(dat2 %>% filter(sex=="F"), aes(x=mean.t.chrom, y=edge.weight, color=patch, fill=patch)) +
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

#get the legend
legend_plot=get_legend(modplot1m)

plot_grid(modplot1m_nolegend, modplot1f_nolegend,  legend_plot, nrow=1, rel_widths=c(2,2,1))
```

![](Analysis_Document_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

# Phenotype Network Plots (seen in Figure 2)

Custom function for subsetting data & getting filtered correlation matrix.
These plots will have the same layout so that each trait shows up in the same place for all networks. The nodes will be colored based on assignments into clusters (i.e., using the same method as seen in the modularity analysis section 1).



``` r
get_pop_cormat <- function(pop,which_sex,traits){
  d_cor<- d %>% 
    filter(population==pop & sex==which_sex) %>% 
    select(all_of(traits_col)) %>% 
    cor(.,use="pairwise.complete",method = "spear")
  d_cor[diag(d_cor)]<-NA
  
  #Filter algorithm
  # Here, simply ≥|0.2|
  d_cor_bad_indx<-which(abs(d_cor)<0.2)
  d_cor[d_cor_bad_indx]<-0
  
  d_cor
}
```

Network plotting function using `qgraph` and also define network layout to follow the hypothetical plot in Figure 1.

``` r
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
```

Set up parameters for trait means for each sex/population, and what the trait labels and shapes are. 

``` r
rawmeansM<-d %>% group_by(population) %>% filter(sex=="M") %>% summarise_at(traits_col,mean,na.rm=T) 

rawmeansF<-d %>% group_by(population) %>% filter(sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)


t.lab=c("TBri", "THue", "TChr", "RBri", "RHue", "RChr", "BBri", "BHue", "BChr", "VBri", "VHue", "VChr")
shps=c("triangle", "triangle", "triangle", "circle", "circle", "circle", "square", "square", "square", "diamond", "diamond", "diamond")

pops=unique(d$population)
```

### Plot networks for males for all populations

``` r
par(mfrow=c(7,4),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#color palette
nodepal=brewer.pal(5,"YlGnBu")[c(1,5,2,4,3)]

for (i in 1: length(pops)){
  cur_pop<-pops[i]
  mat<-get_pop_cormat(cur_pop,"M",traits_col)
  g=graph_from_adjacency_matrix(abs(mat), diag=FALSE, weighted=T, mode="undirected")
  nodecolor=nodepal[membership(cluster_fast_greedy(g))]

  Q(abs(mat),color=nodecolor,border.color="gray20",labels=t.lab,shape=shps,posCol="#181923",negCol=1,vsize=15,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout=net_layout,rescale=TRUE, maximum=1)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)
  
}
```

![](Analysis_Document_files/figure-html/unnamed-chunk-33-1.png)<!-- -->
### Plot networks for females for all populations

``` r
par(mfrow=c(7,4),mar=rep(3,4),xpd=T,oma=rep(1,4),ps=18)

#color palette
nodepal=brewer.pal(5,"YlGnBu")[c(1,5,2,4,3)]

for (i in 1: length(pops)){
  cur_pop<-pops[i]
  mat<-get_pop_cormat(cur_pop,"F",traits_col)
  g=graph_from_adjacency_matrix(abs(mat), diag=FALSE, weighted=T, mode="undirected")
  nodecolor=nodepal[membership(cluster_fast_greedy(g))]
  
  Q(mat,color=nodecolor,border.color="gray20",labels=t.lab,shape=shps,posCol="#181923",negCol=1,vsize=15,lab.col="#181923",lab.font=2,lab.scale=F,label.cex=.7,label.scale.equal=T,layout=net_layout,rescale=TRUE, maximum=1)
  
  mtext(cur_pop,3,line=.6,at=-1.4,adj=0,col="#181923",cex=.6,font=2)

}
```

![](Analysis_Document_files/figure-html/unnamed-chunk-34-1.png)<!-- -->
