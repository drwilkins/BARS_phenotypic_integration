require(pacman)
p_load(tidyverse,qgraph,igraph,devtools,patchwork,ggrepel,ggiraph,glue,ggnetwork,gtools,colourvalues,PHENIX,dplyr,rsample,pbapply,parallel)
remotes::install_github("galacticpolymath/galacticEdTools")
require(galacticEdTools)


# Function Definitions -----------------------------------------------------
# CALCULATE DICHROMATISM
# Code from Wilkins et al. "Analysis of female song provides insight..." AnBeh 2020
# Calculates dimorphism as (Mean F - Mean M)/Pooled SD

#df= a data frame/tibble 
#columns= columns you want to calculate sex differences on
#strata_var= column of binary categories (default="sex)
#strata_levels= case-insensitive 2 levels of strata, ordered Female first (default=c("f","m"))

dimorphCal <- function(df, columns,strata_var="sex",strata_levels=c("f","m"),keep_cols=NULL) {
    #make case-insensitive levels for filtering
    stratum1<-c(toupper(strata_levels[1]),tolower(strata_levels[1]))
    stratum2<-c(toupper(strata_levels[2]),tolower(strata_levels[2]))
    
    #Calculate pooled standard deviations for all traits & strata_levels
    pooledSD <- sapply(1:length(columns), function(i)
    {
      df_i<-df[,c({{strata_var}}, columns[i])]
      sd.f =df_i %>% dplyr::filter(eval(as.symbol(strata_var))%in%stratum1) %>% summarize_if(is.double,sd,na.rm=T) 
      sd.m = df_i %>% dplyr::filter(eval(as.symbol(strata_var))%in%stratum2) %>% summarize_if(is.double,sd,na.rm=T) 
      sqrt((sd.f[columns[i]] ^ 2 + sd.m[columns[i]] ^ 2) / 2)
    })
  
  dfmeans <-as_tibble(df) %>% dplyr::select(all_of(c(strata_var,columns,keep_cols))) %>% 
    group_by(!!!syms(strata_var),!!!syms(keep_cols)) %>% 
    summarize_if(is.double, mean, na.rm =T)
  
  diffs<-(dfmeans[which(unlist(dfmeans[,strata_var]) %in% stratum1), columns] - 
      dfmeans[which(unlist(dfmeans[,strata_var]) %in% stratum2), columns]) / pooledSD
  #output combined info w/ keep_cols if supplied
  tibble(dfmeans[1,keep_cols],diffs)
    
}

# PRIMARY ANALYSIS
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
    
  





# 1. Setup -------------------------------------------------------------------
#Import all data
d0<-read_csv("Data/all_populations.csv")
nrow(d0)

traits<-c('tail.mean','t.avg.bright','t.hue','t.chrom','r.avg.bright','r.hue','r.chrom','b.avg.bright','b.hue','b.chrom','v.avg.bright','v.hue','v.chrom')
#Just the color traits
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


# 2.  Analyze -------------------------------------------------------------
# 

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
res<-readRDS("results_10k_bootstraps.RDS")

# res_no_hyb <-boot_analy(
#   df = d %>% filter(hybrid_zone=="no"),
#   columns = traits_col,
#   toi = c("r.chrom", "t.chrom"),
#   boot_n = 1000,
#   strata_levels = c("F", "M"),
#   keep_cols = c("location", "lat", "long")
# )
# 


# 3. Graph ----------------------------------------------------------------------

# Fig 1. Graph of Phenotypic Integration ~ Avg. Color for Breast and Throat -----
showtext::showtext_opts(dpi=300) #important b/c of a bug
mytheme<-galacticEdTools::theme_galactic(
     base.theme = "linedraw",
     grid.wt.maj = .1,
     grid.wt.min = 0,
     grid.col = "gray70",
    text.cex = 0.7,
    pad.outer = rep(5,4)
   )+theme(strip.text = element_text(size=12))

#breast patch graph
(G_r<-res$mean_traits %>%  
      ggplot(aes(x = avg_r.chrom, y = PINT.c)) +
  stat_ellipse(col="gray60",size=.5) +
  mytheme+
  geom_point(size=3,pch=21,col="black", aes(fill = avg_r.chrom)) +
  scale_fill_gradient(
    limits = range(res$mean_traits$avg_r.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
    ggrepel::geom_text_repel(aes(label =location),col="black", max.overlaps = 15,size=2.5,segment.size=.25,force = 15,min.segment.length = .1,box.padding = 0.9,seed = 100)+
  xlab("Breast | Average Population Darkness (Chroma)")+
  ylab("Phenotypic Integration")
  )

#throat patch graph with PINT (Wagner 1984 method for phenotypic integration)
(G_t<-res$mean_traits %>%  
    group_by(sex) %>% 
      ggplot(aes(x = avg_t.chrom, y = PINT.c,
                 group=sex)) +
  mytheme+
  stat_ellipse(col="gray60",size=.5)+
  geom_point(size=3,pch=21,col="black", aes(fill = avg_t.chrom)) +
  scale_fill_gradient(
    limits = range(res$mean_traits$avg_t.chrom),
    low = "#FFFFCC",
    high = "#CC6600",
    guide = "none"
  ) + 
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
    ggrepel::geom_text_repel(aes(label =location),col="black", max.overlaps = 15,size=2.5,segment.size=.25,force = 15,min.segment.length = .1,box.padding = 0.8,seed = 100)+
  xlab("Throat | Average Population Darkness (Chroma)")+
  ylab("Phenotypic Integration")
  )

#patchwork syntax
(G_combined<-G_t/G_r)
ggsave("figs/Fig 1. PINT ~ breast + throat chroma.png",dpi=300,width=8,height=8)


# Fig.  2. Graph of PINT ~ Sex Difference in Breast Chroma ----------------
res$mean_traits %>% left_join(., res$mean_sex_diffs %>% select(population,SD_r.chrom)) %>% 
  ggplot(aes(x=SD_r.chrom,y=PINT.c))+
  mytheme+
  stat_ellipse(col="gray60",size=.5)+
  geom_point()+
  facet_wrap( ~ sex,labeller =as_labeller(c(M="Males",F="Females") )) + 
    ggrepel::geom_text_repel(aes(label =location),col="black", max.overlaps = 15,size=2.5,segment.size=.25,force = 14,min.segment.length = 0.1,box.padding = 0.8,seed = 100)+
  labs(x=expression(atop(bold(Dichromatism~"in"~Breast~Chroma),"\u2190 Darker Males than Females")),
       y=expression(bold("Phenotypic Integration")))
ggsave("figs/Fig 2. Phenotypic Integration ~ Dichromatism in Breast Chroma.png",dpi=300,width=8,height=4)


# Fig. 3.  Phenotype networks for key populations -------------------------


# Function Definitions ----------------------------------------------------
#Function to calculate raw phenotypic integration (we've only done bootstraps)
# do each sex separately;
# traits_of_interest= vector of traits to use to calculate phenotypic integration
# trait_to_average= which trait to average and output (that you want to compare to phenotypic integration)
# ordered=logical; should results be ordered descending by phenotypic integration?
####>>>>>>>>>>>>>>>>>>>>>
PI_for_pops <-
  function(d,
           which_pops,
           which_sex,
           traits_of_interest = traits_col,
           trait_to_average,
           ordered=TRUE) {
    d1 <-
      d %>%
      dplyr::filter(sex == which_sex, population %in% which_pops)
    #iterate over target populations for sex==which_sex
    res <- lapply(which_pops, function(x_i) {
      d2 <- d1 %>%
        dplyr::filter(population == x_i) 
      res_i<-d2%>%
        dplyr::select(all_of(traits_of_interest)) %>%
        na.omit() %>%
        PHENIX::pint()
      #output tibble
      dplyr::tibble(
        sex = which_sex,
        population = x_i,
        N_PINT = res_i$N,
        PINT = res_i$PINT,
        PINT.c = res_i$PINT.c,
        "avg_{trait_to_average}":=mean(unlist(d2[,trait_to_average]),na.rm=T)
      )
    }) %>% bind_rows()
    if(ordered){
    res %>% dplyr::arrange(desc(PINT.c))
    }else{res}
    
  }

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

## Make custom plot function
Q<-function(COR,
            ...) {
  G <-
    qgraph(
      COR,
      mar=rep(3,4),
      diag = F,
      fade = F,
      label.norm = "0000",
      border.color = "gray20",
      shape = "triangle",
      posCol = "#181923",
      negCol = 1,
      vsize = 15,
      label.cex = 1.1,
      label.font=2,
      label.scale.equal = T,
      layout = "circle",
      rescale=T,
      aspect=T,
      ...
    )
return(G)}
#<<<<<<<<<<<<<<<
#


#look at all included populations to figure out which pops to showcase
  #Males
  res$mean_traits %>% filter(sex=="M")%>% 
    arrange(desc(PINT.c)) %>% print(n=56)
  #Lowest for male= "egypt", but there's something whack going on there (males have extremely low and females extremely high phenotypic integration, which is atypical and possibly from reinforcement), so we'll use Taiwan, Middle= "morocco",High="baotu"
  
  # #Females
  # res$mean_traits %>% filter(sex=="F")%>% 
  # arrange(desc(PINT.c)) %>% print(n=56)
  # 

#which ones are we interested in?
pops_of_interest<-c("baotu","taiwan","morocco")

# Calculate PINT for each pop by sex,
# ordered by descending PINT.c
#  --------------------------------------

male_pops <- PI_for_pops(d,pops_of_interest,"M",traits_col,trait_to_average = "r.chrom")
female_pops <- PI_for_pops(d,pops_of_interest,"F",traits_col,trait_to_average = "r.chrom")
#order females like males
female_pops<-female_pops[match(male_pops$population,female_pops$population),]

#short labels for phenonets
traits_col
net_labs<-c("TBri","THue","TChr","RBri","RHue","RChr","BBri","BHue","BChr","VBri","VHue","VChr")
#Get means for traits in each population for each sex
rawmeansM<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="M") %>% summarise_at(traits_col,mean,na.rm=T)

rawmeansF<-d %>% group_by(population) %>% filter(population %in% pops_of_interest,sex=="F") %>% summarise_at(traits_col,mean,na.rm=T)


### Generate male networks figure
png("figs/Fig 3. 6_Networks_ordered.png",width=11,height=6,units="in",res=300)
par(xpd=T,oma=rep(1,4),ps=18,mar=rep(3,4))
#create somewhat complex layout to have titles and graphs together
l<-layout(matrix(c(1,7,2,8,3,9,4,10,5,11,6,12),nrow=2),widths=rep(rep(c(0.3,0.7),3),2))
# layout.show(l)
#Calculate quantiles for each population's color values to color nodes
  scalar<-sapply(names(rawmeansM)[-1],function(x) as.numeric(gtools::quantcut(unlist(rawmeansM[,x]),q=50 ))) 
  #make 50 quantiles for matching color scores
  rownames(scalar)<-rawmeansM$population
  scalar[,c(1:2,4:5,7:8,10:11)] <-51- scalar[,c(1:2,4:5,7:8,10:11)]  #reverse brightness & hue measures so lower values are darker
  #define color ramp with 50 gradations
  nodepal<-colorRampPalette(c("#FFFFCC","#CC6600"),interpolate="spline")(50) 

#Plot Male phenonets
for (i in 1: nrow(male_pops)){
  cur_pop<-male_pops$population[i]
  mat<-get_pop_cormat(cur_pop,"M",traits_col)
  nodecolor<-nodepal[scalar[as.character(cur_pop),]]
  
  #Plot info before network
  plot.new()
  long_name_i<-d %>% distinct(population,.keep_all = T) %>% 
              filter(population==cur_pop) %>% select(location) %>% unlist
  mtext(long_name_i,
        3,line=.6,at=-1.4,
        adj=0,col="#181923",cex=.7,font=2)
  mtext(paste0("Phenotypic Integration: ",
               male_pops$PINT.c[i]),
        3,line=-.75,at=-1.4,
        adj=0,col="#181923",cex=.45,font=1)
  mtext(paste0("Mean Breast Chroma: ",
               round(male_pops$avg_r.chrom[i],3)),
        3,line=-1.75,at=-1.4,
        adj=0,col="#181923",cex=.45,font=1)
  
  Q(mat,color=nodecolor,labels=net_labs)
  
  
}
  #Plot Female Phenonets
  for (i in 1: nrow(female_pops)){
  cur_pop<-female_pops$population[i]
  mat<-get_pop_cormat(cur_pop,"F",traits_col)
  nodecolor<-nodepal[scalar[as.character(cur_pop),]]
 # groupings<-list(throat=1:3,breast=4:6,belly=7:9,vent=10:12)
  #Plot info before network
  plot.new()
  long_name_i<-d %>% distinct(population,.keep_all = T) %>% 
              filter(population==cur_pop) %>% select(location) %>% unlist
  mtext(long_name_i,
        3,line=.6,at=-1.4,
        adj=0,col="#181923",cex=.7,font=2)
  mtext(paste0("Phenotypic Integration: ",
               female_pops$PINT.c[i]),
        3,line=-.75,at=-1.4,
        adj=0,col="#181923",cex=.45,font=1)
  mtext(paste0("Mean Breast Chroma: ",
               round(female_pops$avg_r.chrom[i],3)),
        3,line=-1.75,at=-1.4,
        adj=0,col="#181923",cex=.45,font=1)
  #plot female phenonet
  Q(mat,color=nodecolor,labels=net_labs)
 
  

  
}
dev.off()

