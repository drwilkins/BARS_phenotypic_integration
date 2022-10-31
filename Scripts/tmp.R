require(dplyr); require(readr)
temp<-read_csv("Data/temp.csv")

new <- temp %>% distinct(pref_name,.keep_all = T) %>% select(pref_name,`hybrid zone?`,population) %>% rename(location=pref_name,
                                                                                                     hybrid_zone=`hybrid zone?`)
write_csv(new,"Data/location<->population_key.csv")
