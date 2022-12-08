library(tidyverse)
library(lubridate)

gen_nlist <- function(log_df) {
  players <- unique(log_df$player)
  acts <- unique(log_df$action)
  num_acts <- length(acts)
  num_nodes <- num_acts*length(players)
  
  nlabel <- vector("character", length = num_nodes)
  for (i in seq_along(players)){
    nlabel[(1+num_acts*(i-1)):(i*num_acts)] <- str_c(players[i], " & ", acts)
  }
  
  nid <- seq(1, num_nodes)
  ngroup <- vector("character", length = num_nodes)
  for (i in seq_along(players)) {
    ngroup[(1+num_acts*(i-1)):(i*num_acts)] <- players[i]
  }
  return(tibble(Id = nid, Label = nlabel, Group = ngroup))
}

gen_elist <- function(log_df, nlist) {
  num_edges <- nrow(log_df)-1
  efrom_list <- vector(mode = "character", length = num_edges)
  eto_list <- vector(mode = "character", length = num_edges)
  d <- vector(mode = "numeric", length = num_edges)
  
  for (i in 1:num_edges) {
    efrom <- str_c(log_df$player[i], " & ", log_df$action[i])
    efrom_idx <- which(str_detect(nlist$Label, efrom))
    efrom_list[i] <- nlist$Id[efrom_idx]
    
    eto <- str_c(log_df$player[i+1], " & ", log_df$action[i+1])
    eto_idx <- which(str_detect(nlist$Label, eto))
    eto_list[i] <- nlist$Id[eto_idx]
    
    d[i] <- log_df$time_stamp[i+1]-log_df$time_stamp[i]
  }
  return(tibble(Source = efrom_list, Target = eto_list, Duration = d))
}

##### Simulated data -----
log_df <- read_csv("./sim_data/assembly_line_4players.csv")
nlist <- gen_nlist(log_df)
elist <- gen_elist(log_df, nlist)


##### Data import and preparation -----
p1 <- read_delim("./data/data_haedong_orig.txt", 
                 delim = ",", col_names = FALSE)
colnames(p1) <- c("action", "time")
p1$time <- str_remove(p1$time, " PM")
p1$time <- mdy_hms(p1$time)
p1 <- p1 %>% mutate(player = "player1")
p1 <- p1 %>% select(action, player, time)

p2 <- read_delim("./data/data_trevor_orig.txt",
                 delim = ",", col_names = FALSE)
colnames(p2) <- c("action", "time")
p2$time <- str_remove(p2$time, " PM")
p2$time <- mdy_hms(p2$time)
p2 <- p2 %>% mutate(player = "player2")
p2 <- p2 %>% select(action, player, time)

log_df <- bind_rows(p1, p2)
log_df <- log_df %>% arrange(time)
log_df$time <- str_c(hour(log_df$time),":",minute(log_df$time),":",second(log_df$time))

##### Node & edge lists (undirected) -----
# node list
nlist <- gen_nlist(log_df)
elist <- gen_elist(log_df, nlist)
# write_csv(nlist, "./data/nlist.csv")
# rite_csv(el_dynamic, "./data/elist_dynamic.csv")
