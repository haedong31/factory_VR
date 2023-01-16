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

gen_elist2 <- function(alog, tframe, num_players) {
  players <- str_c("player",seq(1,num_players))
  col_names <- "t"
  for (i in 1:num_players) {
    col_names <- c(col_names, str_c(c("s","t"), i))
  }
  acts <- unique(alog$action)
  nid <- seq(1, length(acts))
  nlist <- tibble(ID = nid, Label = acts)
  
  num_iters <- ceiling(tail(alog$time_stamp, 1)/tframe)
  single_stk <- rep(0, length(players))
  names(single_stk) <- players
  for (i in 1:num_iters) {
    # for a time window
    t0 <- (i-1)*tframe
    t1 <- i*tframe
    sub_alog <- alog[findInterval
                    (alog$time_stamp, c(t0,t1), 
                      left.open = TRUE, rightmost.closed = TRUE) == 1L,]
    
    r <- vector('list', length = length(col_names))
    names(r) <- col_names
    if (nrow(sub_alog) == 0) {
      # no observations in the time window -> no update (coded by 0)
      r[[1]] <- t1
      r[-1] <- 0
      if (i == 1) {
        elist <- add_row(elist, !!!r)
      } else {
        elist <- add_row(elist, !!!r)
      }
    } else {
      # maximum number of edges 
      num_edges <- rep(0, length(players))
      for (j in seq_along(players)) {
        sub_alog_pj <- sub_alog[sub_alog$player == players[j],]
        if (single_stk[j] != 0) {
          num_edges[j] <- nrow(sub_alog_pj)
        } else {
          num_edges[j] <- max(nrow(sub_alog_pj)-1, 0)
        }
      }
      max_num_edges <- max(num_edges)
      r[[1]] <- rep(t1, max_num_edges)
      
      for (j in seq_along(players)) {
        efrom <- rep(0, max_num_edges)
        eto <- rep(0, max_num_edges)
        
        if (any(sub_alog$player == players[j])) {
          # for player j
          sub_alog_pj <- sub_alog[sub_alog$player == players[j],]
          num_acts <- nrow(sub_alog_pj)
          
          aid <- match(sub_alog_pj$action, acts)
          if (num_acts == 1) {
            # singleton case
            if (single_stk[j] == 0) {
              single_stk[j] <- aid
              r[[2*j]] <- efrom
              r[[2*j+1]] <- eto
            } else {
              efrom[1] <- single_stk[j]
              eto[1] <- aid
              single_stk[j] <- 0
              r[[2*j]] <- efrom
              r[[2*j+1]] <- eto
            }
          } else {
            # exist >= 2 events
            if(single_stk[j] == 0) {
              efrom[1:num_edges[j]] <- aid[1:(num_acts-1)]
              eto[1:num_edges[j]] <- aid[2:num_acts]
              r[[2*j]] <- efrom
              r[[2*j+1]] <- eto
            } else {
              efrom[1:num_edges[j]] <- c(single_stk[j], aid[1:length(aid)-1])
              eto[1:num_edges[j]] <- aid
              single_stk[j] <- 0
              r[[2*j]] <- efrom
              r[[2*j+1]] <- eto
            }
          }
        } else {
          # empty for player j
          r[[2*j]] <- efrom
          r[[2*j+1]] <- eto
        }
      }
      if (i == 1) {
        elist <- tibble(!!!r)
      } else {
        elist <- add_row(elist, !!!r)
      }
    }
  }
  o <- vector('list', length = 2)
  o[[1]] <- nlist
  o[[2]] <- elist
  return(o)
}

##### Simulation data (8 orders) / separate graphs -----
# assembly
alog <- read_csv("./sim_data/assembly_line_4players8.csv")
tframe <- 1
num_players <- 4
g <- gen_elist2(alog, tframe, num_players)
nlist <- g[[1]]
elist <- g[[2]]
write_csv(nlist, "./sim_data/assembly_line_4players8_nlist.csv")
write_csv(elist, "./sim_data/assembly_line_4players8_elist.csv")

# craft
alog <- read_csv("./sim_data/craft_4players8.csv")
tframe <- 1
num_players <- 4
g <- gen_elist2(alog, tframe, num_players)
nlist <- g[[1]]
elist <- g[[2]]
write_csv(nlist, "./sim_data/craft_4players8_nlist.csv")
write_csv(elist, "./sim_data/craft_4players8_elist.csv")

##### Simulated data 1 -----
log_df <- read_csv("./sim_data/assembly_line_4players.csv")
nlist <- gen_nlist(log_df)
elist <- gen_elist(log_df, nlist)
write_csv(nlist, "./sim_data/assembly_line_4players_nlist.csv")
write_csv(elist, "./sim_data/assembly_line_4players_elist.csv")

##### Simulated data 2 -----
log_df <- read_csv("./sim_data/craft_4players.csv")
nlist <- gen_nlist(log_df)
elist <- gen_elist(log_df, nlist)
write_csv(nlist, "./sim_data/craft_4players_nlist.csv")
write_csv(elist, "./sim_data/craft_4players_elist.csv")

##### Real (toy) data -----
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
log_df <- rename(log_df, time_stamp = time)

nlist <- gen_nlist(log_df)
elist <- gen_elist(log_df, nlist)
# write_csv(nlist, "./data/nlist.csv")
# write_csv(el_dynamic, "./data/elist_dynamic.csv")
