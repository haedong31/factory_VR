library(tidyverse)
library(lubridate)
library(igraph)

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

df <- bind_rows(p1, p2)
df <- df %>% arrange(time)

##### Node & edge lists (undirected) -----
# node list
num_acts <- length(unique(df$action))
players <- unique(df$player)
num_nodes <- num_acts * length(players)
nid <- seq(1, num_nodes)
nlabel <- vector("character", length = num_nodes)
for (i in seq_along(players)) {
  nlabel[(1+num_acts*(i-1)):(i*num_acts)] <- players[i]
}
nlist <- tibble(Id = nid, Label = nlabel)

# edge list
num_edges <- nrow(df)-1
efrom <- vector(mode = "character", length = num_edges)
eto <- vector(mode = "character", length = num_edges)
eid <- vector("list", length = num_nodes^2)
for (i in 1:num_nodes) {
  eid[(1+num_nodes*(i-1)):(i*num_nodes)] <- 0
  names(eid)[(1+num_nodes*(i-1)):(i*num_nodes)] <- str_c(as.character(i),nid)
}
w <- vector("numeric", length = num_edges)
t <- df$time
t <- str_c(hour(t),":",minute(t),":",second(t))

for (i in 1:num_edges) {
  if (df$action[i]=="Spawn Pieces" & df$player[i]==players[1]) {
    efrom[[i]] <- "1"
  } else if (df$action[i]=="Teleport Part" & df$player[i]==players[1]) {
    efrom[[i]] <- "2"
  } else if (df$action[i]=="Spawn Pieces" & df$player[i]==players[2]) {
    efrom[[i]] <- "3"
  } else {
    efrom[[i]] <- "4"
  }
  
  if (df$action[i+1]=="Spawn Pieces" & df$player[i+1]==players[1]) {
    eto[[i]] <- "1"
  } else if (df$action[i+1]=="Teleport Part" & df$player[i+1]==players[1]) {
    eto[[i]] <- "2"
  } else if (df$action[i+1]=="Spawn Pieces" & df$player[i+1]==players[2]) {
    eto[[i]] <- "3"
  } else {
    eto[[i]] <- "4"
  }
  
  e1 <- efrom[[i]]
  e2 <- eto[[i]]
  if (e1 == e2) {
    eid[[str_c(e1, e2)]] <- eid[[str_c(e1, e2)]] + 1
  } else{
    eid[[str_c(e1, e2)]] <- eid[[str_c(e1, e2)]] + 1
    eid[[str_c(e2, e1)]] <- eid[[str_c(e2, e1)]] + 1
  }
  w[i] <- eid[[str_c(e1, e2)]]
}

el <- tibble(Source = efrom, Target = eto, Weight = w, Time = t[1:num_edges])

##### Network stats -----
