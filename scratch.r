library(tidyverse)
library(lubridate)

# input arguments
alog <- read_csv("./sim_data/assembly_line_4players8.csv")
tframe <- 1
num_players <- 4

players <- str_c("player",seq(1,num_players))
# col_names <- c("t","s1","t1","s2","t2","s3","t3","s4","t4")
col_names <- "t"
for (i in 1:num_players) {
  col_names <- c(col_names, str_c(c("s","t"), i))
}
acts <- unique(alog$action)
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
      etbl <- add_row(etbl, !!!r)
    } else {
      etbl <- add_row(etbl, !!!r)
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
      etbl <- tibble(!!!r)
    } else {
      etbl <- add_row(etbl, !!!r)
    }
  }
}