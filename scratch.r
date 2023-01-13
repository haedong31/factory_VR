library(tidyverse)
library(lubridate)

alog <- read_csv("./sim_data/assembly_line_4players8.csv")
tframe <- 1
players <- str_c("player",seq(1,4))
acts <- unique(alog$action)
col_names <- c("t","s1","t1","w1","s2","t2","w2",
               "s3","t3","w3","s4","t4","w4")
num_iters <- ceiling(tail(alog$time_stamp, 1)/tframe)

for (i in 1:num_iters) {
  t0 <- (i-1)*tframe
  t1 <- i*tframe
  sub_alog <- alog[findInterval(alog$time_stamp, c(t0,t1), left.open = TRUE) == 1L,]
  
  r <- vector('list', length = length(col_names))
  names(r) <- col_names
  if (nrow(sub_alog) == 0) {
    r[[1]] <- t1
    if (i == 1) {
      r[-1] <- 0
      etbl <- add_row(etbl, !!!r)
    } else {
      r[-1] <- etbl[i-1,-1] %>% as.list() # no update
      etbl <- add_row(etbl, !!!r)
    }
  } else {
    for (j in seq_along(players)) {
      if (any(sub_alog$player == players[j])) {
        sub_alog_pj <- sub_alog[sub_alog$player == players[j],]
        num_acts <- nrow(sub_alog_pj)
        if (num_acts == 1) {
          # save singletons
          r[[1]] <- t1
        } else {
          # generate edges
          r[[1]] <- rep(t1, num_acts-1)
          efrom <- vector('numeric', length = (num_acts-1))
          eto <- vector('numeric', length = (num_acts-1))
          
          aid <- match(sub_alog_pj$action, acts)
          efrom <- aid[1:(num_acts-1)]
          eto <- aid[2:num_acts]
          r[[2*j]] <- efrom # source
          r[[2*j+1]] <- eto # target
        }
      }
    }
  }
}