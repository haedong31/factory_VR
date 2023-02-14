# -*- coding: utf-8 -*-
#%% Classes and functions
import pandas as pd
import numpy as np
import RNG
from pathlib import Path

class ActionLogger(object):
    def __init__(self):
        self.clock = 0.0
        self.n = 0
        self.alog = []
        self.tlog = []

    def write_act(self,a,t):
        self.clock += t
        self.n += 1
        self.alog.append(a)
        self.tlog.append(self.clock)

def learning_curve(a,b,x):
    return a*np.power(x,-b)

def make_wheelset(logger,exp_means,rnd_seed):
    # sourcing
    for _ in range(7):
        logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))

    # (rim + tire) x 2
    for _ in range(2):
        logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed)) # rim
        logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed)) # tire
        logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # rim + tire

    # (wheels + axle)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed)) # axle
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # axle + wheel 1
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # axle + wheel 2
    
    # (axle + axle)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed)) # connecting plate
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # axle + axle

def make_chassis_base(logger,exp_means,rnd_seed):
    # sourcing
    for _ in range(4):
        logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))
    
    # (long plate + long plate)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

    # (connecting plates) x 2
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def make_front(logger,exp_means,rnd_seed):
    for _ in range(5):
        logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))

    # (front bumper)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # frist two blocks
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # last blocok

    # (engine bay)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # slope 1
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed)) # slope 2 

def finish_station2(logger,exp_means,rnd_seed):
    # (wheelsets + base)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))        

    # assemble base and front module
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def make_steering(logger,exp_means,rnd_seed):
    for _ in range(2):
        logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def finish_station3(logger,exp_means,rnd_seed):
    # windshield
    logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

    # (chassis + steering module)
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def make_body_frame(logger,exp_means,rnd_seed):
    for _ in range(7):
        logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def finish_station4(logger,exp_means,rnd_seed):
    # ceiling
    logger.write_act("SpawnPart",RNG.Expon(exp_means[0],rnd_seed))
    logger.write_act("TeleportPart",RNG.Expon(exp_means[1],rnd_seed))
    logger.write_act("Assemble",RNG.Expon(exp_means[2],rnd_seed))

def craft_car(logger,exp_means,rnd_seed):
    make_wheelset(logger,exp_means,rnd_seed)
    make_wheelset(logger,exp_means,rnd_seed)
    make_chassis_base(logger,exp_means,rnd_seed)
    make_front(logger,exp_means,rnd_seed)
    finish_station2(logger,exp_means,rnd_seed)
    make_steering(logger,exp_means,rnd_seed)
    finish_station3(logger,exp_means,rnd_seed)
    make_body_frame(logger,exp_means,rnd_seed)
    finish_station4(logger,exp_means,rnd_seed)
    
#%% Number of orders, initial mean times, and learning curve parameters
exp_num = 'exp3'
save_dir = Path.cwd()/'sim_data'/exp_num
save_dir.mkdir(parents=True,exist_ok=True)
num_order = 8

# initial mean times of actions
spawn_meant0 = 0.5
teleport_meant0 = 0.5
assemble_meant0 = 0.7
pass_meant0 = 0.25

# learning curves for craft and mass production
np.random.seed(241)
lc_rates = [0.6,0.7]
r1 = np.random.uniform(lc_rates[0]-lc_rates[0]*0.1,lc_rates[0]+lc_rates[0]*0.1,4)
r2 = np.random.uniform(lc_rates[1]-lc_rates[1]*0.1,lc_rates[1]+lc_rates[1]*0.1,4)
b1 = -(np.log(r1)/np.log(2))
b2 = -(np.log(r2)/np.log(2))

#%% Mass production 
# Station 1 -  Wheelsets x2
logger1 = ActionLogger()
for i in range(num_order):
    t = [learning_curve(spawn_meant0,b1[0],i+1), 
         learning_curve(teleport_meant0,b1[0],i+1), 
         learning_curve(assemble_meant0,b1[0],i+1)]
    make_wheelset(logger1,t,1)
    make_wheelset(logger1,t,1)    
    logger1.write_act("Pass",learning_curve(pass_meant0,b1[0],i+1))

# Station 2 - Chassis
upstream_parts_idx = np.where(np.array(logger1.alog)=="Pass")[0]
upstream_parts_time = np.array(logger1.tlog)[upstream_parts_idx]

logger2 = ActionLogger()
num_received_up = 0
num_base = 0
num_complete = 0

# first chassis base
t = [learning_curve(spawn_meant0,b1[1],i+1), 
     learning_curve(teleport_meant0,b1[1],i+1), 
     learning_curve(assemble_meant0,b1[1],i+1)]
make_chassis_base(logger2,t,2)
num_base += 1

while True:
    # need to receive wheelsets from upstream
    if (logger2.clock<upstream_parts_time[num_received_up]) or (num_base<=num_received_up):
        # keep making bases until getting wheelsets from upstream
        t = [learning_curve(spawn_meant0,b1[1],num_base+1), 
             learning_curve(teleport_meant0,b1[1],num_base+1), 
             learning_curve(assemble_meant0,b1[1],num_base+1)]
        make_chassis_base(logger2,t,2)
        num_base += 1
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[1],num_complete+1), 
             learning_curve(teleport_meant0,b1[1],num_complete+1), 
             learning_curve(assemble_meant0,b1[1],num_complete+1)]
        make_front(logger2,t,2)
        finish_station2(logger2,t,2)
        logger2.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[1],num_complete+1),2))
        num_complete += 1
    # completed all necessary chassis bases
    if num_base == num_order:
        break

# empty out chassis bases in inventory
while num_received_up < num_order:
    if logger2.clock < upstream_parts_time[num_received_up]:
        # wating wheelsets from upstream and complete chassis
        logger2.clock = upstream_parts_time[num_received_up]
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[1],num_complete+1), 
             learning_curve(teleport_meant0,b1[1],num_complete+1), 
             learning_curve(assemble_meant0,b1[1],num_complete+1)]
        make_front(logger2,t,2)        
        finish_station2(logger2,t,2)
        logger2.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[1],num_complete+1),2))
        num_complete += 1

# Station 3 - Windshield & steering
upstream_parts_idx = np.where(np.array(logger2.alog)=="Pass")[0]
upstream_parts_time = np.array(logger2.tlog)[upstream_parts_idx]

logger3 = ActionLogger()
num_received_up = 0
num_steer = 0
num_complete = 0

# first steering
t = [learning_curve(spawn_meant0,b1[2],num_steer+1), 
     learning_curve(teleport_meant0,b1[2],num_steer+1), 
     learning_curve(assemble_meant0,b1[2],num_steer+1)]
make_steering(logger3,t,3)
num_steer += 1

while True:
    # need to receive chassis from upstream
    if (logger3.clock<upstream_parts_time[num_received_up]) or (num_steer<=num_received_up):
        # keep making front body frame until getting chassis from upstream
        t = [learning_curve(spawn_meant0,b1[2],num_steer+1), 
             learning_curve(teleport_meant0,b1[2],num_steer+1), 
             learning_curve(assemble_meant0,b1[2],num_steer+1)]
        make_steering(logger3,t,3)
        num_steer += 1
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[2],num_complete+1), 
             learning_curve(teleport_meant0,b1[2],num_complete+1), 
             learning_curve(assemble_meant0,b1[2],num_complete+1)]
        finish_station3(logger3,t,3)
        logger3.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[2],num_complete+1),3))
        num_complete += 1

    # completed all necessary front body frames
    if num_steer == num_order:
        break

# empty out steering modules in inventory
while num_received_up < num_order:
    if logger3.clock < upstream_parts_time[num_received_up]:
        # waiting chassis coming from upstream
        logger3.clock = upstream_parts_time[num_received_up]
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[2],num_complete+1), 
             learning_curve(teleport_meant0,b1[2],num_complete+1), 
             learning_curve(assemble_meant0,b1[2],num_complete+1)]
        finish_station3(logger3,t,3)
        logger3.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[2],num_complete+1),3))
        num_complete += 1

# Station 4 -  Body & ceiling
upstream_parts_idx = np.where(np.array(logger3.alog)=="Pass")[0]
upstream_parts_time = np.array(logger3.tlog)[upstream_parts_idx]

logger4 = ActionLogger()
num_received_up = 0
num_body = 0
num_complete = 0

# first body frame
t = [learning_curve(spawn_meant0,b1[2],num_body+1), 
     learning_curve(teleport_meant0,b1[2],num_body+1), 
     learning_curve(assemble_meant0,b1[2],num_body+1)]
make_body_frame(logger4,t,4)
num_body += 1

while True:
    if (logger4.clock<upstream_parts_time[num_received_up]) or (num_body<=num_received_up):
        # keep making body frames until getting parts from upstream
        t = [learning_curve(spawn_meant0,b1[2],num_body+1), 
             learning_curve(teleport_meant0,b1[2],num_body+1), 
             learning_curve(assemble_meant0,b1[2],num_body+1)]
        make_body_frame(logger4,t,4)
        num_body += 1
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[2],num_complete+1), 
             learning_curve(teleport_meant0,b1[2],num_complete+1), 
             learning_curve(assemble_meant0,b1[2],num_complete+1)]
        finish_station4(logger4,t,4)
        logger4.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[3],num_complete+1),4))
        num_complete += 1

    # completed all necessary body frames
    if num_body == num_order:
        break

# empty out body frames in inventory
while num_received_up < num_order:
    if logger4.clock < upstream_parts_time[num_received_up]:
        # waiting parts coming from upstream
        logger4.clock = upstream_parts_time[num_received_up]
    else:
        num_received_up += 1
        t = [learning_curve(spawn_meant0,b1[2],num_complete+1), 
             learning_curve(teleport_meant0,b1[2],num_complete+1), 
             learning_curve(assemble_meant0,b1[2],num_complete+1)]
        finish_station4(logger4,t,4)
        logger4.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b1[3],num_complete+1),4))
        num_complete += 1

# player,action,time_stamp
log_df1 = pd.DataFrame({'player':"player1", 'action':logger1.alog, 'time_stamp':logger1.tlog})
log_df2 = pd.DataFrame({'player':"player2", 'action':logger2.alog, 'time_stamp':logger2.tlog})
log_df3 = pd.DataFrame({'player':"player3", 'action':logger3.alog, 'time_stamp':logger3.tlog})
log_df4 = pd.DataFrame({'player':"player4", 'action':logger4.alog, 'time_stamp':logger4.tlog})
log_df = pd.concat([log_df1,log_df2,log_df3,log_df4])
log_df = log_df.sort_values(by=['time_stamp'],ascending=True)
log_df.to_csv(save_dir/('mass'+str(num_order)+'.csv'),index=False)

#%% Craft production
logger1 = ActionLogger()
logger2 = ActionLogger()
logger3 = ActionLogger()
logger4 = ActionLogger()

for i in range(int(num_order/4)):
    t = [learning_curve(spawn_meant0,b2[0],i+1), 
         learning_curve(teleport_meant0,b2[0],i+1), 
         learning_curve(assemble_meant0,b2[0],i+1)]
    craft_car(logger1,t,5)
    logger1.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b2[0],i+1),5))

    # player 2
    t = [learning_curve(spawn_meant0,b2[1],i+1), 
         learning_curve(teleport_meant0,b2[1],i+1), 
         learning_curve(assemble_meant0,b2[1],i+1)]
    craft_car(logger2,t,6)
    logger2.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b2[1],i+1),6))

    # player 3
    t = [learning_curve(spawn_meant0,b2[2],i+1), 
         learning_curve(teleport_meant0,b2[2],i+1), 
         learning_curve(assemble_meant0,b2[2],i+1)]    
    craft_car(logger3,t,7)
    logger3.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b2[2],i+1),7))

    # player 4
    t = [learning_curve(spawn_meant0,b2[3],i+1), 
         learning_curve(teleport_meant0,b2[3],i+1), 
         learning_curve(assemble_meant0,b2[3],i+1)]    
    craft_car(logger4,t,8)
    logger4.write_act("Pass",RNG.Expon(learning_curve(pass_meant0,b2[3],i+1),8))
    
log_df1 = pd.DataFrame({'player':"player1", 'action':logger1.alog, 'time_stamp':logger1.tlog})
log_df2 = pd.DataFrame({'player':"player2", 'action':logger2.alog, 'time_stamp':logger2.tlog})
log_df3 = pd.DataFrame({'player':"player3", 'action':logger3.alog, 'time_stamp':logger3.tlog})
log_df4 = pd.DataFrame({'player':"player4", 'action':logger4.alog, 'time_stamp':logger4.tlog})
log_df = pd.concat([log_df1,log_df2,log_df3,log_df4])
log_df = log_df.sort_values(by=['time_stamp'],ascending=True)
log_df.to_csv(save_dir/('craft'+str(num_order)+'.csv'),index=False)
    