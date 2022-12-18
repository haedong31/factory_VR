import pandas as pd
import numpy as np
import RNG
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# num_order = np.random.randint(low=5,high=11,size=1)[0] # ~ Unif[1,11)
num_order = 8

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

def exp_decay(x,a,b):
    return a*np.exp(-b*x)

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

    # inspection
    # logger.write_act("Inspection",RNG.Expon(insp_meant,rnd_seed))

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

    # inspection
    # logger.write_act("Inspection",RNG.Expon(insp_meant,rnd_seed))

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

    # inspection & pass to downstream
    # logger.write_act("Inspection",RNG.Expon(insp_meant,rnd_seed))

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

    # inspection & pass to downstream
    # logger.write_act("Inspection",RNG.Expon(insp_meant,rnd_seed))

def craft_car(logger,rnd_seed):
    make_wheelset(logger,rnd_seed)
    make_wheelset(logger,rnd_seed)
    make_chassis_base(logger,rnd_seed)
    make_front(logger,rnd_seed)
    finish_station2(logger,rnd_seed)
    make_steering(logger,rnd_seed)
    finish_station3(logger,rnd_seed)
    make_body_frame(logger,rnd_seed)
    finish_station4(logger,rnd_seed)
    
# unit of time: minuate
spawn_meant0 = 0.5
teleport_meant0 = 0.5
assemble_meant0 = 1.0
pass_meant0 = 0.25
# insp_meant = 1.5

spawn_meant_param,_ = curve_fit(exp_decay,np.array([0,num_order-1]),np.array([spawn_meant0,spawn_meant0*0.5]))
teleport_meant_param,_ = curve_fit(exp_decay,np.array([0,num_order-1]),np.array([teleport_meant0,teleport_meant0*0.5]))
assemble_meant_param,_ = curve_fit(exp_decay,np.array([0,num_order-1]),np.array([assemble_meant0,assemble_meant0*0.5]))
pass_meant_param,_ = curve_fit(exp_decay,np.array([0,num_order-1]),np.array([pass_meant0,pass_meant0*0.5]))

x = np.arange(num_order*2)
spawn_meant = exp_decay(x,spawn_meant_param[0],spawn_meant_param[1])
teleport_meant = exp_decay(x,teleport_meant_param[0],teleport_meant_param[1])
assemble_meant = exp_decay(x,assemble_meant_param[0],assemble_meant_param[1])
pass_meant = exp_decay(x,pass_meant_param[0],pass_meant_param[1])

logger1 = ActionLogger()
for i in range(num_order):
    t = [spawn_meant[2*i],teleport_meant[2*i],assemble_meant[2*i]]
    make_wheelset(logger1,t,1)
    t = [spawn_meant[2*i+1],teleport_meant[2*i+1],assemble_meant[2*i+1]]
    make_wheelset(logger1,t,1)
    # pass wheelsets to downstream
    logger1.write_act("Pass",RNG.Expon(pass_meant[i],1))
    
x = np.arange(num_order)
spawn_meant = exp_decay(x,spawn_meant_param[0],spawn_meant_param[1])
teleport_meant = exp_decay(x,teleport_meant_param[0],teleport_meant_param[1])
assemble_meant = exp_decay(x,assemble_meant_param[0],assemble_meant_param[1])
pass_meant = exp_decay(x,pass_meant_param[0],pass_meant_param[1])

upstream_parts_idx = np.where(np.array(logger1.alog)=="Pass")[0]
upstream_parts_time = np.array(logger1.tlog)[upstream_parts_idx]

logger2 = ActionLogger()
num_received_up = 0
num_base = 0
num_complete = 0

# first chassis base
t = [spawn_meant[num_base],teleport_meant[num_base],assemble_meant[num_base]]
make_chassis_base(logger2,t,1)
num_base += 1

while True:
    # need to receive wheelsets from upstream
    if (logger2.clock<upstream_parts_time[num_received_up]) or (num_base<=num_received_up):
        # keep making bases until getting wheelsets from upstream
        t = [spawn_meant[num_base],teleport_meant[num_base],assemble_meant[num_base]]
        make_chassis_base(logger2,t,1)
        num_base += 1
    else:
        num_received_up += 1
        t = [spawn_meant[num_complete],teleport_meant[num_complete],assemble_meant[num_complete]]
        make_front(logger2,t,1)
        finish_station2(logger2,t,1)
        logger2.write_act("Pass",RNG.Expon(pass_meant[num_complete],1))
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
        t = [spawn_meant[num_complete],teleport_meant[num_complete],assemble_meant[num_complete]]
        make_front(logger2,t,1)        
        finish_station2(logger2,t,1)
        logger2.write_act("Pass",RNG.Expon(pass_meant[num_complete],2))
        num_complete += 1