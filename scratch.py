import pandas as pd
import numpy as np
import RNG

num_order = np.random.randint(low=5,high=11,size=1)[0]#  ~ Unif[1,11)

# unit of time: minuate
spawn_meant = 0.25
teleport_meant = 0.5
assemble_meant = 1.0
insp_meant = 1.5
pass_meant = 1.0

class ActionLogger(object):
    def __init__(self):
        self.clock = 0.0
        self.alog = []
        self.tlog = []

    def write_act(self,a,t):
        self.clock += t
        self.alog.append(a)
        self.tlog.append(self.clock)

logger1 = ActionLogger()
for _ in range(num_order):
    for _ in range(2):
        # sourcing
        for _ in range(7):
            logger1.write_act("SpawnPart",RNG.Expon(spawn_meant,1))

        # manufacturing
        # (rim + tire) x 2
        for _ in range(2):
            logger1.write_act("TeleportPart",RNG.Expon(teleport_meant,1)) # rim
            logger1.write_act("TeleportPart",RNG.Expon(teleport_meant,1)) # tire
            logger1.write_act("Assemble",RNG.Expon(assemble_meant,1)) # rim + tire

        # (wheels + axle)
        logger1.write_act("TeleportPart",RNG.Expon(teleport_meant,1)) # axle
        logger1.write_act("Assemble",RNG.Expon(assemble_meant,1)) # axle + wheel 1
        logger1.write_act("Assemble",RNG.Expon(assemble_meant,1)) # axle + wheel 2
        
        # (axle + axle)
        logger1.write_act("TeleportPart",RNG.Expon(teleport_meant,1)) # connecting plate
        logger1.write_act("Assemble",RNG.Expon(assemble_meant,1)) # axle + axle

        # inspection
        logger1.write_act("Inspection",RNG.Expon(insp_meant,1))
    # pass wheelsets to downstream
    logger1.write_act("Pass",RNG.Expon(pass_meant,1))

    upstream_parts_idx = np.where(np.array(logger1.alog)=="Pass")[0]
upstream_parts_time = np.array(logger1.tlog)[upstream_parts_idx]

logger2 = ActionLogger()
num_base_inv = 0
num_front_inv = 0
num_received_wheelset = 0
num_base = 0
num_front = 0
num_complete = 0

# first chassis base (long plate + long plate)
for _ in range(4):
    logger2.write_act("SpawnPart",RNG.Expon(spawn_meant,2))
logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
num_base_inv += 1
num_base += 1

while True:
    # need to receive wheelsets from the upstream
    if logger2.clock < upstream_parts_time[num_received_wheelset]:
        # keep making bases until getting wheelsets from upstream
        for _ in range(4):
            logger2.write_act("SpawnPart",RNG.Expon(spawn_meant,2))        
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        num_base += 1
    else:
        num_received_wheelset += 1
        for _ in range(5):
            logger2.write_act("SpawnPart",RNG.Expon(spawn_meant,2))

        # (wheelsets + base)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))        
        
        # (front bumper)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # (engine bay)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # assemble base and front module
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # inspection & pass chassis to downstream
        logger2.write_act("Inspection",RNG.Expon(insp_meant,2))
        logger2.write_act("Pass",RNG.Expon(pass_meant,2))
        num_complete += 1

    # completed all necessary chassis bases
    if num_base == num_order:
        break

# empty out chassis bases in inventory
while num_received_wheelset < num_order:
    if (logger2.clock < upstream_parts_time[num_received_wheelset]):
        # wating wheelsets from upstream and complete chassis
        logger2.clock = upstream_parts_time[num_received_wheelset]
    else:
        num_received_wheelset += 1
        for _ in range(5):
            logger2.write_act("SpawnPart",RNG.Expon(spawn_meant,2))

        # (wheelsets + base)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))        
        
        # (front bumper)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # (engine bay)
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))
        logger2.write_act("TeleportPart",RNG.Expon(teleport_meant,2))
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # assemble base and front module
        logger2.write_act("Assemble",RNG.Expon(assemble_meant,2))

        # inspection & pass chassis to downstream
        logger2.write_act("Inspection",RNG.Expon(insp_meant,2))
        logger2.write_act("Pass",RNG.Expon(pass_meant,2))
        num_complete += 1