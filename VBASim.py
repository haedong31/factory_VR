# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 18:08:03 2016
@author: yl17
Updated to R2 by bln Oct 4 2017
Updated to R3 by bln Oct 11 2017

"""

import BasicClasses

def VBASimInit(calendar,queues,ctstats,dtstats,resources,clock):
    # Function to initialize VBASim.Python
    # Typically called before the first replication and between replications
    
    #Empty event calendar
    while (calendar.N() > 0):
        EV = calendar.Remove()
        
    #Empty queues
    #On first call, append the CStats created by FIFOQueue
    for Q in queues:
        if Q.WIP not in ctstats:
            ctstats.append(Q.WIP)
        while Q.NumQueue() > 0:
            En = Q.Remove(clock)

    #Reinitialize Resources
    #On first call, append the CStats created by FIFOQueue and Resource
    for Re in resources:
        Re.Busy = 0.0
        if Re.NumBusy not in ctstats:
            ctstats.append(Re.NumBusy)  
    
    #Clear statistics
    for CT in ctstats:
        CT.Clear(clock)
        CT.Xlast = 0.0   # added bln 
        
    for DT in dtstats:
        DT.Clear()
    



def Schedule(calendar,EventType, EventTime, clock):
    #Schedule future events of EventType to occur at time Clock + EventTime
    
    addedEvent = BasicClasses.EventNotice()
    addedEvent.EventType = EventType
    addedEvent.EventTime = clock + EventTime
    calendar.Schedule(addedEvent)
    

    
def SchedulePlus(calendar,EventType, EventTime, TheObject, clock):
    #Schedule future events of EventType to occur at time Clock + EventTime
    #and pass with the event notice TheObject
    
    addedEvent = BasicClasses.EventNotice()
    addedEvent.EventType = EventType
    addedEvent.EventTime = clock + EventTime
    addedEvent.WhichObject = TheObject
    calendar.Schedule(addedEvent)
    
    
def ClearStats(ctstats,dtstats,clock):
    #Clear statistics in TheDTStats and TheCTStats
    
    for CT in ctstats:
        CT.Clear(clock)
    for DT in dtstats:
        DT.Clear()