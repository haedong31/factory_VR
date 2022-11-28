# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 21:41:39 2016
@author: yl17
Updated to R2 bln Oct 4 2017
"""
import math

class Activity():
    def __init__(self):
        self.WhichActivity = 0
        self.WhichNode = 0

    
class CTStat():
    def __init__(self):
        self.Area = 0.0
        self.Tlast = 0.0
        self.TClear = 0.0
        self.Xlast = 0.0
        
    def Record(self,X,clock):
        self.Area = self.Area + self.Xlast * (clock - self.Tlast)
        self.Tlast = clock
        self.Xlast = X

    def Mean(self,clock):
        mean = 0.0
        if (clock - self.TClear) > 0.0:
           mean = (self.Area + self.Xlast * (clock - self.Tlast)) / (clock - self.TClear)
        return mean
    
    def Clear(self,clock):
        self.Area = 0.0
        self.Tlast = clock
        self.TClear = clock


class DTStat():
    def __init__(self):
        self.Sum = 0.0
        self.SumSquared = 0.0
        self.NumberOfObservations = 0.0
    
    def Record(self,X):
        self.Sum = self.Sum + X
        self.SumSquared = self.SumSquared + X * X
        self.NumberOfObservations = self.NumberOfObservations + 1
        
    def Mean(self):
        mean = 0.0
        if self.NumberOfObservations > 0.0:
            mean = self.Sum / self.NumberOfObservations
        return mean

    def StdDev(self):
        stddev = 0.0
        if self.NumberOfObservations > 1.0:
            stddev = math.sqrt((self.SumSquared - self.Sum**2 / self.NumberOfObservations) / (self.NumberOfObservations - 1))
        return stddev
            
    def N(self):
        return self.NumberOfObservations
    
    def Clear(self):
        self.Sum = 0.0
        self.SumSquared = 0.0
        self.NumberOfObservations = 0.0


class Entity():
    def __init__(self,clock,Type):
        self.CreateTime = clock
        self.Type = Type
        self.Loc = [-1.0,-1.0]
        self.aType = -1
        self.Station = -1
            
        
class InventoryEntity():
    def __init__(self,clock,currentstock):
        self.CreateTime = clock
        self.CurrentStock = currentstock
        
        
class EventNotice():
    def __init__(self):
        self.EventTime = 0.0
        self.EventType = ""
        self.WhichObject = ()
        
        
class EventCalendar():
    def __init__(self):
        self.ThisCalendar = []
    
    def Schedule(self,addedEvent):
        if len(self.ThisCalendar) == 0:
            self.ThisCalendar.append(addedEvent)
        elif self.ThisCalendar[-1].EventTime <= addedEvent.EventTime:
            self.ThisCalendar.append(addedEvent)
        else:
            for rep in range(0,len(self.ThisCalendar),1):
                if self.ThisCalendar[rep].EventTime > addedEvent.EventTime:
                    break
            self.ThisCalendar.insert(rep,addedEvent)
    
    def Remove(self):
        if len(self.ThisCalendar) > 0:
            remove = self.ThisCalendar[0]
            self.ThisCalendar.remove(self.ThisCalendar[0])
            return remove
        
    def N(self):
        return len(self.ThisCalendar)
        
        
    
class FIFOQueue():
    def __init__(self):
        self.WIP = CTStat()
        self.ThisQueue = []
        
    def NumQueue(self):
        return len(self.ThisQueue)
        
    def Add(self,X,clock):
        self.ThisQueue.append(X)
        numqueue = self.NumQueue()
        self.WIP.Record(float(numqueue),clock)    
    
    def Remove(self,clock):
        if len(self.ThisQueue) > 0:
            remove = self.ThisQueue[0]
            self.ThisQueue.remove(self.ThisQueue[0])
            self.WIP.Record(float(self.NumQueue()),clock)   
            return remove
        
    def Mean(self,clock):
        return self.WIP.Mean(clock)
        
    def Clear(self,clock):
        self.WIP.Clear(clock)
        
class Resource():
    def __init__(self):
        self.Busy = 0
        self.NumberOfUnits = 0
        self.NumBusy = CTStat()
        
    def Seize(self, Units, clock):
        diff = self.NumberOfUnits - Units - self.Busy
        if diff >= 0:
            self.Busy = self.Busy + Units
            self.NumBusy.Record(float(self.Busy),clock)
            seize = True
        else:
            seize = False
        return seize
        
    def Free(self, Units, clock):
        diff = self.Busy - Units
        if diff < 0:
            free = False
        else:
            self.Busy = self.Busy - Units
            self.NumBusy.Record(float(self.Busy), clock)
            free = True
        return free
    
    def Mean(self,clock):
        return self.NumBusy.Mean(clock)
        
    def SetUnits(self, Units):
        self.NumberOfUnits = Units
        
    def Clear(self,clock):
        self.NumBusy.Clear(clock)
        
        
class PriorityQueue():
    def __init__(self):
        self.WIP = CTStat()
        self.ThisQueue = []
        
    def NumQueue(self):
        return len(self.ThisQueue)
        
    def Add(self,X,clock):
        if len(self.ThisQueue) == 0:
            self.ThisQueue.append(X)
        else:
            if self.ThisQueue[-1].Priority >= X.Priority:
                self.ThisQueue.append(X)
            else:
                for rep in range(0,len(self.ThisQueue),1):
                    if self.ThisQueue[rep].Priority < X.Priority:
                        break
                self.ThisQueue.insert(rep,X)
        numqueue = self.NumQueue()
        self.WIP.Record(float(numqueue),clock)    
    
    def Remove(self,clock):
        if len(self.ThisQueue) > 0:
            remove = self.ThisQueue[0]
            self.ThisQueue.remove(self.ThisQueue[0])
            return remove
        numqueue = self.NumQueue()
        self.WIP.Record(float(numqueue),clock)        
        
    def Mean(self,clock):
        return self.WIP.Mean(clock)
        
