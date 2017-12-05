import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.mlab as mlab
import matplotlib.ticker as mtick
import itertools
import math
import random
import networkx as nx
import shapefile as shp
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import sys
import time

from scipy.spatial import Delaunay
from scipy.interpolate import interp1d 

def cycle_counter(SOC):
    
    keep = SOC[0]
    
    cycle_count = 0
    
    for i in range(1, len(SOC)):
        
        if SOC[i] < SOC[i-1]:
            #print 'i = ', i
            cycle_count = (SOC[i-1]-keep) + cycle_count
            keep = SOC[i]
            #print 'cycle ', cycle_count
            #print 'keep ', keep
    return cycle_count

#------- calculates distance in km between two coordinates (input must be in degrees) --------
def geoDist(c1, c2):
    R = 6371
    lon1, lat1 = math.radians(c1[0]), math.radians(c1[1])
    lon2, lat2 = math.radians(c2[0]), math.radians(c2[1])
    a = math.sin((lat1-lat2)/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin((lon1-lon2)/2)**2
    d = 2*R*math.asin(math.sqrt(a))

    return d

def makeTree(points, clusterHomes, parcCoordDict):
    tri = Delaunay(points)
    find_neighbors = lambda x, tri: list(set(vertex for triang in tri.vertices if x in triang for vertex \
    in triang if vertex !=x))
    adj_list = [None] * len(clusterHomes)
    for n in range(len(clusterHomes)):
        adj_list[n] = find_neighbors(n, tri)
        adj_list[n].sort()
        adj_list[n] = [clusterHomes[i] for i in adj_list[n]] # get the id of the clusterhome
        # format adj_list in NetworkX adjlist format:
        adj_list[n].insert(0, clusterHomes[n]) # add in the particular home as the first list element
        adj_list[n] = str(adj_list[n]).replace(',', '')
        adj_list[n] = str(adj_list[n]).strip('[]')

    G = nx.read_adjlist(adj_list)
    #print G.nodes
    nodeList = np.zeros(len(clusterHomes), dtype=int)
    for i in range(len(clusterHomes)):
        n = G.nodes()[i].replace('u', '')
        #print G.nodes()[i]
        nodeList[i] = n
    nodeList = nodeList.tolist()
    nodeDict = { key: value for (key, value) in zip(G.nodes(), nodeList) }
    G = nx.relabel_nodes(G, nodeDict)

    edge_dist = zip(G.edges(), [geoDist(parcCoordDict.get(edge[0]), parcCoordDict.get(edge[1])) for edge in G.edges()])
    distDict = { key: value for (key, value) in edge_dist }
    nx.set_edge_attributes(G, 'distances', distDict)
    G = nx.minimum_spanning_tree(G, 'distances')
    
    return G


def cycle_counter(SOC):
    
    keep = SOC[0]
    
    cycle_count = 0
    
    for i in range(1, len(SOC)):
        
        if SOC[i] < SOC[i-1]:
            #print 'i = ', i
            cycle_count = (SOC[i-1]-keep) + cycle_count
            keep = SOC[i]
            #print 'cycle ', cycle_count
            #print 'keep ', keep
    return cycle_count

def calculate_IRR(capCost, yearlyIncome, lifetime, yearlyCosts, i_if, c_if, r):
    #print capCost,  yearlyIncome, lifetime, yearlyCosts, i_if, c_if, r
    cumulativeCosts = float( -capCost )
    NPV = []
    NPV.append( float( -capCost ) )
    Costs_list = []
    Costs_list.append( cumulativeCosts )
    for i in np.arange(1,lifetime):
        # at the end of year i
        cumulativeCosts = yearlyIncome*(1+i_if)**i - yearlyCosts*(1+c_if)**i
        Costs_list.append( cumulativeCosts )
        NPV.append( cumulativeCosts/((1+r)**i) )
    cumulativeCosts = yearlyIncome*(1+i_if)**i - yearlyCosts*(1+c_if)**i
    Costs_list.append( cumulativeCosts*(lifetime-i) )
    NPV.append( (cumulativeCosts*(lifetime-i)) / ((1+r)**lifetime) )
        
    #print Costs_list
    b = np.sum(NPV)
    a = np.irr( np.asarray( Costs_list ) )    
        
    return a, b

def cycleCharging(noUnits, storageProperties, load, price, billingOption, Delta_t):

    seriesLength=len(load)
    demand = -load
    original_demand = np.zeros((seriesLength))
    original_demand[ np.where(load<0)[0] ] = -load[ np.where(load<0)[0] ]
    
    ### ------------------------- ESS properties ----------------- ###
    maxSOC = noUnits*storageProperties[0]/Delta_t# kWh
    # scale maxSOC appropriately
    maxChg = noUnits*storageProperties[1] # kW
    maxDisChg = noUnits*storageProperties[2] # kW
    etaChg = storageProperties[3]
    etaDisChg = storageProperties[4]
    ### ---------------------------------------------------------- ###

    ### ------------------------- CHOOSE POLICY OPTION ----------- ###
    policyOption = billingOption[0] 
    marketArbitrage = billingOption[2]
    ### ---------------------------------------------------------- ###
    
    priceSolarExport = billingOption[1]
    solarExports = np.zeros((seriesLength))
    exports = np.zeros((seriesLength))
    buy_price = np.zeros((seriesLength))
    for i in range(seriesLength):
        buy_price[i] = price[i]
        if policyOption == 1:
            if demand[i] < 0:
                solarExports[i]=-demand[i]
                exports[i] = solarExports[i]
                demand[i] = 0
    # now alter the price in the periods where there is solar available (if req.)
    buy_price[ np.where(solarExports>0)[0] ] = priceSolarExport
    
    # get the storage profiles
    SOC = np.zeros((seriesLength))
    deltaSOC = np.zeros((seriesLength))

    # calculate initial costs
    oldPrice = np.zeros((seriesLength))
    for i in range(seriesLength):
        oldPrice[i]=(price[i]*demand[i] - solarExports[i]*buy_price[i])*Delta_t
    oldTotalPrice = np.sum( oldPrice )

    # cycle through the series
    for i in range(seriesLength):
        bottleneck = np.zeros((3))
        if exports[i]>0:
            # establish the min bottleneck for charging
            bottleneck[0] = maxSOC - SOC[i-1]
            bottleneck[1] = maxChg
            bottleneck[2] = exports[i]*etaChg
            # act on the min bottleneck
            actual_bottleneck = np.min(bottleneck) 
            exports[i] = exports[i] - actual_bottleneck/etaChg
            SOC[i] = SOC[i-1] + actual_bottleneck
            deltaSOC[i] = 0 + actual_bottleneck
        else:
            # establish the min bottleneck for discharging
            bottleneck[0] = SOC[i-1]
            bottleneck[1] = -maxDisChg
            bottleneck[2] = demand[i]*etaDisChg
            # operate on teh actual bottleneck
            actual_bottleneck = np.min(bottleneck)
            demand[i] = demand[i] - actual_bottleneck*etaDisChg
            SOC[i] = SOC[i-1] - actual_bottleneck
            deltaSOC[i] = 0 - actual_bottleneck
            
    newPrice = np.zeros((seriesLength))
    for i in range(seriesLength):
        newPrice[i] = (price[i]*demand[i] - exports[i]*buy_price[i])*Delta_t
    newTotalPrice = np.sum ( newPrice )
    #print oldTotalPrice - newTotalPrice
    
    #now calculate the new demand
    newDemand = demand-exports
    SOC = SOC*Delta_t
    
    return newDemand, oldTotalPrice, newTotalPrice, SOC

def ESSRGModule(noUnits, storageProperties, load, price, billingOption, Delta_t):

    seriesLength=len(load)
    demand = -load
    original_demand = np.zeros((seriesLength))
    original_demand[ np.where(load<0)[0] ] = -load[ np.where(load<0)[0] ]
    
    ### ------------------------- ESS properties ----------------- ###
    maxSOC = noUnits*storageProperties[0]/Delta_t# kWh
    # scale maxSOC appropriately
    maxChg = noUnits*storageProperties[1] # kW
    maxDisChg = noUnits*storageProperties[2] # kW
    etaChg = storageProperties[3]
    etaDisChg = storageProperties[4]
    ### ---------------------------------------------------------- ###

    ### ------------------------- CHOOSE POLICY OPTION ----------- ###
    policyOption = billingOption[0] 
    marketArbitrage = billingOption[2]
    ### ---------------------------------------------------------- ###
    
    priceSolarExport = billingOption[1]
    solarExports = np.zeros((seriesLength))
    exports = np.zeros((seriesLength))
    buy_price = np.zeros((seriesLength))
    for i in range(seriesLength):
        buy_price[i] = price[i]
        if policyOption == 1:
            if demand[i] < 0:
                solarExports[i]=-demand[i]
                exports[i] = solarExports[i]
                demand[i] = 0
    # now alter the price in the periods where there is solar available (if req.)
    buy_price[ np.where(solarExports>0)[0] ] = priceSolarExport
    
    # get the storage profiles
    SOC = np.zeros((seriesLength))
    deltaSOC = np.zeros((seriesLength))

    # calculate initial costs
    oldPrice = np.zeros((seriesLength))
    for i in range(seriesLength):
        oldPrice[i]=(price[i]*demand[i] - solarExports[i]*buy_price[i])*Delta_t
    oldTotalPrice = np.sum( oldPrice )

    remove = np.ones((seriesLength))

    while np.any(remove) == True:
        matrix = np.zeros(( 3, len ( np.where(remove==True)[0] ) ))
        matrix[0,:] = np.where(remove==True)[0]
        matrix[1,:] = price[ np.where(remove==True)[0] ]
        matrix[2,:] = buy_price[ np.where(remove==True)[0] ]
        indici = np.where( matrix[1,:]==np.max( matrix[1,:] ) )[0]

        # try and preferentially select maxh according to minh
        priceMatch = np.where( matrix[1,:]==matrix[2,:] )[0]
        indPref = []
        for x in indici:
            if x in priceMatch:
                indPref.append(x)
        if not(indPref):
            maxh = matrix[0,indici[0]]
        else:
            maxh = matrix[0,indPref[0]]
            
        maxh = maxh.astype(int)
        
        # find the last hour before maxh when storage was full
        r1 = np.where( SOC[0:maxh+1] == maxSOC )[0]
        if r1.size == 0:
            r1 = 0
        else:
            r1 = np.where( SOC[0:maxh+1] == maxSOC )[0][-1]+1
        # find the first hour after maxh when storage is empty
        r2 = np.where( SOC[maxh:] == 0 )[0]
        if r2.size == 0:
            r2 = len(SOC)
        else:
            r2 = np.where( SOC[maxh:] == 0 )[0][0]+maxh-1
            
        # find the minh in the time range
        range_price = buy_price[r1:r2+1]
        range_remove = remove[r1:r2+1]
        
        # if there is no hour in the range then remove maxh and skip to the end
        if np.any(range_remove) == False:
            remove[maxh]=False
        else:
            matrix = np.zeros(( 2, len ( np.where(range_remove==True)[0] ) ))
            matrix[0,:] = np.where(range_remove==True)[0]
            matrix[1,:] = range_price[ np.where(range_remove==True)[0] ]
            #print matrix
            indici = np.where( matrix[1,:]==np.min( matrix[1,:] ) )[0]
            minh = matrix[0,indici[0]] + r1
            
            minh = minh.astype(int)
        
            MoC = buy_price[minh]/(etaChg*etaDisChg)
            #print 'Minh = ', minh
        
            if MoC<price[maxh] and [minh!=maxh]:
            
                bottleneck = np.zeros((5))
                bottleneck[0] = deltaSOC[maxh]-maxDisChg
                bottleneck[1] = maxChg - deltaSOC[minh]
                if maxh > minh:
                    bottleneck[2] = maxSOC - np.max(SOC[minh:maxh+1])
                else:
                    bottleneck[2] = np.min(SOC[maxh:minh+1])
                    
                ############### ARBITRAGE ###########################    
                # initially cannot output more than local consumption
                if marketArbitrage == False:
                    if demand[maxh]>0:
                        bottleneck[3] = demand[maxh]/etaDisChg
                    else:
                        bottleneck[3] = 0
                else:
                    bottleneck[3] = maxSOC
                    
                # if there is export initally cannot store more than export
                if exports[minh] > 0:
                    bottleneck[4] = exports[minh]*etaChg
                else:
                    bottleneck[4] = maxSOC
                #print bottleneck
                actual_bottleneck = np.min(bottleneck)  
                
                # Update the aggregated demand
                demand[maxh] = demand[maxh] - actual_bottleneck*etaDisChg
                # if not buying from solar increase the demand at minh
                if exports[minh] <= 0:
                    demand[minh] = demand[minh] + actual_bottleneck/etaChg
                
                # update the storage charging schedule
                exports[minh] = exports[minh] - actual_bottleneck/etaChg
                deltaSOC[minh] = deltaSOC[minh] + actual_bottleneck
                deltaSOC[maxh] = deltaSOC[maxh] - actual_bottleneck
                if maxh>minh:
                    SOC[minh:maxh] = SOC[minh:maxh] + actual_bottleneck
                else:
                    SOC[maxh:minh] = SOC[maxh:minh] - actual_bottleneck
            
                # account for rounding errors
                if demand[maxh] < 0 + 0.0001 and demand[maxh] > 0 - 0.0001:
                    demand[maxh] = 0
                if exports[minh] < 0 + 0.0001:
                    exports[minh] = 0
                    # update the price
                    buy_price[minh] = price[minh]
                
                #print SOC
                #print deltaSOC
            
                # check if at the charge or discharge operation is at capacity at either
                # maxh or minh and remove that hour from the price distribution.
                if deltaSOC[maxh] <= maxDisChg:
                    remove[maxh] = False
                if deltaSOC[minh] >= maxChg:
                    remove[minh] = False
                
                ######################## ARBITRAGE #####################
                # check whether exports are allowed or not
                if marketArbitrage == False:
                    # if not remove periods where all demand is met
                    if demand[maxh] <= 0:
                        remove[maxh] = False
                # else: do nothing
            
            else:
                # if there is no price incentive remove the hours
                remove[maxh] = False
                remove[minh] = False
            
    newPrice = np.zeros((seriesLength))
    for i in range(seriesLength):
        newPrice[i] = (price[i]*demand[i] - exports[i]*buy_price[i])*Delta_t
    newTotalPrice = np.sum ( newPrice )
    #print oldTotalPrice - newTotalPrice
    
    #print 'state of charge'
    #print SOC
    #print 'delta SOC'
    #print np.sum(deltaSOC[deltaSOC>0]), np.sum(deltaSOC[deltaSOC<0])
    #print 'old demand'
    #print original_demand
    
    #now calculate the new demand
    newDemand = demand-exports
    SOC = SOC*Delta_t
    
    return newDemand, oldTotalPrice, newTotalPrice, SOC