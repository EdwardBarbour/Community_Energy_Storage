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
from scipy.spatial import Delaunay

def DCpowflow(N, nxgraph): 

    #print 'Calculating DC power flows...'

    loadDict = nx.get_node_attributes(nxgraph, 'loads')
    # print loadDict.items()

    # sort nodes and convert to numpy array
    loads = np.array([loadDict.get(i) for i in sorted(loadDict.keys())])
    reactDict = nx.get_edge_attributes(nxgraph, 'reactances')
    edges = reactDict.keys()

    # find nodes outside of the island to exclude	
    missing = np.setdiff1d(np.arange(N), loadDict.keys())

    # determine balancing power per node
    fromGrid = 0
    toGrid = 0
    if np.absolute(np.sum(loads[loads < 0])) > np.sum(loads[loads > 0]):
        fromGrid = fromGrid + -np.sum(loads[loads < 0])\
        * (1 - np.sum(loads[loads > 0])/np.absolute(np.sum(loads[loads < 0]))) #/ len(nxgraph.nodes())
    else:
        toGrid = toGrid + np.sum(loads[loads > 0])\
        * (1 - np.absolute(np.sum(loads[loads < 0]))/np.sum(loads[loads > 0])) #/ len(nxgraph.nodes())

    # balance load and generation
    # (if L > G, decrease load; if G > L, decrease generation)
    if np.absolute(np.sum(loads[loads < 0])) > np.sum(loads[loads > 0]):
        loads[loads < 0] = loads[loads < 0] / np.absolute(np.sum(loads[loads < 0]))*np.sum(loads[loads > 0])
    else:
        loads[loads > 0] = loads[loads > 0] / np.sum(loads[loads > 0])*np.absolute(np.sum(loads[loads < 0]))

    B = np.zeros((N, N))
    for edge in edges:
        B[edge[0], edge[1]] = -1/reactDict[edge]
        B[edge[1], edge[0]] = -1/reactDict[edge] 

    # delete empty rows and columns 	
    B = np.delete(np.delete(B, missing, 0), missing, 1)  

    reacsum = -np.sum(B, 0)
    
    for i in range(len(B)):
        B[i, i] = reacsum[i]

    theta = np.linalg.solve(B[1:, 1:], loads[1:])
    theta = np.insert(theta, 0, 0)

    nodeList = sorted(loadDict.keys())
    edgesNew = [(nodeList.index(edge[0]), nodeList.index(edge[1])) for edge in edges]
    
    flows = [-B[edge[0], edge[1]] * (theta[edge[0]]-theta[edge[1]]) for edge in edgesNew]
    flowDict = { key: value for (key, value) in zip(edges, flows) }

    return (flowDict, fromGrid, toGrid)
    
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
    
    
def get_grid_Power(G, nodeLabelDict, flowDict, PCC):
    
    end = nodeLabelDict.keys()[0]
    gridPower = np.zeros((12))
    gridP = -PCC
    for neigh in G.neighbors(end):
        try: 
            checkFlow = flowDict[end, neigh]
        except:
            checkFlow = -flowDict[neigh, end]
        if checkFlow > 0:
            gridP = gridP + checkFlow
    if gridP > 0:
        gridPower[0] = 0
    else:
        gridPower[0] = gridP

    for nodeNum in np.arange(1,12):

        start = nodeLabelDict.keys()[nodeNum]
        #print nodeLabelDict[start], nodeLabelDict[end]

        for path in nx.all_simple_paths(G, source=start, target=end):
            #print path
            testFlow = np.zeros((len(path)-1))
            for i in range(len(testFlow)):
                #print path[i]
                #print path[i+1]
                try: 
                    testFlow[i] = flowDict[path[i], path[i+1]]
                except:
                    testFlow[i] = -flowDict[path[i+1], path[i]]
            #print testFlow
            if np.any(testFlow>=0):
                gridPower[nodeNum]=0

            else:
                # power drawn from grid is min on path minus outflow
                gridP = np.max(testFlow)
                #print gridP
                for neigh in G.neighbors(start):
                    #print nodeLabelDict[start], nodeLabelDict[neigh]
                    try: 
                        checkFlow = flowDict[start, neigh]
                    except:
                        checkFlow = -flowDict[neigh, start]
                    #print checkFlow
                    if start and neigh in path:
                        #print 'in path to grid'
                        gridP = gridP + 0
                    else:
                        if checkFlow > 0:
                            gridP = gridP + checkFlow
                if gridP > 0:
                    gridPower[nodeNum] = 0
                else:
                    gridPower[nodeNum] = gridP
                    
    return gridPower
    
def get_grid_Export(G, nodeLabelDict, flowDict, PCC):

    gridPower = np.zeros((12))
    end = nodeLabelDict.keys()[0]
    gridPower = np.zeros((12))
    gridP = -PCC
    for neigh in G.neighbors(end):
        try: 
            checkFlow = flowDict[end, neigh]
        except:
            checkFlow = -flowDict[neigh, end]
        if checkFlow < 0:
            gridP = gridP + checkFlow
    if gridP < 0:
        gridPower[0] = 0
    else:
        gridPower[0] = gridP

    for nodeNum in np.arange(1,12):

        start = nodeLabelDict.keys()[nodeNum]
        #print nodeLabelDict[start], nodeLabelDict[end]

        for path in nx.all_simple_paths(G, source=start, target=end):
            #print path
            testFlow = np.zeros((len(path)))
            for i in range(len(testFlow)-1):
                #print path[i]
                #print path[i+1]
                try: 
                    testFlow[i] = flowDict[path[i], path[i+1]]
                except:
                    testFlow[i] = -flowDict[path[i+1], path[i]]
            testFlow[len(testFlow)-1] = -PCC
            if np.any(testFlow<=0):
                gridPower[nodeNum]=0

            else:
                # power exported to grid is min on path minus inflow INCLUDING PCC
                gridP = np.min(testFlow)
                #print gridP
                for neigh in G.neighbors(start):
                    #print nodeLabelDict[start], nodeLabelDict[neigh]
                    try: 
                        checkFlow = flowDict[start, neigh]
                    except:
                        checkFlow = -flowDict[neigh, start]
                    #print checkFlow
                    if start and neigh in path:
                        #print 'in path to grid'
                        gridP = gridP + 0
                    else:
                        if checkFlow < 0:
                            gridP = gridP + checkFlow
                if gridP < 0:
                    gridPower[nodeNum] = 0
                else:
                    gridPower[nodeNum] = gridP
                    
    return gridPower
    
def calculate_Payback(capCost, yearlyIncome, lifetime, replacementInterval, partsCost):
    cumulativeCosts = float( -capCost-partsCost )
    pperiod = float( lifetime )
    for i in range(1,lifetime):
        # at the end of year i
        cumulativeCosts = cumulativeCosts + yearlyIncome
        if cumulativeCosts >= 0:
            pperiod = (i) - (cumulativeCosts/yearlyIncome)
            #print i, cumulativeCosts/yearlyIncome
            break
        # if at the end of year i cumcosts are negative and part needs replaced then do    
        if i%replacementInterval == 0:
            cumulativeCosts = cumulativeCosts - partsCost
        
    return math.ceil(pperiod*100)/100

def calculate_Payback2(capCost, yearlyIncome, lifetime, replacementInterval, partsCost, yearlyCosts):
    cumulativeCosts = float( -capCost-partsCost )
    pperiod = float( lifetime )
    for i in range(1,lifetime):
        # at the end of year i
        cumulativeCosts = cumulativeCosts + yearlyIncome - yearlyCosts
        if cumulativeCosts >= 0:
            pperiod = (i) - (cumulativeCosts/yearlyIncome)
            #print i, cumulativeCosts/yearlyIncome
            break
        # if at the end of year i cumcosts are negative and part needs replaced then do    
        if i%replacementInterval == 0:
            cumulativeCosts = cumulativeCosts - partsCost
        
    return math.ceil(pperiod*100)/100
    
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