import re
import sys
import scipy.stats.stats as stat
import lib.charikarAlgorithm as charikarAlgorithm

#Print info about cmd command if the call is wrong
def printInfo():
    print('Usage: python3 managerList.py PARAM [FILTERS]... -files [FILES]...')
    print('FILES can be a list of .csv or .zip')
    print('PARAM:')
    print('\t-vitis\tLists of vitis')
    print('\t-human\tLists of human. Need to be follow by:')
    print('\t\t-fantom\t\tLists from Fantom DB')
    print('\t\t-TCGA\t\tLists from TCGA DB')
    print('FILTERS:')
    print('\t-a\t\t\tAutosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS')
    print('\t-f [NUMBER]\t\tIgnored genes with frel<=NUMBER')
    print('\tONLY FOR VITIS GENES:')
    print('\t\t-t [PATTERN,...]\tTake genes that in \'functional annotation\' or \'Network1\' column there is at least one pattern')
    print('\tONLY FOR HUMAN GENES:')
    print('\t\t-i [\'comp\'/\'notcomp\']\tIgnored edges between isoforms of same gene')
    sys.exit(-1)

#Read the filter parameters and gene list name files
#Return list of lists of gene already filtered
def readParameters(input):
    listFilter = []
    listFiles = []
    i = 1
    while i < len(input):
        #pattern '^-\w$' is letter of a filter
        if re.search(r'^-\w$', input[i]):
            filterParam = (input[i],)
            i += 1
            endloop = False
            #loop to take all parameters of a filter
            while i < len(input) and endloop == False:
                #stop if find the next filter or add the params to the tuple of the filter
                if re.search(r'^-', input[i]):
                    endloop = True
                else:
                    filterParam = filterParam + (input[i],)
                    i += 1
            if filterParam[0] == input[i-1] and (filterParam[0] != '-a' and filterParam[0] != '-i'):
                #return an error if the filter has zero parameters and is not '-a'
                print('ERROR: incorrect number of parameters')
                printInfo()
            else:
                #add tuple to the list
                listFilter.append(filterParam)
        #if find pattern '^-files$', the following parameters are files
        elif re.search(r'^-files$', input[i]):
            i += 1
            while i < len(input):
                listFiles.append(input[i])
                i += 1
        elif re.search(r'^-vitis$|^-human$|^-fantom$|^-TCGA$', input[i]):
            i+=1
        #if is else, return an error
        else:
            print("ERROR: "+input[i])
            printInfo()
            i += 1
    #return an error if the list of file is empty
    if len(listFiles) == 0:
        print('ERROR: no files')
        printInfo()

    #print('FILTER: '+str(listFilter))
    #print('FILES: '+str(listFiles))
    return (listFilter, listFiles)

#Build matrix with only the name of the genes, to simplify the search of interaction
#Return a list of lists
def buildMatrixOnlyName(mGenes):
    mName = []
    for g in mGenes:
        i = 0
        lName = []
        for n in g:
            if i != 0:
                #name of the node is in second position of the tuple
                lName.append(n[1])
            else:
                #first is always the name of the start node
                lName.append(n)
            i += 1
        mName.append(lName)
    return mName

#Build the complete graph of interaction between genes
#Return list of tuples (nodeA, nodeB, frel)
def buildGraph(mGenes):
    i = 0
    j = 0
    #list of tuples (nodeA, nodeB, Frel)
    mGenesName = buildMatrixOnlyName(mGenes)
    graph = []
    for g in mGenes:
        for t in g:
            while j < len(mGenes) and i != j:
                # if len(t) == 5:
                #     #find indirect interaction
                #     if t[1] in mGenesName[j]:
                #         nodeB = mGenesName[j]
                #         lGene = mGenes[j]
                #         tGene = lGene[nodeB.index(t[1])]
                #         frel = tGene[2]
                #         nodeB = nodeB[0]
                #         if len(tGene) == 5:
                #             graph.append((g[0], t[1], frel))
                #             graph.append((t[1], nodeB, frel))
                # else:
                #find direct interaction
                if t in mGenesName[j]:
                    nodeB = mGenesName[j]
                    lGene = mGenes[j]
                    tGene = lGene[nodeB.index(t)]
                    frel = tGene[2]
                    nodeB = nodeB[0]
                    graph.append((t, nodeB, frel))
                j += 1
            j = 0
        i += 1
        j = 0
    return graph

#Trasform directed graph in undirected
#Return list of tuples (nodeA, nodeB, frel)
def manageDuplicates(completeGraph):
    minGraph = []
    i = 0
    j = 1
    edgeDone = []
    for edge1 in completeGraph:
        duplicate = False
        #check the two way of the edge
        edgeIN = edge1[0]+edge1[1]
        edgeOUT = edge1[1]+edge1[0]
        while j < len(completeGraph) and duplicate == False:
            edge2 = completeGraph[j]
            if edge1[0] == edge2[1] and edge1[1] == edge2[0] or edge1[0] == edge2[0] and edge1[1] == edge2[1]:
                duplicate = True
                #check if is already insert
                if edgeIN not in edgeDone or edgeOUT not in edgeDone:
                    edgeDone.append(edgeIN)
                    edgeDone.append(edgeOUT)
                    #insert only an edge with the average of frel
                    minGraph.append((edge1[0], edge1[1], round(((edge1[2]+edge2[2])/2), 4)))
            j += 1
        #check if exist edge that are not duplicates
        if duplicate == False and edgeIN not in edgeDone:
            minGraph.append(edge1)
        i += 1
        j = i+1
    return minGraph

#Find the most dense area of network
def findCoreGraph(completeGraph):
    #trasform name gene in idNode
    idNode = {}
    i = 1
    for k in completeGraph:
        if k[0] not in idNode:
            idNode[k[0]] = i
            i += 1
        if k[1] not in idNode:
            idNode[k[1]] = i
            i += 1
    #call Charikar's Algorithm to find the core
    lIdCore = charikarAlgorithm.findCoreNetwork(completeGraph)
    #transform idNode to name of gene
    lCore = []
    for elem in lIdCore:
        lCore.append(list(idNode.keys())[list(idNode.values()).index(elem)])
    #create list with edges between nodes of core
    core = []
    max = len(lCore)
    i = 0
    j = 1
    while i < len(lCore):
        wNode1 = lCore[i]
        while j < len(lCore):
            wNode2 = lCore[j]
            for edge in completeGraph:
                if wNode1 == edge[0] and wNode2 == edge[1] or wNode1 == edge[1] and wNode2 == edge[0]:
                    core.append(edge)
            j += 1
        i += 1
        j = i
    return core

#function to ord list of edges
def ord(t):
    return t[0]

#Calculate Pearson correlation
def pearsonCorrelation(edges, file):
    #calculate Pearson correlation
    listNameGene = []
    for k in edges:
        if k[0] not in listNameGene:
            listNameGene.append(k[0])
        if k[1] not in listNameGene:
            listNameGene.append(k[1])
    f = open('import_doc/'+file, 'r')
    lineCorr = []
    for line in f:
        listLine = line.split(',')
        listLine[-1] = listLine[-1][0:-1]
        for name in listNameGene:
            if name in listLine[0]:
                listFloat = [listLine[0]]
                for k in listLine[1:]:
                    listFloat.append(float(k))
                lineCorr.append(listFloat)
    f.close()
    edgeCorrValue = []
    for k in edges:
        line1 = []
        line2 = []
        i = 0
        while i < len(lineCorr):
            if k[0] in lineCorr[i][0]:
                line1 = lineCorr[i]
            if k[1] in lineCorr[i][0]:
                line2 = lineCorr[i]
            i += 1
        edgeCorrValue.append((line1[0], line2[0], (stat.pearsonr(line1[1:], line2[1:]))[0]))
    return edgeCorrValue


#print file CSV with edges of graph
def printCSV(nameFile, edges):
    f = open(nameFile+'.csv', 'w')
    f.write('NodeA,NodeB,frel\n')
    for k in edges:
        f.write(k[0]+','+k[1]+','+str(k[2])+'\n')
    f.close();
    print('Create: \''+nameFile+'.csv\'', flush=True)
