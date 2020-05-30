# #################################################################################################
#   Adapted version from: https://github.com/varungohil/Densest-Subgraph-Discovery                #
#   Authors  #  Mridul Sharma (mridul.sharma@iitgn.ac.in) , Heer Ambavi (heer.ambavi@iitgn.ac.in) #
###################################################################################################
##
#density of subgraph
def maxdensity(dic):
    a=0
    if(len(dic)==0):
        return a
    for i in dic.keys():
        #a+=len(dic[i])
        a+=(len(dic[i])*dicW[i])
    a=a/(2.0*len(dic))
    return a

#counter is minimum degree in subgraph
def f_pop(counter):
    tobedel=[]
    for i in deg.keys():
        if(deg[i]==0 or deg[i]*dicW[i] == counter):
            tobedel.append(i)
    #print dic
    for k in tobedel:
        for j in dic[k]:
            (dic[j]).remove(k)
            deg[j]=deg[j]-1
        del dic[k]
        del deg[k]

#averago of weight of each node
def avgWeight(idNode, completeGraph):
    dicW.clear()
    weight = []
    for k in idNode.keys():
        if idNode[k] in dic.keys():
            weight.append([(k, w) for (u, v, w) in completeGraph if (u == k and idNode[v] in dic.keys()) or (v == k and idNode[u] in dic.keys())])
    for k in weight:
        if len(k) > 0:
            avg = 0.0
            for t in k:
                avg += t[1]
            avg = avg/len(k)
            dicW[idNode[k[0][0]]] = avg

edges = []
dicW = {}
deg={}
dic={}

def findCoreNetwork(completeGraph):
    idNode = {}
    i = 1
    for k in completeGraph:
        if k[0] not in idNode:
            idNode[k[0]] = i
            i += 1
        if k[1] not in idNode:
            idNode[k[1]] = i
            i += 1
    for k in completeGraph:
        edges.append(str(idNode[k[0]])+' '+str(idNode[k[1]]))

    for i in edges:
        edge=list(map(int,i.split()))
        if edge[0] not in dic.keys():
            dic[edge[0]]=[]
            deg[edge[0]]=0
        dic[edge[0]].append(edge[1])
        deg[edge[0]]+=1

        if edge[1] not in dic.keys():
            dic[edge[1]]=[]
            deg[edge[1]]=0
        dic[edge[1]].append(edge[0])
        deg[edge[1]]+=1

    avgWeight(idNode, completeGraph)
    #main
    maxdens=0
    subgraph=[]
    while(len(deg)>0):
        mindeg=float('inf')
        for i in deg.keys():
            #add weight
            mindeg = min(deg[i]*dicW[i], mindeg)
        avgWeight(idNode, completeGraph)
        f_pop(mindeg)

        if (maxdensity(dic)>maxdens):
            maxdens=maxdensity(dic)
            subgraph=list(dic.keys())
    return subgraph
