#Apply filter on Frel
#Remove all gene with Frel <= perc
def filterFrel(listGenes, perc):
    newListGenes = []
    for g in listGenes:
        try:
            if g[2] > perc:
                newListGenes.append(g)
        except:
            #try/except if the line is the first (name gene)
            newListGenes.append(g)
    return newListGenes

#Apply filter on rank
#Remove all gene with Frel > rank
def filterRank(listGenes, rank):
    newListGenes = []
    for g in listGenes:
        try:
            if g[0] <= rank:
                newListGenes.append(g)
        except:
            #try/except if the line is the first (name gene)
            newListGenes.append(g)
    return newListGenes

#Return all tuples with at least one pattern inside 'functional annotation' or  'Network1' columns
def filterType(listGenes, patterns):
    newListGenes = [listGenes[0]]
    i = 1
    while i < len(listGenes):
        g = listGenes[i]
        toSave = False
        for p in patterns:
            try:
                if p in g[1]:
                    toSave = True
                if p in g[3]:
                    toSave = True
                if p in g[4]:
                    toSave = True
            except:
                pass
                #try/except if the line is the first (name gene)
                #print('EXCEPT')
                #newListGenes.append(g)
        if toSave:
            newListGenes.append(g)
        i += 1
    return newListGenes
