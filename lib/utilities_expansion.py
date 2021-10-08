import lib.utilities as ut
import lib.filters as filters
import re
import sys
import zipfile
from io import BytesIO
import os

#Print info about cmd command if the call is wrong
def printInfo():
    print('Usage: python3 integrateCoupleGenes.py PARAM TYPEA [FILTERS]... -files [GENES] [FILES]... [ISOFORM]')
    print('PARAM:')
    print('\t-vitis\tLists of vitis')
    print('\t-fantom\tLists of fantom DB')
    print('\t-TCGA\tLists of TCGA DB')
    print('TYPEA:')
    print('\t-frel\tBuild expansion network based on FREL. Required filter \'-f\'')
    print('\t-rank [INT]\tBuild expansion network based on RANK. Take top genes')
    print('\t-shared\tBuild expansion network based on SHARED GENES.')
    print('FILTERS:')
    print('\t-a\t\t\tAutosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS')
    print('\t-c\t\t\tAdd edges between associated genes')
    print('\t-e\t\t\tPrint Venn Diagram and Histogram for complete analysis')
    print('\t-f [NUMBER]\t\tIgnored genes with frel<=NUMBER')
    print('GENES: file .csv with the genes to analyze. Example: \'CoupleGeneToIntegrate/coupleGene.csv\'')
    print('FILES can be a list of .csv or .zip')
    print('ISOFORM: file .csv from the execution of ManagerList.py with the composition of edge gene-gene (tag \'-i comp\').\n\tTo use only with \'-fantom\'')
    sys.exit(-1)

#Manager to read and filter files
def readFilesGenes(listFiles, coupleGenes, listfilter, vitis, TCGAdb, listBioNameUpdate):
    dictGeneToAnalyze = {}
    #Read only file of a gene to analyze
    for f in coupleGenes:
        for elem in f:
            dictGeneToAnalyze[elem] = 0
    matrixGenes = []
    extensionFiles = listFiles[0][-4:]
    for f in listFiles:
        if f[-4:] == extensionFiles and (f[-4:] != '.csv' or f[-4:] != '.zip'):
            if f[-4:] == '.zip':
                archive = zipfile.ZipFile(f, 'r')
                fileArchive = archive.namelist()
                for namefilezip in fileArchive:
                    if '.zip' in namefilezip:
                        #Read if is a subzip
                        subarchive = zipfile.ZipFile(BytesIO(archive.read(namefilezip)), 'r')
                        fileSubArchive = subarchive.namelist()
                        for namefilesubzip in fileSubArchive:
                            if 'expansion' in namefilesubzip:
                                csvText = str(subarchive.read(namefilesubzip))
                                csvText = csvText.split(r'"')
                                csvText = csvText[0]
                                csvListText = csvText.split(r'\n')
                                nameGene = ((re.search(r'-\w*\s', csvListText[0])).group())[1:-1]
                                if nameGene in dictGeneToAnalyze.keys():
                                    csvTemp = open(namefilesubzip, 'w')
                                    for l in csvListText:
                                        csvTemp.write(l+'\n')
                                    csvTemp.close()
                                    listGenes = readFilesHuman(namefilesubzip, TCGAdb)
                                    os.remove(namefilesubzip)
                                    # #Filter lists
                                    for filter in listfilter:
                                        listGenes = applyFilter(listGenes, filter)
                                    matrixGenes.append(listGenes)
                    elif 'expansion' in namefilezip:
                        #Read if is a list of human
                        csvText = str(archive.read(namefilezip))
                        csvText = csvText.split(r'"')
                        csvText = csvText[0]
                        csvListText = csvText.split(r'\n')
                        nameGene = ((re.search(r'-\w*\s', csvListText[0])).group())[1:-1]
                        if nameGene in dictGeneToAnalyze.keys():
                            csvTemp = open(namefilezip, 'w')
                            for l in csvListText:
                                csvTemp.write(l+'\n')
                            csvTemp.close()
                            listGenes = readFilesHuman(namefilezip, TCGAdb)
                            os.remove(namefilezip)
                            #Filter lists
                            for filter in listfilter:
                                listGenes = applyFilter(listGenes, filter)
                            matrixGenes.append(listGenes)
                    elif 'csv' in namefilezip:
                        #Read file csv inside an archive zip
                        csvText = str(archive.read(namefilezip))
                        if csvText[1] == '"':
                            csvText = csvText.split(r'"')
                        else:
                            csvText = csvText.split('\'')
                        csvListText = []
                        nameGene = ''
                        #different split if is vitis or human
                        if vitis:
                            csvText = csvText[1]
                            csvListText = csvText.split(r'\n')
                            try:
                                nameGene = listBioNameUpdate[((csvListText[0].split(r','))[3]).upper()]
                            except:
                                listBioNameUpdate[((csvListText[0].split(r','))[3]).upper()] = ((csvListText[0].split(r','))[3]).upper()
                                nameGene = ((csvListText[0].split(r','))[3]).upper()
                        else:
                            csvText = csvText[0]
                            csvListText = csvText.split(r'\n')
                            nameGene = (((re.search(r'-\w*\s', csvListText[0])).group())[1:-1]).upper()
                        #if is a gene to analyze read it
                        if nameGene in dictGeneToAnalyze.keys():
                            csvTemp = open(namefilezip, 'w')
                            for l in csvListText:
                                csvTemp.write(l+'\n')
                            csvTemp.close()
                            listGenes = []
                            if vitis:
                                listGenes = (ut.readFilesVitis(namefilezip,True))[0]
                                if nameGene == 'GT-001':
                                    listGenes[0] = nameGene #TODEL_GT-001
                            else:
                                listGenes = readFilesHuman(namefilezip, TCGAdb)
                            os.remove(namefilezip)
                            # #Filter lists
                            for filter in listfilter:
                                listGenes = applyFilter(listGenes, filter)
                            matrixGenes.append(listGenes)
            else:
                #Read gene files .csv
                fileRead = open(f, 'r')
                csvText = fileRead.read()
                csvText = csvText.split(r'"')
                csvText = csvText[0]
                csvListText = csvText.split(r'\n')
                nameGene = ((re.search(r'-\w*\s', csvListText[0])).group())[1:-1]
                if nameGene in dictGeneToAnalyze.keys():
                    listGenes = readFilesHuman(f, TCGAdb)
                    #Filter lists
                    for filter in listfilter:
                        listGenes = applyFilter(listGenes, filter)
                    matrixGenes.append(listGenes)
        else:
            print('ERROR: FILES HAVE DIFFERENT EXTENSION. File need to have the same extension. All .csv or all .zip')
            sys.exit(-1)

    return matrixGenes

#Create edges list of common genes
# (a, b, c, d) --> (a, c, d)
def buildEdges(listCommonGenes):
    edges = []
    for e in listCommonGenes[1:]:
        edges.append((e[0], e[2], float(e[3])))
    return edges

#Create string from the combination of name genes
def buildNamefile(l):
    nameF = ''
    for g in sorted(l[0]):
        if nameF == '':
            nameF = g
        else:
            nameF += '_'+g
    return nameF

#Switch the filter to the correct function
#Return list of genes filtered
def applyFilter(listGenes, filter):
    if filter[0] == '-f' and len(filter) == 2:
        listGenes = filters.filterFrel(listGenes, float(filter[1]))
    elif filter[0] == '-rank' and len(filter) == 2:
        listGenes = filters.filterRank(listGenes, int(filter[1]))
    elif (filter[0] == '-a' or filter[0] == '-c' or filter[0] == '-e') and len(filter) >= 1:
        #already managed
        pass
    else:
        print('ERROR: incorrect parameters filter '+str(filter))
        printInfo()
    return listGenes

#Read file .csv
#Return list of tuples with position [0] name of the gene, in pos 2+ all gene associated
def readFilesHuman(filename, TCGAdb):
    print('Open file: '+filename, flush=True)
    try:
        f = open(filename, 'rU')
        text = f.read()
        f.close()
    except:
        #except if file does not exist
        print('ERROR: FILE NOT FOUND. File \''+filename+'\' does not exist')
        sys.exit(-1)
    #split rows
    listLine = text.split('\n')
    #split column based on symbol ','
    listCells = []
    for line in listLine:
        listCells.append(line.split(','))
    #Name of the gene is in cell in row 0, column 3
    nameGene = ((((listCells[0][0].split(' '))[3]).split('-'))[0])
    if TCGAdb:
        nameGene = ((((listCells[0][0].split(' '))[3]).split('-'))[1])  ##todel
    #add at list of genes
    #list_Genes.append(nameGene.upper())
    listTuples = [nameGene.upper()]
    for i in listCells:
        if len(i) > 4:
            try:
                #add new tuple: (rank, node, Frel)
                #possible add more parameters
                tuple = (int(i[0]), i[1].upper(), float(i[3]))
                listTuples.append(tuple)
            except:
                pass
    return listTuples

#Read which gene need to compare
def readFiles(file):
    print('Open file: '+file, flush=True)
    try:
        f = open(file, 'rU')
        text = f.read()
        f.close()
    except:
        #except if file does not exist
        print('ERROR: FILE NOT FOUND. File \''+file+'\' does not exist')
        sys.exit(-1)
    #split rows
    listLine = text.split('\n')
    listTuples = []
    i = 0
    while i < len(listLine):
        try:
            if listLine[i] != '':
                listTuples.append(listLine[i].split(','))
        except:
            pass
        i += 1
    return listTuples

#Read name of associated genes
#Create a list with associated genes with no repetition
def nameAssociateGene(g, nameUpdate, TCGAdb):
    alreadyRead = {}
    listGenes = []
    for elem in g[1:]:
        if elem[2] not in g[0] and elem[2] not in alreadyRead.keys():
            listGenes.append(elem[2])
            alreadyRead[elem[2]] = 1
    if not TCGAdb:
        i = 0
        while i < len(listGenes):
            listGenes[i] = (list(nameUpdate.keys())[list(nameUpdate.values()).index(listGenes[i])])
            i += 1
    return listGenes

#Print in a file the numbers of gene shared
def printNumberVenn(listCommonGenes, nameDir):
    i = 0
    while i < len(listCommonGenes[0]):
        nameF = buildNamefile(listCommonGenes[0][i])
        #create dir for each couple of genes
        nameF = nameF.replace("<", "_")
        nameF = nameF.replace(">", "_")
        nameDirGenes = nameDir+str(i)+'/'
        f = open(nameDirGenes+'venn_number.txt', 'w')
        f.write('---> NUMBER GENES IN \''+nameF+'\': '+str(int(len([u for u in listCommonGenes[0][i][1:] if not (u[0] in listCommonGenes[0][i][0] and u[2] in listCommonGenes[0][i][0])])/len(listCommonGenes[0][i][0])))+'\n')
        listKey = sorted([u for u in listCommonGenes[1][i].keys()], key=len)
        #listKey = sorted([u for u in listCommonGenes[1][i].keys() if len(u.split(',')) == 1], key=len)
        #listKey = sorted([(re.findall(r'\w+', u)) for u in listCommonGenes[1][i].keys()], key=len)
        for key in listKey[:-1]:
            nameFile = ''
            for g in sorted(key.split('\'')):
                if len(g) > 4:
                    if nameFile == '':
                        nameFile = g
                    else:
                        nameFile += '_'+g

            f.write('---> NUMBER GENES IN \''+nameFile+'\': '+str(int(len(listCommonGenes[1][i][str(key)])/len(nameFile.split('_'))))+'\n')
        f.close()
        print('Create: \''+nameDirGenes+'venn_number.txt\'', flush=True)
        i += 1

#Trasform directed graph in undirected
#Return list of tuples (nodeA, rank, nodeB, frel)
def manageDuplicates(completeGraph):
    minGraph = []
    i = 0
    j = 1
    edgeDone = []
    for edge1 in completeGraph:
        duplicate = False
        #check the two way of the edge
        edgeIN = edge1[0]+edge1[2]
        edgeOUT = edge1[2]+edge1[0]
        while j < len(completeGraph) and duplicate == False:
            edge2 = completeGraph[j]
            if edge1[0] == edge2[2] and edge1[2] == edge2[0] or edge1[0] == edge2[0] and edge1[2] == edge2[2]:
                duplicate = True
                #check if is already insert
                if edgeIN not in edgeDone or edgeOUT not in edgeDone:
                    edgeDone.append(edgeIN)
                    edgeDone.append(edgeOUT)
                    #insert only an edge with the average of frel
                    minGraph.append((edge1[0], 'avg_'+str(int((edge1[1]+edge2[1])/2)), edge1[2], round(((edge1[3]+edge2[3])/2), 4)))
            j += 1
        #check if exist edge that are not duplicates
        if duplicate == False and edgeIN not in edgeDone:
            minGraph.append(edge1)
        i += 1
        j = i+1
    return minGraph

#Function for sorted list
def ord2(tuple):
    return tuple[2]

#Find common genes in lists of genes
def findCommonGenes(couples, listFiles):
    #Create dict based on name to know the index of the gene
    listDictGenesFiles = {}
    for l in listFiles:
        dictGenesFiles = {}
        i = 1
        while i < len(l):
            dictGenesFiles[l[i][1]] = i
            i += 1
        listDictGenesFiles[l[0]] = (listFiles.index(l), dictGenesFiles)
    #find common genes
    listCommonGenes = []
    listForVenn = []
    for c in couples:
        innerListGenes = [c]
        edgeBetweenGenesLGN = []
        innerListForVenn = {}
        dictGeneToSave = {}
        i = 0
        #create a dictionary with inside lists with in position 0 the genes that belong to
        # and follow by all genes in that list
        while i < len(listFiles):
            if listFiles[i][0] in c:
                for elem in listFiles[i][1:]:
                    nGene = elem[1]
                    if nGene in c:
                        edgeBetweenGenesLGN.append((listFiles[i][0], elem[0], elem[1], elem[2]))
                    dictGeneToSave[nGene] = []
                    for g in c:
                        if nGene in ((listDictGenesFiles[g])[1]).keys():
                            dictGeneToSave[nGene] = dictGeneToSave[nGene]+[g]
            i += 1
        #create innerListForVenn[] for every possible combination of the genes in couples
        keys = []
        for g in c:
            keys.append([g])
        i = 0
        j = 1
        k = 0
        while i < len(c):
            while j < len(c):
                tmpL = [c[i]]
                while j < len(c):
                    tmpL.append(c[j])
                    keys.append(tmpL.copy())
                    j+=1
                k+=1
                j=i+1+k
            i+=1
            j=i+1
        #Add each possible combination of genes to innerListForVenn
        for k in keys:
            innerListForVenn[str(sorted(k))] = []
        edgeBetweenGenesLGN = manageDuplicates(edgeBetweenGenesLGN)
        innerListGenes = innerListGenes+edgeBetweenGenesLGN
        for key in dictGeneToSave.keys():
            if len(dictGeneToSave[key]) == len(c):
                for g in c:
                    try:
                        tmpCouple = listFiles[listDictGenesFiles[g][0]][((listDictGenesFiles[g])[1])[key]]
                        innerListGenes.append((g, tmpCouple[0], tmpCouple[1], tmpCouple[2]))
                        innerListForVenn[str(sorted(dictGeneToSave[key]))] += [(g, tmpCouple[0], tmpCouple[1], tmpCouple[2])]
                    except:
                        pass
            i = 1
            #divide nodes for venn diagram
            while len(c)-i >= 0:
                if len(dictGeneToSave[key]) == len(c)-i:
                    for g in c:
                        try:
                            tmpCouple = listFiles[listDictGenesFiles[g][0]][((listDictGenesFiles[g])[1])[key]]
                            innerListForVenn[str(sorted(dictGeneToSave[key]))] += [(g, tmpCouple[0], tmpCouple[1], tmpCouple[2])]
                        except:
                            pass
                i+=1
        listForVenn.append(innerListForVenn)
        innerListGenes = [innerListGenes[0]]+sorted(innerListGenes[1:], key=ord2)
        listCommonGenes.append(innerListGenes)

    return (listCommonGenes, listForVenn)

#Find common isoform in genes from Fantom DB
def findCommonGenesFantom(couples, listFiles, isoformInEdge):
    #Create dict based on name to know the index of the gene
    listDictGenesFiles = {}
    for l in listFiles:
        dictGenesFiles = {}
        i = 1
        while i < len(l):
            dictGenesFiles[l[i][1]] = i
            i += 1
        listDictGenesFiles[l[0]] = (listFiles.index(l), dictGenesFiles)
    #find common genes
    listCommonGenes = []
    listForVenn = []
    #Based on the file that say which edges isoform-isoform is in edge gene-gene
    for c in couples:
        #find isoform we need to use
        edgesNodes = [(u[0], u[1]) for u in isoformInEdge]
        i = 0
        j = 1
        isoformToSearch = []
        edgeBetweenGenesLGN = []
        numberEdgesBetweenGeneCouple = 0
        while i < len(c):
            # numberEdgesBetweenGeneCouple += len(c)-i-1
            while j < len(c):
                if (c[i], c[j]) in edgesNodes:
                    isoformToSearch.append((isoformInEdge[edgesNodes.index((c[i], c[j]))])[2:])
                    numberEdgesBetweenGeneCouple += 1
                if (c[j], c[i]) in edgesNodes:
                    isoformToSearch.append((isoformInEdge[edgesNodes.index((c[j], c[i]))])[2:])
                    numberEdgesBetweenGeneCouple += 1
                j += 1
            i += 1
            j = i+1
        # if numberEdgesBetweenGeneCouple == len(isoformToSearch):
        innerListGenes = [c]
        innerListForVenn = {}
        dictGeneToSave = {}
        listNameIsoform = {}
        #save name isoform of lists in a dictionary to improve the performance
        for l in isoformToSearch:
            for n in l:
                tmpIsoform = n.split('<-->')
                for iso in tmpIsoform:
                    listNameIsoform[iso] = ((re.search(r'@\w*', iso)).group())[1:]
        for f in listFiles:
            if f[0] in listNameIsoform:
                for elem in f[1:]:
                    nGene = elem[1]
                    if (nGene.split('@'))[1] in c and (f[0].split('@'))[1] != (elem[1].split('@'))[1]:
                        edgeBetweenGenesLGN.append(((f[0].split('@'))[1], elem[0], (elem[1].split('@'))[1], elem[2]))
                    dictGeneToSave[nGene] = []
                    for g in c:
                        if nGene in ((listDictGenesFiles[list(listNameIsoform.keys())[list(listNameIsoform.values()).index(g)]])[1]).keys():
                            dictGeneToSave[nGene] = dictGeneToSave[nGene]+[g]
        #for key in dictGeneToSave.keys():
        #create innerListForVenn[] for every possible combination of the genes in couples
        keys = []
        for g in c:
            keys.append([g])
        i = 0
        j = 1
        k = 0
        while i < len(c):
            while j < len(c):
                tmpL = [c[i]]
                while j < len(c):
                    tmpL.append(c[j])
                    keys.append(tmpL.copy())
                    j+=1
                k+=1
                j=i+1+k
            i+=1
            j=i+1
        for k in keys:
            innerListForVenn[str(sorted(k))] = []
        edgeBetweenGenesLGN = manageDuplicates(edgeBetweenGenesLGN)
        innerListGenes = innerListGenes+edgeBetweenGenesLGN
        #Divide genes to the correct innerList
        for key in dictGeneToSave.keys():
            if len(dictGeneToSave[key]) == len(c):
                for g in c:
                    try:
                        tmpCouple = listFiles[listDictGenesFiles[[list(listNameIsoform.keys())[list(listNameIsoform.values()).index(g)]][0]][0]][((listDictGenesFiles[[list(listNameIsoform.keys())[list(listNameIsoform.values()).index(g)]][0]])[1])[key]]
                        innerListGenes.append((g, tmpCouple[0], tmpCouple[1], tmpCouple[2]))
                        innerListForVenn[str(sorted(dictGeneToSave[key]))] += [(g, tmpCouple[0], tmpCouple[1], tmpCouple[2])]
                    except:
                        pass
            #divide nodes for venn diagram
            i = 1
            while len(c)-i >= 0:
                if len(dictGeneToSave[key]) == len(c)-i:
                    for g in c:
                        try:
                            tmpCouple = listFiles[listDictGenesFiles[[list(listNameIsoform.keys())[list(listNameIsoform.values()).index(g)]][0]][0]][((listDictGenesFiles[[list(listNameIsoform.keys())[list(listNameIsoform.values()).index(g)]][0]])[1])[key]]
                            innerListForVenn[str(sorted(dictGeneToSave[key]))] += [(g, tmpCouple[0], tmpCouple[1], tmpCouple[2])]
                        except:
                            pass
                i+=1
        listForVenn.append(innerListForVenn)
        innerListGenes = [innerListGenes[0]]+sorted(innerListGenes[1:], key=ord2)
        listCommonGenes.append(innerListGenes)
    return (listCommonGenes, listForVenn)

#Build Edges for frel and rank
def buildEdgesFrelRank(listCouple, listFiles):
    listsEdges = []
    for lgn in listCouple:
        innerListEdges = [lgn]
        for l in listFiles:
            if l[0] in lgn:
                for genes in l[1:]:
                    innerListEdges.append((l[0], genes[0], genes[1], genes[2]))
        listsEdges.append(innerListEdges)
    return listsEdges

#Build Edges for frel and rank
def buildEdgesFrelRankIsoform(listCouple, listFiles, isoformInEdge):
    listsEdges = []
    for c in listCouple:
        #find isoform we need to use
        edgesNodes = [(u[0], u[1]) for u in isoformInEdge]
        i = 0
        j = 1
        isoformToSearch = []
        numberEdgesBetweenGeneCouple = 0
        while i < len(c):
            numberEdgesBetweenGeneCouple += len(c)-i-1
            while j < len(c):
                if (c[i], c[j]) in edgesNodes:
                    isoformToSearch.append((isoformInEdge[edgesNodes.index((c[i], c[j]))])[2:])
                if (c[j], c[i]) in edgesNodes:
                    isoformToSearch.append((isoformInEdge[edgesNodes.index((c[j], c[i]))])[2:])
                j += 1
            i += 1
            j = i+1
        listNameIsoform = {}
        #save name isoform of lists in a dictionary to improve the performance
        for l in isoformToSearch:
            for n in l:
                tmpIsoform = n.split('<-->')
                for iso in tmpIsoform:
                    listNameIsoform[iso] = ((re.search(r'@\w*', iso)).group())[1:]
        #TODO FIX
        innerListEdges = [c]
        for l in listFiles:
            if l[0] in listNameIsoform.keys():
                for genes in l[1:]:
                    innerListEdges.append((listNameIsoform[l[0]], genes[0], genes[1], genes[2]))
        listsEdges.append(innerListEdges)
    return listsEdges

#Manage line with multiple genes
def manageBR(l):
    for u in l[1:]:
        if '<BR>' in u[1]:
            names = u[1].split('<BR>')
            tmp = list(u)
            tmp[1] = names[0]
            l[l.index(u)] = tuple(tmp)
            for name in names[1:]:
                tmp = list(u)
                tmp[1] = name
                l.append(tuple(tmp))
    return l

#
def printCSV(edgesGraph, listForVenn, nameDir,listBioNameUpdate):
    #Read information of Vitis genes
    f = open('import_doc/NewAnnotVitisnet3.csv', 'r')
    text = f.readlines()
    lineIntro = str(text[0].split(',')[0])+','+str(text[0].split(',')[1])+','+str(text[0].split(',')[2])+','+str(text[0].split(',')[3])+','+str(text[0].split(',')[4])+','+str(text[0].split(',')[5][:-1])
    i = 1
    dictStrToWrite = {}
    while i < len(text):
        u = text[i].split(',')
        dictStrToWrite[str(u[0]).upper()] = str(u[0])+','+str(u[1])+','+str(u[2])+','+str(u[3])+','+str(u[4])+','+str(u[5][:-1])
        i += 1

    for k in edgesGraph:
        nameF = buildNamefile(k)
        #create dir for each couple of genes
        nameDirGenes = nameDir+str(edgesGraph.index(k))+'/'
        nameF = nameF.replace("<", "_")
        nameF = nameF.replace(">", "_")
        os.mkdir(nameDirGenes)
        #Write file .csv
        f = open(nameDirGenes+'edges_graph'+'.csv', 'w')
        f.write(str(nameF)+',rank,frel,'+lineIntro+'\n')
        for elem in k[1:]:
            try:
                f.write(str(elem[0])+','+str(elem[1])+','+str(elem[3])+','+dictStrToWrite[(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(elem[2])])]+'\n')
            except:
                f.write(str(elem[0])+','+str(elem[1])+','+str(elem[3])+','+str(elem[2])+'\n')#
                #f.write(str(elem[0])+','+str(elem[1])+','+str(elem[3])+','+elem[2]+'\n')
            #
            # string = str(elem[0])+','+str(elem[1])+','+str(elem[2])+','+str(elem[3])+'\n'
            # f.write(string)
        f.close();

    for k in listForVenn:
        listKey = sorted([u for u in k.keys()], key=len)
        nameF = ''
        lenMax = len(sorted([u for u in listKey if len(u.split(',')) == 1]))
        for g in sorted([u for u in listKey if len(u.split(',')) == 1]):
            if nameF == '':
                nameF = g.split('\'')[1]
            else:
                nameF += '_'+g.split('\'')[1]
        #Print .csv for each combination of genes of LGN. They contain the list of genes associated to that combination
        #FIX
        i = 2
        dictNumFile = {}
        while i < lenMax:
            dictNumFile[i] = 1
            i += 1
        for key in listKey:
            if len(key.split(',')) != lenMax:
                if len(key.split(',')) == 1:
                    nameFile = ''
                    for g in sorted(key.split('\'')):
                        if len(g) > 4:
                            if nameFile == '':
                                nameFile = g
                            else:
                                nameFile += '_'+g
                    nameF = nameF.replace("<", "_")
                    nameF = nameF.replace(">", "_")
                    nameFile = nameFile.replace("<", "_")
                    nameFile = nameFile.replace(">", "_")
                    f = open(nameDir+str(listForVenn.index(k))+'/'+nameFile+'.csv', 'w')
                else:
                    f = open(nameDir+str(listForVenn.index(k))+'/intersectionOF'+str(len(key.split(',')))+'genes'+str(dictNumFile[len(key.split(','))])+'.csv', 'w')
                    dictNumFile[len(key.split(','))] = dictNumFile[len(key.split(','))]+1
                f.write(nameFile+',rank,frel,'+lineIntro+'\n')
                for elem in k[str(key)]:
                    try:
                        f.write(str(elem[0])+','+str(elem[1])+','+str(elem[3])+','+dictStrToWrite[(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(elem[2])])]+'\n')
                    except:
                        f.write(str(elem[0])+','+str(elem[1])+','+str(elem[3])+','+str(elem[2])+'\n')

                f.close()
