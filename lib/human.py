import sys
import re
import os
import copy as cp
from io import BytesIO
import datetime as dt
import time
import zipfile
import lib.filters as filters
import lib.graphic as graphic
import lib.utilities as ut
import lib.vitis as vit

autoSaveImg = False
ignoreEdgesIsoform = False
comprimeNode = False
typeDB = False #False = TCGA, True = Fantom
min_frel = 0
list_Genes = []
listBioNameUpdate = {}

#Build list of lists of gene
#Return a list of lists
def buildMatrixGenesHuman(listFilter, listFiles, fantom):
    global list_Genes
    global typeDB
    typeDB = fantom
    #BUILD MATRIX OF GENES
    matrixGenes = []
    extensionFiles = listFiles[0][-4:]
    for f in listFiles:
        if f[-4:] == extensionFiles and (f[-4:] != '.csv' or f[-4:] != '.zip'):
            if f[-4:] == '.zip':
                archive = zipfile.ZipFile(f, 'r')
                fileArchive = archive.namelist()
                for namefilezip in fileArchive:
                    if '.zip' in namefilezip:
                        subarchive = zipfile.ZipFile(BytesIO(archive.read(namefilezip)), 'r')
                        fileSubArchive = subarchive.namelist()
                        for namefilesubzip in fileSubArchive:
                            if 'expansion' in namefilesubzip:
                                csvText = str(subarchive.read(namefilesubzip))
                                csvText = csvText.split(r'"')
                                csvText = csvText[0]
                                csvListText = csvText.split(r'\n')
                                csvTemp = open(namefilesubzip, 'w')
                                for l in csvListText:
                                    csvTemp.write(l+'\n')
                                csvTemp.close()
                                listGenes = ut.readFilesVitis(namefilesubzip, False)
                                list_Genes = listGenes[1]
                                listGenes = listGenes[0]
                                os.remove(namefilesubzip)
                                #Filter lists
                                for filter in listFilter:
                                    listGenes = applyFilter(listGenes, filter)
                                matrixGenes.append(listGenes)
                    elif 'expansion' in namefilezip:
                        csvText = str(archive.read(namefilezip))
                        if csvText[1] == '"':
                            csvText = csvText.split(r'"')
                        else:
                            csvText = csvText.split('\'')
                        csvText = csvText[1]
                        csvListText = csvText.split(r'\n')
                        csvTemp = open(namefilezip, 'w')
                        for l in csvListText:
                            csvTemp.write(l+'\n')
                        csvTemp.close()
                        listGenes = ut.readFilesVitis(namefilezip, typeDB)
                        list_Genes = listGenes[1]
                        listGenes = listGenes[0]
                        os.remove(namefilezip)
                        #Filter lists
                        for filter in listFilter:
                            listGenes = applyFilter(listGenes, filter)
                        matrixGenes.append(listGenes)
                    elif 'csv' in namefilezip:
                        csvText = str(archive.read(namefilezip))
                        if csvText[1] == '"':
                            csvText = csvText.split(r'"')
                        else:
                            csvText = csvText.split('\'')
                        csvText = csvText[1]
                        csvListText = csvText.split(r'\n')
                        csvTemp = open(namefilezip, 'w')
                        for l in csvListText:
                            csvTemp.write(l+'\n')
                        csvTemp.close()
                        listGenes = ut.readFilesVitis(namefilezip, typeDB)
                        list_Genes = listGenes[1]
                        listGenes = listGenes[0]
                        os.remove(namefilezip)
                        #Filter lists
                        for filter in listFilter:
                            listGenes = applyFilter(listGenes, filter)
                        matrixGenes.append(listGenes)
            else:
                #Read gene files .csv
                listGenes = ut.readFilesVitis(f, typeDB)
                list_Genes = listGenes[1]
                listGenes = listGenes[0]
                #Filter lists
                for filter in listFilter:
                    listGenes = applyFilter(listGenes, filter)
                matrixGenes.append(listGenes)
        else:
            print('ERROR: FILES HAVE DIFFERENT EXTENSION. File need to have the same extension. All .csv or all .zip')
            sys.exit(-1)

    #READ UPDATE NAME GENE
    f = open('import_doc/couple_name_gene.csv', 'r')
    text = f.readlines()
    listLineName = []
    i = 0
    while i < len(text):
        listLineName.append(text[i].split(','))
        i += 1
    for l in listLineName:
        if l[0] != '':
            listBioNameUpdate[l[0]] = l[0]+'_'+l[1][:-1]
        for n in list_Genes:
            if n == l[0]:
                if l[0] != '':
                    if comprimeNode:
                        if (re.search(r'@\w*', l[1])).group() not in list_Genes:
                            list_Genes[list_Genes.index(n)] = (re.search(r'@\w*', l[1])).group()
                        else:
                            del(list_Genes[list_Genes.index(n)])
                    else:
                        list_Genes[list_Genes.index(n)] = l[0]+'_'+l[1][:-1]
    f.close()

    matrixGenesOld = []
    #comprime node if request
    if comprimeNode:
        #deep copy used to check edges in pearson correlation in comprime node version
        matrixGenesOld = cp.deepcopy(matrixGenes)
        #update name with name of only gene
        for l in matrixGenes:
            l[0] = (re.search(r'@\w*', listBioNameUpdate[l[0]])).group()
            i = 1
            while i < len(l):
                tmp = l[i]
                try:
                    l[i] = (tmp[0], (re.search(r'@\w*', listBioNameUpdate[tmp[1]])).group(), tmp[2])
                    i += 1
                except:
                    del(l[i])
        #unify list of same gene but different isoform
        i = 0
        j = 1
        while i < len(matrixGenes):
            while j < len(matrixGenes):
                if matrixGenes[i][0] == matrixGenes[j][0]:
                    matrixGenes[i] = matrixGenes[i]+matrixGenes[j][1:]
                    del(matrixGenes[j])
                else:
                    j += 1
            i += 1
            j = i+1
        #remove duplicate
        for l in matrixGenes:
            i = 1
            j = 2
            while i < len(l):
                while j < len(l):
                    if l[i][1] == l[j][1]:
                        #l[i] = (l[i][0], l[i][1], max(l[i][2], l[j][2]))
                        l[i] = (int((l[i][0]+l[j][0])/2), l[i][1], round(((l[i][2]+l[j][2])/2), 4))
                        del(l[j])
                    else:
                        j += 1
                i += 1
                j = i+1
            i = 1
            while i < len(l):
                if l[i][1] == l[0]:
                    del(l[i])
                else:
                    i += 1

    #return lists of gene lists filtered
    return (matrixGenes, matrixGenesOld, comprimeNode)

#Switch the filter to the correct function
#Return list of genes filtered
def applyFilter(listGenes, filter):
    if filter[0] == '-f' and len(filter) == 2:
        listGenes = filters.filterFrel(listGenes, float(filter[1]))
        global min_frel
        min_frel = float(filter[1])
    elif filter[0] == '-t':
        print('ERROR: incorrect parameters filter. -t only for Vitis')
        ut.printInfo()
    elif filter[0] == '-a' and len(filter) == 1:
        global autoSaveImg
        autoSaveImg = True
    elif filter[0] == '-i' and len(filter) == 2:
        global ignoreEdgesIsoform
        ignoreEdgesIsoform = True
        if filter[1] == 'comp':
            global comprimeNode
            comprimeNode = True
        elif filter[1] != 'notcomp':
            print('ERROR: incorrect parameters filter -i')
            ut.printInfo()
    #elif other filter
    else:
        print('ERROR: incorrect parameters filter')
        ut.printInfo()
    return listGenes

#Read file .csv
#Return list of tuples with position [0] name of the gene, in pos 2+ all gene associated
# def readFilesHuman(filename):
#     print('Open file: '+filename, flush=True)
#     try:
#         f = open(filename, 'rU')
#         text = f.read()
#         f.close()
#     except:
#         #except if file does not exist
#         print('ERROR: FILE NOT FOUND. File \''+filename+'\' does not exist')
#         sys.exit(-1)
#     #split rows
#     listLine = text.split('\n')
#     #split column based on symbol ','
#     listCells = []
#     for line in listLine:
#         listCells.append(line.split(','))
#     #Name of the gene is in cell in row 0, column 3
#     if typeDB:
#         #nameGene = (listCells[0][3].split('-'))[0]
#         nameGene = ((((listCells[0][0].split(' '))[3]).split('-'))[0]) #FANTOM
#     else:
#         nameGene = ((((listCells[0][0].split(' '))[3]).split('-'))[1])  #TCGA
#     #add at list of genes
#     list_Genes.append(nameGene.upper())
#     listTuples = [nameGene.upper()]
#     for i in listCells:
#         if not len(i) < 5:
#             try:
#                 #add new tuple: (rank, node, Frel)
#                 #possible add more parameters
#                 tuple = (int(i[0]), i[1].upper(), float(i[3]))
#                 listTuples.append(tuple)
#             except:
#                 pass
#     return listTuples

#
def indexDictGene(g, listStr):
    listIndex = []
    for elem in listStr:
        if g in elem:
            listIndex.append(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(elem)])
    return listIndex

#Print in output the graphs
def printOutput(coreGraph, graphGenes, graphGenesOld):
    if not os.path.exists('networkOutput'):
        os.mkdir('networkOutput')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('networkOutput/'+printTime)
    print('Creating directory: \'networkOutput/'+printTime+'\'', flush=True)

    #Pearson correlation
    print('Calculating Pearson Correlation...', flush=True)
    #only on complete graph because all edges in core are in the complete graph
    if comprimeNode:
        fIsoEdges = open('networkOutput/'+printTime+'/isoformInEdge.txt', 'w')
        fIsoEdgesCSV = open('networkOutput/'+printTime+'/isoformInEdge.csv', 'w')
        listEdgesOriginalGraph = [(a,b) for (a,b,c) in graphGenesOld]
        prePearson = {}
        for e in graphGenes:
            numberIsoformPerEdges = 0
            fIsoEdges.write('EDGE: '+str(e[0])+' - '+str(e[1])+'. Componed by:\n')
            fIsoEdgesCSV.write(str(e[0][1:])+','+str(e[1][1:]))
            listIndex1 = indexDictGene(e[0], listBioNameUpdate.values())
            listIndex2 = indexDictGene(e[1], listBioNameUpdate.values())
            for k in [(a, b) for a in listIndex1 for b in listIndex2]:
                if k in listEdgesOriginalGraph:
                    fIsoEdges.write('> '+str(listBioNameUpdate[k[0]])+' - '+str(listBioNameUpdate[k[1]])+'\n')
                    fIsoEdgesCSV.write(','+str(listBioNameUpdate[k[0]])+'<-->'+str(listBioNameUpdate[k[1]]))
                    numberIsoformPerEdges += 1
                    prePearson[k] = 1
            fIsoEdges.write('NUMBER OF EDGES: '+str(numberIsoformPerEdges)+'\n\n')
            fIsoEdgesCSV.write('\n')
        fIsoEdges.close()
        fIsoEdgesCSV.close()
        pearsonComplete = ut.pearsonCorrelation(list(prePearson.keys()), 'hgnc_cc_zero_filtered_mat.csv')
        i = 0
        while i < len(pearsonComplete):
            tmp = pearsonComplete[i]
            pearsonComplete[i] = ((re.search(r'@\w*', listBioNameUpdate[tmp[0]])).group(), (re.search(r'@\w*', listBioNameUpdate[tmp[1]])).group(), tmp[2])
            i += 1
        pearsonComplete = ut.manageDuplicates(pearsonComplete)
    else:
        pearsonComplete = [] #TCGA
        if typeDB:
            pearsonComplete = ut.pearsonCorrelation(graphGenes, 'hgnc_cc_zero_filtered_mat.csv') #FANTOM
        #update name of genes
        i = 0
        while i < len(coreGraph):
            tmp = coreGraph[i]
            if tmp[0] in listBioNameUpdate.keys() and tmp[1] in listBioNameUpdate.keys():
                coreGraph[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2])
            i += 1
        i = 0
        while i < len(graphGenes):
            tmp = graphGenes[i]
            if tmp[0] in listBioNameUpdate.keys() and tmp[1] in listBioNameUpdate.keys():
                graphGenes[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2])
            if typeDB:
                tmp = pearsonComplete[i] #FANTOM
                if tmp[0] in listBioNameUpdate.keys() and tmp[1] in listBioNameUpdate.keys():
                    pearsonComplete[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2]) #FANTOM
            else:
                if tmp[0] in listBioNameUpdate.keys() and tmp[1] in listBioNameUpdate.keys():
                    pearsonComplete[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], 1) #TCGA
                else:       #TCGA
                    pearsonComplete.append((tmp[0], tmp[1], 1))            #TCGA
            i += 1

    print('Calculating Pearson Correlation complete', flush=True)

    #RETURN CORE
    nameFileCore = 'networkOutput/'+printTime+'/Core_Graph'
    coreGraph = sorted(coreGraph, key=ut.ord)
    ut.printCSV(nameFileCore, coreGraph)
    #draw graph in a image
    graphic.drawGraph('H', coreGraph, nameFileCore+'_Circular', pearsonComplete, autoSaveImg, [], 1-min_frel, comprimeNode, False, typeDB)
    graphic.drawGraph('H', coreGraph, nameFileCore, pearsonComplete, autoSaveImg, [], 1-min_frel, comprimeNode, True, typeDB)

    #RETURN GRAPH
    nameFileCompleteGraph = 'networkOutput/'+printTime+'/Complete_Graph'
    graphGenes = sorted(graphGenes, key=ut.ord)
    ut.printCSV(nameFileCompleteGraph, graphGenes)
    #draw graph in a image
    graphic.drawGraph('H', graphGenes, nameFileCompleteGraph+'_Circular', pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, comprimeNode, False, typeDB)
    graphic.drawGraph('H', graphGenes, nameFileCompleteGraph, pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, comprimeNode, True, typeDB)

def removeIsoformEdges(edges):
    if ignoreEdgesIsoform and not comprimeNode:
        newEdges = []
        for e in edges:
            if e[0] in listBioNameUpdate.keys() and e[1] in listBioNameUpdate.keys():
                gene1 = re.search(r'@\w*', listBioNameUpdate[e[0]])
                gene2 = re.search(r'@\w*', listBioNameUpdate[e[1]])
                if gene1.group() != gene2.group():
                    newEdges.append(e)
        return newEdges
    return edges
