import sys
import re
import os
import datetime as dt
import time
import zipfile
import lib.filters as filters
import lib.graphic as graphic
import lib.utilities as ut

autoSaveImg = False
min_frel = 0
list_Genes = []

#Build list of lists of gene
#Return a list of lists
def buildMatrixGenesVitis(listFilter, listFiles):
    global list_Genes
    matrixGenes = []
    extensionFiles = listFiles[0][-4:]
    for f in listFiles:
        if f[-4:] == extensionFiles and (f[-4:] != '.csv' or f[-4:] != '.zip'):
            if f[-4:] == '.zip':
                archive = zipfile.ZipFile(f, 'r')
                fileArchive = archive.namelist()
                for namefilezip in fileArchive:
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
                    listGenes = ut.readFilesVitis(namefilezip, True)
                    list_Genes = listGenes[1]
                    listGenes = listGenes[0]
                    os.remove(namefilezip)
                    #Filter lists
                    for filter in listFilter:
                        listGenes = applyFilter(listGenes, filter)
                    matrixGenes.append(listGenes)
            else:
                #Read gene files .csv
                listGenes = ut.readFilesVitis(namefilezip, True)
                list_Genes = listGenes[1]
                listGenes = listGenes[0]
                #Filter lists
                for filter in listFilter:
                    listGenes = applyFilter(listGenes, filter)
                matrixGenes.append(listGenes)
        else:
            print('ERROR: FILES HAVE DIFFERENT EXTENSION. File need to have the same extension. All .csv or all .zip')
            sys.exit(-1)
    #return lists of gene lists filtered
    return matrixGenes

#Switch the filter to the correct function
#Return list of genes filtered
def applyFilter(listGenes, filter):
    if filter[0] == '-f' and len(filter) == 2:
        listGenes = filters.filterFrel(listGenes, float(filter[1]))
        global min_frel
        min_frel = float(filter[1])
    elif filter[0] == '-t' and len(filter) >= 2:
        listGenes = filters.filterType(listGenes, filter[1:])
    elif filter[0] == '-a' and len(filter) == 1:
        global autoSaveImg
        autoSaveImg = True
    elif filter[0] == '-i':
        print('ERROR: incorrect parameters filter. -i only for Human')
        ut.printInfo()
    #elif other filter
    else:
        print('ERROR: incorrect parameters filter')
        ut.printInfo()
    return listGenes


#Print in output the graphs
def printOutput(coreGraph, graphGenes):
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
    pearsonComplete = ut.pearsonCorrelation(graphGenes, 'vv_exprdata_2.csv')

    #UPDATE NAME GENE
    f = open('import_doc/NewAnnotVitisnet3.csv', 'r')
    text = f.readlines()
    listLineName = []
    i = 1
    while i < len(text):
        listLineName.append(text[i].split(','))
        i += 1
    listBioNameUpdate = {}
    for l in listLineName:
        if l[3] != '':
            listBioNameUpdate[l[0].upper()] = l[3]
        elif l[2] != '':
            listBioNameUpdate[l[0].upper()] = l[2]
        else:
            listBioNameUpdate[l[0].upper()] = l[0].upper()
        for n in list_Genes:
            if n == l[0].upper():
                if l[3] != '':
                    list_Genes[list_Genes.index(n)] = l[3]
                elif l[2] != '':
                    list_Genes[list_Genes.index(n)] = l[2]
                else:
                    list_Genes[list_Genes.index(n)] = l[0].upper()
    f.close()
    #update name of genes
    i = 0
    while i < len(coreGraph):
        tmp = coreGraph[i]
        coreGraph[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2])
        i += 1
    i = 0
    while i < len(graphGenes):
        tmp = graphGenes[i]
        graphGenes[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2])
        tmp = pearsonComplete[i]
        pearsonComplete[i] = (listBioNameUpdate[tmp[0]], listBioNameUpdate[tmp[1]], tmp[2])
        i += 1

    #RETURN CORE
    nameFileCore = 'networkOutput/'+printTime+'/Core_Graph'
    coreGraph = sorted(coreGraph, key=ut.ord)
    ut.printCSV(nameFileCore, coreGraph)
    #draw graph in a image
    graphic.drawGraph('V', coreGraph, nameFileCore+'_Circular', pearsonComplete, autoSaveImg, [], 1-min_frel, False, False, True)
    graphic.drawGraph('V', coreGraph, nameFileCore, pearsonComplete, autoSaveImg, [], 1-min_frel, False, True, True)

    #RETURN GRAPH
    nameFileCompleteGraph = 'networkOutput/'+printTime+'/Complete_Graph'
    graphGenes = sorted(graphGenes, key=ut.ord)
    ut.printCSV(nameFileCompleteGraph, graphGenes)
    #draw graph in a image
    graphic.drawGraph('V', graphGenes, nameFileCompleteGraph+'_Circular', pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, False, False, True)
    graphic.drawGraph('V', graphGenes, nameFileCompleteGraph, pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, False, True, True)
