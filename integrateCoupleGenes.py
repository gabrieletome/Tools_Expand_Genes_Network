import sys
import re
import zipfile
from io import BytesIO
import os
import datetime as dt
import lib.vitis as vit
import lib.human as hum
import lib.utilities_expansion as utex
import lib.components_graph as comp
import lib.graphic as graphic
import lib.utilities as ut

listBioNameUpdate = {}
TCGAdb = True
vitis = False
autoSaveImg = False
completeGraph = False
min_frel = 0
nameDir = 'tmpName'

#Manager to read and filter files
def readFilesGenes(listFiles, coupleGenes, listfilter):
    global min_frel
    global autoSaveImg
    global completeGraph
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
                                    listGenes = utex.readFilesHuman(namefilesubzip, TCGAdb)
                                    os.remove(namefilesubzip)
                                    # #Filter lists
                                    for filter in listfilter:
                                        listGenes = utex.applyFilter(listGenes, filter)
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
                            listGenes = utex.readFilesHuman(namefilezip, TCGAdb)
                            os.remove(namefilezip)
                            #Filter lists
                            for filter in listfilter:
                                listGenes = utex.applyFilter(listGenes, filter)
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
                            else:
                                listGenes = utex.readFilesHuman(namefilezip, TCGAdb)
                            os.remove(namefilezip)
                            # #Filter lists
                            for filter in listfilter:
                                listGenes = utex.applyFilter(listGenes, filter)
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
                    listGenes = utex.readFilesHuman(f, TCGAdb)
                    #Filter lists
                    for filter in listfilter:
                        listGenes = utex.applyFilter(listGenes, filter)
                    matrixGenes.append(listGenes)
        else:
            print('ERROR: FILES HAVE DIFFERENT EXTENSION. File need to have the same extension. All .csv or all .zip')
            sys.exit(-1)
    for f in listfilter:
        if f[0] == '-f':
            min_frel = float(f[1])
        elif f[0] == '-a':
            autoSaveImg = True
        elif f[0] == '-c':
            completeGraph = True
    return matrixGenes

#Read gene list name files
#Return list of files
def readParameters(input):
    listFiles = []
    listFilter = []
    i = 1
    while i < len(input):
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
            if filterParam[0] == input[i-1] and (filterParam[0] != '-a' and filterParam[0] != '-c'):
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
        elif re.search(r'^-fantom$', input[i]) or re.search(r'^-TCGA$', input[i]) or re.search(r'^-vitis$', input[i]):
            i += 1
        #if is else, return an error
        else:
            print("ERROR: "+input[i])
            i += 1
    #return an error if the list of file is empty
    if len(listFiles) < 2:
        print('ERROR: no files')
        utex.printInfo()

    #print('FILTER: '+str(listFilter))
    #print('FILES: '+str(listFiles))
    return (listFiles,listFilter)

#print file CSV with edges of graph
def printCSV(listCommonGenes):
    #create directory
    if not os.path.exists('commonGenesOutput'):
        os.mkdir('commonGenesOutput')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('commonGenesOutput/'+printTime)
    global nameDir
    nameDir = 'commonGenesOutput/'+printTime+'/'
    print('Creating directory: \''+nameDir+'\'', flush=True)
    for k in listCommonGenes:
        nameF = utex.buildNamefile(k)
        #create dir for each couple of genes
        nameDirGenes = nameDir+nameF+'/'
        os.mkdir(nameDirGenes)
        #Write file .csv
        f = open(nameDirGenes+nameF+'.csv', 'w')
        string = str(nameF)+'\n'
        f.write(string)
        for elem in k[1:]:
            string = str(elem[0])+','+str(elem[1])+','+str(elem[2])+','+str(elem[3])+'\n'
            f.write(string)
        f.close();

#Function to read the update name of human genes
def updateNameHuman():
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
    f.close()

#Function to read the update name of vitis genes
def updateNameVitis():
    #update the dictionary with name of vitis genes
    f = open('import_doc/NewAnnotVitisnet3.csv', 'r')
    text = f.readlines()
    listLineName = []
    i = 1
    while i < len(text):
        listLineName.append(text[i].split(','))
        i += 1
    for l in listLineName:
        if l[3] != '':
            if l[3] not in listBioNameUpdate.values():
                listBioNameUpdate[l[0].upper()] = l[3]
            else:
                listBioNameUpdate[l[0].upper()] = l[3]+'_'+l[0].upper()
        elif l[2] != '':
            if l[2] not in listBioNameUpdate.values():
                listBioNameUpdate[l[0].upper()] = l[2]
            else:
                listBioNameUpdate[l[0].upper()] = l[2]+'_'+l[0].upper()
        else:
            listBioNameUpdate[l[0].upper()] = l[0].upper()

#Main function
def main():
    if len(sys.argv) >= 2:
        if sys.argv[1] == '--help':
            utex.printInfo()
        else:
            #Read first parameter and set the correct boolean to true
            global TCGAdb
            if sys.argv[1] == '-fantom':
                updateNameHuman()
                TCGAdb = False
            elif sys.argv[1] == '-TCGA':
                updateNameHuman()
            elif sys.argv[1] == '-vitis':
                global vitis
                vitis = True
                TCGAdb = False
                updateNameVitis()
            else:
                print('ERROR: need to specify \'-fantom\' or \'-TCGA\'')
            #read parameters
            cmd = readParameters(sys.argv)
            #read each file .csv, read first line, check if is our gene
            #If YES -> go further, If NO -> check next file
            #When we have read all expansion of all gene, manage lists
            listCouple = utex.readFiles(cmd[0][0])
            #find common genes in the lists readed
            listCommonGenes = []
            isoformInEdge = []
            if TCGAdb or vitis:
                #read files if are TCGA or Vitis
                listFiles = readFilesGenes(cmd[0][1:], listCouple, cmd[1])
                if vitis:
                    #upload name if is Vitis
                    for l in listFiles:
                        l[0] = listBioNameUpdate[l[0]]
                        i = 1
                        while i < len(l):
                            try:
                                l[i] = (l[i][0], listBioNameUpdate[l[i][1]], l[i][2])
                            except:
                                pass
                            i += 1
                #find common genes
                listCommonGenes = utex.findCommonGenes(listCouple, listFiles)
            else:
                #read files if is Fantom
                listFiles = readFilesGenes(cmd[0][1:-1], listCouple, cmd[1])
                #upload name
                for l in listFiles:
                    l[0] = listBioNameUpdate[l[0]]
                    i = 1
                    while i < len(l):
                        l[i] = (l[i][0], listBioNameUpdate[l[i][1]], l[i][2])
                        i += 1
                #read which isoforms compone the edges
                isoformInEdge = utex.readFiles(cmd[0][-1])
                #find common genes
                listCommonGenes = utex.findCommonGenesFantom(listCouple, listFiles, isoformInEdge)

            #print CSV with genes share between every gene of LGN
            printCSV(listCommonGenes[0])
            #Draw the Venn diagram, Histogram
            utex.printNumberVenn(listCommonGenes, nameDir)
            graphic.printVenn(listCommonGenes[1], listCouple, nameDir)
            graphic.printHistogram(listCommonGenes[0], listFiles, nameDir, TCGAdb or vitis, isoformInEdge)

            pearsonComplete = []
            for l in listCommonGenes[0]:
                #Write genes to download and add the edges
                if completeGraph:
                    #write in a .txt the genes to download
                    listFileToRead = utex.nameAssociateGene(l, listBioNameUpdate, TCGAdb)
                    strToPrint = "To draw edges between associated genes please download lists in file \'"+nameDir+str(utex.buildNamefile(l))+"/listToDownload.txt"
                    f = open(nameDir+str(utex.buildNamefile(l))+'/listToDownload.txt', 'w')
                    for g in listFileToRead:
                        f.write(str(g) + ', ')
                    f.close()
                    print(strToPrint, flush=True)
                    #wait until they gave the path
                    print('Please write the path of file zip:', flush=True)
                    path = input('>>')
                    print('FILE ZIP: '+path, flush=True)
                    if vitis:
                        #read file and find edges between them
                        matrixGenes = vit.buildMatrixGenesVitis([('-f', min_frel)], [path])
                        print('Finding edges between genes...')
                        graphGenes = ut.manageDuplicates(ut.buildGraph(matrixGenes))
                        print('Process complete')
                        #save new edges in the list
                        nameF = utex.buildNamefile(l)
                        #create dir for each couple of genes
                        nameDirGenes = nameDir+nameF+'/'
                        #Write file .csv
                        f = open(nameDirGenes+nameF+'.csv', 'a')
                        i = 0
                        while i < len(graphGenes):
                            try:
                                rank = matrixGenes[[u[0] for u in matrixGenes].index(graphGenes[i][0])]
                                rank = rank[[v for (k,v,w,a,b) in rank[1:]].index(graphGenes[i][1])]
                            except:
                                rank = matrixGenes[[u[0] for u in matrixGenes].index(graphGenes[i][1])]
                                rank = rank[[v for (k,v,w,a,b) in rank[1:]].index(graphGenes[i][0])]
                            f.write(str(listBioNameUpdate[graphGenes[i][0]])+','+str(rank[0])+','+str(listBioNameUpdate[graphGenes[i][1]])+','+str(graphGenes[i][2])+'\n')
                            l.append((listBioNameUpdate[graphGenes[i][0]], rank[0], listBioNameUpdate[graphGenes[i][1]], graphGenes[i][2]))
                            i += 1
                        f.close()
                    else:
                        pass
                        #TODO: find a way to download the lists of human

                #Calculating pearson correlation for each edge
                print('Calculating Pearson correlation '+str(utex.buildNamefile(l))+'...')
                if vitis: #TODO: manage Pearson Correlation if expansion with GT-001
                    listForPearson = [((list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(a)]),(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(c)]),d) for (a,b,c,d) in l[1:]]
                    tmp = ut.pearsonCorrelation(listForPearson, 'vv_exprdata_2.csv')
                    pearson = [(listBioNameUpdate[u],listBioNameUpdate[v],p) for (u,v,p) in tmp]
                    pearsonComplete.append(pearson)
                else:
                    listForPearson = [(a,c,d) for (a,b,c,d) in l[1:]]
                    pearsonComplete.append(listForPearson)
                    #TODO: calculate pearson correlation for human
                print('Pearson Correlation done')

            #Draw graph
            graphic.printCommonGraph(listCommonGenes[0], pearsonComplete, 1-min_frel, nameDir, autoSaveImg)

#Calls the main() function.
if __name__ == '__main__':
    main()
