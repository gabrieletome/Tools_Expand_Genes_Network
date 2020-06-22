#########################################
# Draw graph of NESRA expansion         #
#########################################

import sys
import os
import datetime as dt
import lib.utilities_expansion as utex
import lib.utilities as ut
import lib.graphic as graphic

listBioNameUpdate = {}

#read file
def readFiles(filename):
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
    #listTuples = [listLine[0].split(' ')[3]]
    listTuples = [listLine[0].split('_')[0]] #TODEL
    i = 1
    while i < len(listLine):
        try:
            if listLine[i] != '' and listLine[i].split(',')[0] == 'GT-001': #TODEL GT-001
                #listTuples.append((int(listLine[i].split(',')[0]), listLine[i].split(',')[1].upper(), listLine[i].split(',')[2].upper(), float(listLine[i].split(',')[4]), bool(listLine[i].split(',')[5])))
                listTuples.append((listLine[i].split(',')[0], int(listLine[i].split(',')[1]), listLine[i].split(',')[2], float(listLine[i].split(',')[3]))) #TODEL
        except:
            pass
        i += 1
    return listTuples

#read file
def readFilesInt(filename):
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
    listTuples = []
    i = 1
    while i < len(listLine):
        try:
            if listLine[i] != '': #TODEL GT-001
                listTuples.append((int(listLine[i].split(',')[0]), listLine[i].split(',')[1].upper(), listLine[i].split(',')[2].upper(), float(listLine[i].split(',')[4]), (listLine[i].split(',')[5])))
                #listTuples.append((listLine[i].split(',')[0].upper(), int(listLine[i].split(',')[1]), listLine[i].split(',')[2].upper(), float(listLine[i].split(',')[3]))) #TODEL
        except:
            pass
        i += 1
    return listTuples


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

#Manage line with multiple genes
def manageBR(l, x1, x2):
    for u in l[1:]:
        if '<BR>' in u[x1]:
            names = u[x1].split('<BR>')
            tmp = list(u)
            tmp[x1] = names[0]
            l[l.index(u)] = tuple(tmp)
            for name in names[1:]:
                tmp = list(u)
                tmp[x1] = name
                l.append(tuple(tmp))
        if '<BR>' in u[x2]:
            names = u[x2].split('<BR>')
            tmp = list(u)
            tmp[x2] = names[0]
            l[l.index(u)] = tuple(tmp)
            for name in names[1:]:
                tmp = list(u)
                tmp[x2] = name
                l.append(tuple(tmp))
    return l

#Function for sorted list
def ord2(tuple):
    return tuple[0]

def main():
    if len(sys.argv) >= 2:
        #build matrix of genes
        #listTuple = manageBR(readFiles(sys.argv[-1]), 1, 2)
        listTuple = readFiles(sys.argv[-2]) #TODEL
        listTuple = [listTuple[0]]+sorted(listTuple[1:], key=ord2)
        updateNameVitis()
        #add filter
        filter = sys.argv[1:-2]
        if filter[0] == '-frel':
            #listTuple = [listTuple[0]]+[(a,b,c,d,e) for (a,b,c,d,e) in listTuple[1:] if d > float(filter[1])]
            listTuple = [listTuple[0]]+[(a,b,c,d) for (a,b,c,d) in listTuple[1:] if d > float(filter[1])] #TODEL
        elif filter[0] == '-rank':
            listTuple = [listTuple[0]]+[(a,b,c,d,e) for (a,b,c,d,e) in listTuple[1:] if a < int(filter[1])]
        #edges = [['MYB141','VvERF042','VvERF043','VvERF044']]+[(listBioNameUpdate[u[1]], u[0], listBioNameUpdate[u[2]], u[3]) for u in listTuple[1:]]
        edges = [['GT-001']]+[(u[0], u[1], u[2], u[3]) for u in listTuple[1:]]
        #read file interaction
        listTupleInt = [(u[0], listBioNameUpdate[u[1]], listBioNameUpdate[u[2]], u[3], u[4]) for u in manageBR(readFilesInt(sys.argv[-1]), 1, 2)] #TODEL
        edgesInt = [['MYB141','VvERF042','VvERF043','VvERF044']]+[(u[2], u[0], u[1], u[3]) for u in listTupleInt if u[4] == 'True' or u[1] in [k[2] for k in edges[1:]]]+[(u[1], u[0], u[2], u[3]) for u in listTupleInt if u[2] in [x[2] for x in edges[1:]]]
        edgesInt = [edgesInt[0]]+[(a,b,c,d) for (a,b,c,d) in edgesInt[1:] if d > float(filter[1])]
        print(edgesInt)
        
        pearsonComplete = []
        # listForPearson = [(b,c,d) for (a,b,c,d,e) in listTuple[1:]]
        # tmp = utex.manageBR(ut.pearsonCorrelation(listForPearson, 'vv_exprdata_2.csv'))
        #pearson = [(listBioNameUpdate[u], listBioNameUpdate[v],p) for (u,v,p) in manageBR(tmp, 0, 1)]
        pearson = [(u[0], u[2], 1) for u in edgesInt[1:]]
        pearsonComplete.append(pearson)
        #create directory
        if not os.path.exists('commonGenesOutput'):
            os.mkdir('commonGenesOutput')
        #write graph in a .txt file
        time = dt.datetime.now()
        printTime = time.strftime("%Y%m%d%H%M%S")
        os.mkdir('commonGenesOutput/'+printTime)
        os.mkdir('commonGenesOutput/'+printTime+'/'+utex.buildNamefile(edgesInt))
        os.mkdir('commonGenesOutput/'+printTime+'/'+utex.buildNamefile(edgesInt)+'/'+utex.buildNamefile(edgesInt))
        global nameDir
        nameDir = 'commonGenesOutput/'+printTime+'/'+utex.buildNamefile(edgesInt)+'/'
        print('Creating directory: \''+nameDir+'\'', flush=True)
        graphic.printCommonGraph([edgesInt], pearsonComplete, 1-round((edges[-1][3]), 1), nameDir, True, listBioNameUpdate)

#Calls the main() function.
if __name__ == '__main__':
    main()
