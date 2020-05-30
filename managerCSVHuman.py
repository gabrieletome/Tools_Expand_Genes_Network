import sys
import os
import re

def readFiles(fileHGNC):
    HGNC = readfileHGNC(fileHGNC)
    # list_text_files = []
    # for name in filename:
    #     list_text_files.append(readList(name))
    listCouple = []
    # for l in list_text_files:
    #     for g in l:
    #         g = (g.split('\"'))[1]
    #         print(g)
    #         if g in HGNC.keys():
    #             listCouple.append((g, HGNC[g]))
    for g in HGNC.keys():
        listCouple.append((g, HGNC[g]))
    return listCouple


#Read file .csv
#Return list of tuples with position [0] name of the gene, in pos 2+ all gene associated
def readList(filename):
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
    listTuples = []
    i = 0
    while i < len(listCells):
        try:
            k = 1
            while k < len(listCells[i]):
                if listCells[i][k] != '':
                    listTuples.append(listCells[i][k])
                k += 1
        except:
            pass
        i += 1
    return listTuples


#Read file .csv
#Return list of tuples with position [0] name of the gene, in pos 2+ all gene associated
def readfileHGNC(filename):
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
    dictTuples = {}
    i = 1
    while i < len(listCells):
        try:
            dictTuples[listCells[i][0]] = listCells[i][6]
        except:
            pass
        i += 1
    return dictTuples

#Read gene list name files
#Return list of files
def readParameters(input):
    listFiles = []
    i = 1
    while i < len(input):
        #if find pattern '^-files$', the following parameters are files
        if re.search(r'^-files$', input[i]):
            i += 1
            while i < len(input):
                listFiles.append(input[i])
                i += 1
        #if is else, return an error
        else:
            print("ERROR: "+input[i])
            i += 1
    #return an error if the list of file is empty
    if len(listFiles) == 0:
        print('ERROR: no files')

    #print('FILTER: '+str(listFilter))
    #print('FILES: '+str(listFiles))
    return listFiles

#print file CSV with edges of graph
def printCSV(nameFile, edges):
    f = open(nameFile+'.csv', 'w')
    print(len(edges))
    for k in edges:
        string = k[0]+','+k[1]+'\n'
        f.write(string)
    f.close();
    print('Create: \''+nameFile+'.csv\'', flush=True)

def main():
    if len(sys.argv) >= 1:
        #read parameters
        #cmd = readParameters(sys.argv)
        fileHGNC = 'Lists_Human/hgnc_filtered_anno.csv'
        #build matrix of genes
        listCouple = readFiles(fileHGNC)
        printCSV('Lists_Human/couple_name_gene', listCouple)

#Calls the main() function.
if __name__ == '__main__':
    main()
