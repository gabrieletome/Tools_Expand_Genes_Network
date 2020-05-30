import sys
import re

def readFiles(file):
    print('Open file: '+file[1], flush=True)
    try:
        f = open(file[1], 'rU')
        text = f.read()
        f.close()
    except:
        #except if file does not exist
        print('ERROR: FILE NOT FOUND. File \''+file[1]+'\' does not exist')
        sys.exit(-1)
    #split rows
    listLine = text.split('\n')
    listTuples = []
    i = 0
    while i < len(listLine):
        try:
            if listLine[i] != '':
                listTuples.append(listLine[i])
        except:
            pass
        i += 1
    lines = readfile(file[0], listTuples)
    return lines

#Read file .csv
#Return list of tuples with position [0] name of the gene, in pos 2+ all gene associated
def readfile(filename, listGene):
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
        if listCells[i][0] in listGene:
            listTmp = [listCells[i][0]]
            try:
                k = 1
                while k < len(listCells[i]):
                    if listCells[i][k] != '':
                        listTmp.append(listCells[i][k])
                    k += 1
            except:
                pass
            listTuples.append(listTmp)
        i += 1
    return listTuples

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

def countZero(listLines):
    listCouple = []
    for l in listLines:
        count = 0
        i = 1
        while i < len(l):
            if float(l[i]) == 0:
                count += 1
            i += 1
        if len(listCouple) == 0:
            listCouple.append(('Total column:', i))
        listCouple.append((l[0], count))
    return listCouple

#print file CSV with edges of graph
def printCSV(nameFile, edges):
    f = open(nameFile+'.csv', 'w')
    print(len(edges))
    for k in edges:
        string = k[0]+','+str(k[1])+'\n'
        f.write(string)
    f.close();
    print('Create: \''+nameFile+'.csv\'', flush=True)

def main():
    if len(sys.argv) >= 2:
        #read parameters
        cmd = readParameters(sys.argv)
        #build matrix of genes
        listCouple = countZero(readFiles(cmd))
        printCSV('Lists_Human/zero_matrix_pearson', listCouple)

#Calls the main() function.
if __name__ == '__main__':
    main()
