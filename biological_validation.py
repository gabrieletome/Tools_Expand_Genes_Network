import sys
import os
import datetime as dt
import rpy2.robjects as ro
import pandas as pd
import lib.diffexp_go_analysis as topGO

def createSavingDir():
    #create directory
    if not os.path.exists('outputBiologicalValidation'):
        os.mkdir('outputBiologicalValidation')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('outputBiologicalValidation/'+printTime)
    global nameDir
    nameDir = 'outputBiologicalValidation/'+printTime+'/'
    print('Creating directory: \''+nameDir+'\'', flush=True)
    return nameDir

def print_output_topGO(genes, strValue, results_table, results, nameDir):

    #Write file .csv
    f = open(nameDir+'validation_topGO.txt', 'w')
    f.write(strValue+'\n\n')
    string_tmp = 'LIST GENES: '
    for g in genes[:-1]:
        string_tmp = string_tmp+str(g)+', '
    string_tmp = string_tmp+str(genes[-1])
    f.write(string_tmp+'\n')

    f.write(str(results)+'\n\n')
    pd_dt = ro.conversion.rpy2py(results_table)
    f.write(str(pd_dt))
    f.close();

def main():
    if len(sys.argv) >= 2:
        if sys.argv[1] == '-topGO':
            #read parameters
            nameDir = createSavingDir()
            #genes, strValue, results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt')
            genes, strValue, results_table, results = topGO.topGO_analysis(sys.argv[2], sys.argv[3], nameDir)
            print_output_topGO(genes, strValue, results_table, results, nameDir)
    else:
        print('ERROR: wrong nuomber of parameters')

#Calls the main() function.
if __name__ == '__main__':
    main()
