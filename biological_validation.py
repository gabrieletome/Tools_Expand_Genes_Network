import sys
import os
import datetime as dt
import rpy2.robjects as ro
import pandas as pd
import lib.diffexp_go_analysis as topGO

#Create saving directory
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

#Print txt with all information of topGO
def print_output_topGO(results_table, results, nameDir, ontology):
    #Write file .txt
    f = open(nameDir+'validation_topGO_'+ontology+'.txt', 'w')
    f.write(str(results)+'\n\n')
    pd_dt = ro.conversion.rpy2py(results_table)
    f.write(str(pd_dt))
    f.close();
    print('Create: '+nameDir+'validation_topGO_'+ontology+'.txt', flush=True)

def main():
    if len(sys.argv) >= 2:
        if sys.argv[1] == '-topGO':
            #read parameters
            nameDir = createSavingDir()
            #results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt')
            #Run BP ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], sys.argv[3], nameDir, 'BP')
            print_output_topGO(results_table, results, nameDir, 'BP')
            #Run MF ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], sys.argv[3], nameDir, 'MF')
            print_output_topGO(results_table, results, nameDir, 'MF')
    else:
        print('ERROR: wrong nuomber of parameters')

#Calls the main() function.
if __name__ == '__main__':
    main()
