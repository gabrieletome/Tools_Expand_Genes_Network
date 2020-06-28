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

#Create fasta file
def createFasta(listGenes, completeFasta, dir):
    with open(listGenes) as in_handle:
        genes_list = topGO.parse_input_csv(in_handle)
    f_complete = open(completeFasta, 'r')
    completeF = f_complete.read().split('\n')
    dict_FASTA = {}
    i = 0
    while i < len(completeF):
        if completeF[i][1:] in genes_list:
            dict_FASTA[completeF[i]] = completeF[i+1]
        i += 2
    f_complete.close()
    namefileFASTA = nameDir+'fasta_'+listGenes.split('/')[-1][:-4]+'.fasta'
    print(namefileFASTA)
    f = open(namefileFASTA, 'w')
    for k in dict_FASTA.keys():
        f.write(k+'\n'+dict_FASTA[k]+'\n')
    f.close()
    return namefileFASTA

def main():
    if len(sys.argv) >= 2:
        if sys.argv[1] == '-topGO':
            #create saving directory
            nameDir = createSavingDir()
            #results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt')
            #Run BP ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], sys.argv[3], nameDir, 'BP')
            print_output_topGO(results_table, results, nameDir, 'BP')
            #Run MF ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], sys.argv[3], nameDir, 'MF')
            print_output_topGO(results_table, results, nameDir, 'MF')
        elif sys.argv[1] == '-dreme':
            #create saving directory
            nameDir = createSavingDir()
            #create fasta file from list of genes
            fastaFile = createFasta(sys.argv[2], sys.argv[3], nameDir)
            #execute DREME analisys
            commandDREME = 'dreme -o '+nameDir+'result_DREME -p '+fastaFile+' -n '+sys.argv[3]
            print('Execute: '+commandDREME)
            os.system(commandDREME)
    else:
        print('ERROR: wrong nuomber of parameters')

#Calls the main() function.
if __name__ == '__main__':
    main()
