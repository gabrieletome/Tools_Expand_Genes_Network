#!/usr/bin/env python
#################################################
# Code from: https://github.com/chapmanb/bcbb   #
# Library to use R package topGO in Python      #
# Adapted to work in Python3                    #
#################################################
"""Provide topGO analysis of overrepresented GO annotation terms in a dataset.

Usage:
    stats_go_analysis.py <input CVS> <gene to GO file>
"""
from __future__ import with_statement
import sys
import csv
import collections
import os
import rpy2.robjects as robjects

def main(input_csv, gene_to_go_file):
    topGO_analysis(input_csv, gene_to_go_file, '')

def topGO_analysis(input_csv, gene_to_go_file, nameDir):
    gene_pval = 1e-2
    #go_pval = 0.2
    go_pval = 1.01
    go_term_type = "MF"
    #go_term_type = "BP"
    #topgo_method = "classic" # choice of classic, elim, weight
    topgo_method = 'weight01'

    with open(input_csv) as in_handle:
        genes_w_pvals = parse_input_csv(in_handle)
    with open(gene_to_go_file) as in_handle:
        gene_to_go, go_to_gene = parse_go_map_file(in_handle, genes_w_pvals)
    if len(gene_to_go) == 0:
        raise ValueError("No GO terms match to input genes. "
              "Check that the identifiers between the input and GO file match.")
    go_terms, results_table, results = run_topGO(genes_w_pvals, gene_to_go, go_term_type, gene_pval, go_pval, topgo_method, nameDir)
    genes, strValue = print_go_info(go_terms, go_term_type, go_to_gene)
    return genes, strValue, results_table, results

def print_go_info(go_terms, go_term_type, go_to_gene):
    for final_pval, go_id, go_term in go_terms:
        genes = []
        for check_go in [go_id] + get_go_children(go_id, go_term_type):
            genes.extend(go_to_gene.get(check_go, []))
        genes = sorted(list(set(genes)))
        strValue = "-> %s (%s) : %0.4f" % (go_id, go_term, final_pval)
        print(strValue)
        for g in genes:
            print(g)
        return genes, strValue

def get_go_children(go_term, go_term_type):
    """Retrieve all more specific GO children from a starting GO term.
    """
    robjects.r('''
        library(GO.db)
    ''')
    child_map = robjects.r["GO%sCHILDREN" % (go_term_type)]
    children = []
    to_check = [go_term]
    while len(to_check) > 0:
        new_children = []
        for check_term in to_check:
            new_children.extend(list(robjects.r.get(check_term, child_map)))
        new_children = list(set([c for c in new_children if type(c) == type('str')]))
        children.extend(new_children)
        to_check = new_children
    children = list(set(children))
    return children

def _dict_to_namedvector(init_dict):
    """Call R to create a named vector from an input dictionary.
    """
    return robjects.r.c(**init_dict)

def run_topGO(gene_vals, gene_to_go, go_term_type, gene_pval, go_pval, topgo_method, nameDir):
    """Run topGO, returning a list of pvalues and terms of interest.
    """
    # run topGO with our GO and gene information
    robjects.r('''
        library(topGO)
    ''')
    robjects.r('''
        topDiffGenes = function(allScore) {
          return (allScore < %s)
        }
    ''' % gene_pval)
    params = {"ontology" : go_term_type,
              "annot" : robjects.r["annFUN.gene2GO"],
              "geneSelectionFun" : robjects.r["topDiffGenes"],
              "allGenes" : _dict_to_namedvector(gene_vals),
              "gene2GO" : _dict_to_namedvector(gene_to_go)
              }
    go_data = robjects.r.new("topGOdata", **params)
    results = robjects.r.runTest(go_data, algorithm=topgo_method, statistic="fisher")
    scores = robjects.r.score(results)
    num_summarize = min(50, len(scores.names))
    # extract term names from the topGO summary dataframe
    #results_table = robjects.r.GenTable(go_data, elimFisher=results, orderBy="elimFisher", topNodes=num_summarize)
    results_table = robjects.r.GenTable(go_data, classicFisher=results, topNodes=num_summarize)

    robjects.r.showSigOfNodes(go_data, scores, firstSigNodes = 5, useInfo='all')
    os.remove('Rplots.pdf')
    paramPrint = {'firstSigNodes': 5, 'fn.prefix': nameDir+"graph", 'useInfo': "all", 'pdfSW': True}
    robjects.r.printGraph(go_data, results, **paramPrint)

    GO_ID_INDEX = 0
    TERM_INDEX = 1
    ids_to_terms = dict()
    for index, go_id in enumerate(results_table[GO_ID_INDEX]):
        ids_to_terms[go_id] = results_table[TERM_INDEX][index]
    go_terms = []
    # convert the scores and results information info terms to return
    for index, item in enumerate(scores):
        if item < go_pval:
            go_id = scores.names[index]
            go_terms.append((item, go_id, ids_to_terms.get(go_id, "")))
    go_terms.sort()
    return go_terms, results_table, results

def parse_go_map_file(in_handle, genes_w_pvals):
    gene_list = genes_w_pvals.keys()
    gene_to_go = collections.defaultdict(list)
    go_to_gene = collections.defaultdict(list)
    for line in in_handle:
        parts = line.split("\t")
        gene_id = parts[0]
        go_id = parts[1].strip()
        if gene_id in gene_list:
            gene_to_go[gene_id].append([u for u in go_id.split(',') if u != ''])
            for k in go_id.split(','):
                if k != '':
                    go_to_gene[k].append(gene_id)
                # gene_to_go[gene_id].append(go_id)
                # go_to_gene[go_id].append(gene_id)
    return dict(gene_to_go), dict(go_to_gene)

def parse_input_csv(in_handle):
    reader = csv.reader(in_handle)
    reader.__next__() # header
    all_genes = dict()
    for (gene_name, _, _, pval) in reader:
        all_genes[gene_name] = float(pval)
    return all_genes

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print( __doc__)
        sys.exit()
    main(sys.argv[1], sys.argv[2])
