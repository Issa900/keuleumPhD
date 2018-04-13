#!/usr/bin/env python
# -*- coding: utf-8 -*-


#""Thisfuntion takes as  input a ""
import os, sys, argparse, re
# import matplotlib
# import matplotlib.pyplot as plt
# import numpy as np
# import itertools

def read_count(input_file):
    " count the number of reads on Y chromosome and autosome"
    
    with open (input_file, 'r') as f:
        next(f)
        MYdict={}
        genenb=0
        for line in f:
            genenb+=1
            line = line.strip()
            tmp=line.split()
            Geneid=tmp[0]
            Chr=tmp[1]
            numFeature=len(Chr) # the number of read is considered here as the number of feauture detect on chromosome.
            chr_name=Chr.split(';')[0]
           
            if chr_name in MYdict:
                MYdict[chr_name]+=numFeature
            else:
                MYdict[chr_name]=numFeature
       
        print ('Number of Genes annotated =', genenb) # Can have this value by grep -c 'gene'feature_output.count
        
        # Count the number of reads for the chromosome Y and autosomes.
        ChrY_count = autosomal_count = 0
        for key in MYdict:
            if key=='Chr_Y_B' or key=='Chr_Y_A':
                ChrY_count+=int(MYdict[key])
            else:
                autosomal_count+=int(MYdict[key])
        #print('number of reads on Y Chr =', ChrY_count)
        #print('number of reads on autosomes =', autosomal_count)
        #print('Total of mapped reads=', ChrY_count + autosomal_count)
        
        return MYdict

def match_X_scaffold (list_X_specificGene, gtf_file):
    " Match genes on X chromosomes file found in the gtf file and extract the gene_id of the corresponding by storig it in the dict geneMatch "
    with open (list_X_specificGene) as xGene, open (gtf_file) as gtf:
        
        next(xGene)
        
        geneXlist=[]
        dico_gene={}
        geneMatch={}
        
        for line in xGene:
            line=line.strip()
            if line:
                gene_X=line.split()[0]
                geneXlist.append(gene_X)
        print(len('Number of X chr genes ='geneXlist))
        
        for gene in gtf:
            tmp=gene.split()
            gene_id=tmp[11].strip('";')
            gene_name=tmp[13].strip('";')
            
            if gene_id in dico_gene:
                dico_gene[gene_id].append(genename)
            else:
                dico_gene[gene_id]=gene_name
            
        for x in dico_gene:
            val=dico_gene[x]
            #for x in gene
            if val in geneXlist:
                print(val)
                if x in geneMatch:
                    geneMatch[x].append(val)
                else:
                    geneMatch[x]=val 
        
        return geneMatch
                
            
                

    
# input_file='/home/issa/Documents/prog/LCC.bam.txt'
count_file=sys.argv[1] # Output file of the program FeatureCount(.txt)
gtf=sys.argv[2] # GTF file (Dowload the GFF3 file from Marchantia.info and convert into GTF by gffread pachage
x=sys.argv[3] # List of genes identified on the X chromosome

f1=read_count(count_file)

match_X_scaffold(x, gtf)
