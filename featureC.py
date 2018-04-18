#!/usr/bin/env python
# -*- coding: utf-8 -*-


# 
#*** takes as  input files i) output file of FeatureCount, and ii) a list of specif genes_id
import os, sys, argparse, re
import itertools
from itertools import islice
import operator
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import statistics, math

    
def match_X_scaffold (list_X_specificGene, input_file):
    " Match genes on X chromosomes file found in the gtf file and extract the gene_id of the corresponding by storig it in the dict geneMatch "
    
    with open (input_file, 'r') as f, open (list_X_specificGene) as xGene:
        
        next(xGene)
                
        geneXlist=[]
        for line in xGene:
            line=line.strip()
            if line:
                gene_X=line.split()[0]
                geneXlist.append(gene_X)
       
        all_reads=[]
        chr_dict={}
        gene_dict={}
        geneMatch={}
        
        Y_reads=[]
        X_reads=[]
        A_reads=[]
        
        #nb=0
        genenb=0
        for line in islice(f, 2, None):
            if line.startswith('Geneid'):
                pass
            
            genenb+=1
            line = line.strip()
            tmp=line.split()
            
            Gene_id=tmp[0]
            
            
            Chr=tmp[1]
            chr_name=Chr.split(';')[0]
            exonNum=len(Chr.split(';')) 
            
            matchLen=int(tmp[5])
            
            coverage=int(tmp[-1])+1
            all_reads.append(int(coverage))
            
            if Gene_id in gene_dict:
                gene_dict[Gene_id]+=coverage
            else:
                gene_dict[Gene_id]=coverage
            
            if chr_name in chr_dict:
                chr_dict[chr_name]+=coverage
            else:
                chr_dict[chr_name]=coverage
             
             
            #Count the coverage of reads matched with eacht X genes in the list
            
            if Gene_id in geneXlist:
                #nb+=1
                X_reads.append(coverage)
            
                if Gene_id not in geneMatch:
                    geneMatch[Gene_id]=coverage
                else:
                    geneMatch[Gene_id]+=coverage
                    
            elif chr_name=='Chr_Y_A' or chr_name=='Chr_Y_B':
                Y_reads.append(coverage)
            else:
                A_reads.append(coverage)
           # To count the number of features (exons), grouped into Meta-(gene)
            #if chr_name in MYdict and int(matchLen)>=100:
                #MYdict[chr_name]+=exonNum
            #else:
                #MYdict[chr_name]=exonNum
       
        
        ChrY_count = tmp_autosomal_count = ChrX_count = 0
        
        
        #Count total number of reads on Y and autosome only
        
        for key in chr_dict:
            if key=='Chr_Y_B' or key=='Chr_Y_A':
                #Y_reads.append(int(gene_dict[key]))
                ChrY_count+=int(chr_dict[key])
                print(chr_dict[key])
                
            else:
                tmp_autosomal_count+=int(chr_dict[key])
                
        for gene_name in geneMatch:
            i=int(geneMatch[gene_name])
            ChrX_count+=i
            
        autosomal_count = int(tmp_autosomal_count) - int(ChrX_count)
        
        print('#Number of reads on:')
        print('Y Chromosome =', ChrY_count)
        print('X Chromosome =', ChrX_count)
        print('Autosomes =', autosomal_count)
        print('Total =', ChrY_count + autosomal_count+ChrX_count)
        
        print ('Number of Genes annotated =', genenb) # Can have this value by grep -c 'gene'feature_output.count
        # Count the number of reads for the chromosome Y and autosomes.
        
        readLen=sorted(gene_dict.items(), key=operator.itemgetter(1))
        #print(readLen[0])
        
        #Min, Max and sd of reads mapped
        
        minA=sorted(A_reads)[0]
        maxA=sorted(A_reads)[-1]
        sd_A=statistics.stdev(A_reads)
        moy_A=sum(A_reads)/len(A_reads)
        
        print(len(A_reads))
      
        minY=sorted(Y_reads)[0]
        maxY=sorted(Y_reads)[-1]
        sd_Y=statistics.stdev(Y_reads)
        moy_Y=sum(Y_reads)/len(Y_reads)
        print(len(Y_reads))
        
        minX=sorted(X_reads)[0]
        maxX=sorted(X_reads)[-1]
        sd_X=statistics.stdev(X_reads)
        moy_X=sum(X_reads)/len(X_reads)
        print(len(X_reads))
        
       
# Plot dist of coverage
        
        
       
        #sdval=np.array(sdval)
        #axes=range(len(sdval))
        #fig, ax = plt.subplots()
        #ax.plot(axes, sdval)
        
        plt.plot(range(len(all_reads)), (all_reads))
        plt.grid(True)
        plt.title(u"Distribution of scaffold coverage")
        plt.xlabel('Number of genes annotated', fontsize=12)
        plt.ylabel('Number of read mapped', fontsize=12)
        
        #plt.show()
        
        #plt.savefig(count_file+'.disttribution.pdf')
                
        
    #read gtf_file
    #with open (gtf_file) as gtf:
        #for gene in gtf:
            #tmp=gene.split()
            #gene_id=tmp[11].strip('";')
            #gene_name=str(tmp[13].strip('";'))
            
            #if gene_id in dico_gene.keys():
                #list(dico_gene[gene_id]).append(gene_name)
            
            #else:
                #dico_gene[gene_id]=gene_name
        
        
        #for x in dico_gene:
            #val=dico_gene[x]
            #for gene in geneXlist:
                #if str(val)==str(gene): #[str(i) for i in geneXlist]:
                    #geneMatch[x]=gene
        ##print((geneMatch))
        

    
                

count_file=sys.argv[1] # Output file of the program FeatureCount(.txt)
#gtf=sys.argv[2] # GTF file (Dowload the GFF3 file from Marchantia.info and convert into GTF by gffread pachage
x=sys.argv[2] # List of genes identified on the X chromosome

#f1=read_count(count_file)
ID_Xgenes=match_X_scaffold(x, count_file)
