#!/usr/bin/env python
# -*- coding: utf-8 -*-


# 
#*** takes as  input files i) output file of FeatureCount, and ii) a list of specific genes_id
import os, sys, argparse, re
import itertools
from itertools import islice
import operator
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import seaborn as sns
import statistics, math
import scipy.stats as stats
    
def match_X_scaffold (list_X_specificGene, input_file):
    " Match genes on X chromosomes file found in the gtf file and extract the gene_id of the corresponding by storig it in the dict geneMatch "
    print('\n Analyses of the number of reads mapped on sex chromosome versus autosomes \n \n')
    
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
                
            else:
                tmp_autosomal_count+=int(chr_dict[key])
                
        for gene_name in geneMatch:
            i=int(geneMatch[gene_name])
            ChrX_count+=i
            
        autosomal_count = int(tmp_autosomal_count) - int(ChrX_count)
        
        
        print ('# Total number of genes annotated =\n', genenb)
        
        print('\n# Total number of reads mapped on:', '\n Y Chromosome =', ChrY_count, '\n X Chromosome =', ChrX_count, '\n Autosomes =', autosomal_count, '\nTotal =', ChrY_count + autosomal_count+ChrX_count, '\n\n')
        
        
         # Can have this value by grep -c 'gene'feature_output.count
        # Count the number of reads for the chromosome Y and autosomes.
        
        #readLen=sorted(gene_dict.items(), key=operator.itemgetter(1))
        #print(readLen[0])
        
        #Min, Max and sdt of reads mapped
        print('# Difference in the number of reads mapped\n')
        
    
        minY=sorted(Y_reads)[0]
        maxY=sorted(Y_reads)[-1]
        sd_Y=statistics.stdev(Y_reads)
        moy_Y=sum(Y_reads)/len(Y_reads)
        print('# Y Chromosome \n Minimum =', minY, '\n Maximum =', maxY, '\n SD =', sd_Y)
        
        minX=sorted(X_reads)[0]
        maxX=sorted(X_reads)[-1]
        sd_X=statistics.stdev(X_reads)
        moy_X=sum(X_reads)/len(X_reads)
        print('# X Chromosome \n Minimum =', minX, '\n Maximum =', maxX, '\n SD =', sd_X)
        
        minA=sorted(A_reads)[0]
        maxA=sorted(A_reads)[-1]
        sd_A=statistics.stdev(A_reads)
        moy_A=sum(A_reads)/len(A_reads)
        print('# Autosomes \n Minimum =', minA, '\n Maximum =', maxA, '\n SD =', sd_A)
       
# Plot dist of coverage
        
        
       
        
        #fig, ax = plt.subplots()
        #ax.plot(axes, sdval)
        
        all_reads.sort()
        
        ##Normalization
        
        #mean_r = np.mean(all_reads)
        #std_r=np.std(all_reads)
        #s = stats.norm.pdf(all_reads, mean_r, std_r) 
        #plt.plot(all_reads, s)
        
        plt.plot(range(len(all_reads)), (all_reads), color='red', label= 'sorted data')
        
        plt.grid(True)
        plt.title(u"Distribution of scaffold coverage")
        plt.xlabel('Number of genes annotated', fontsize=12)
        plt.ylabel('Number of read mapped', fontsize=12)
        
        #plt.show()
        
        plt.savefig(count_file+'.sorted.pdf')
                
        
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
#gtf=sys.argv[2] # GTF file (Dowload the GFF3 file from Marchantia.info and convert into GTF by gffread package
x=sys.argv[2] # List of genes identified on the X chromosome

#f1=read_count(count_file)
ID_Xgenes=match_X_scaffold(x, count_file)
