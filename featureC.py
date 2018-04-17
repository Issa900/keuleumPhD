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

#def read_count(input_file):
    #" count the number of reads >= 100bp on Y chromosome and autosomes"
    
    #with open (input_file, 'r') as f:
        #next(f)
        #MYdict={}
        #genenb=0
        #for line in islice(f, 2, None):
            #if line.startswith('Geneid'):
                #pass
            #genenb+=1
            #line = line.strip()
            #tmp=line.split()
            #Geneid=tmp[0]
            #Chr=tmp[1]
            #numFeature=len(Chr) # the number of read is considered here as the number of feauture detect on chromosome.
            #matchLen= int(tmp[5])
            #chr_name=Chr.split(';')[0]
           
            #if chr_name in MYdict and int(matchLen)>=100:
                ##print (matchLen)
                #MYdict[chr_name]+=numFeature
            #else:
                #MYdict[chr_name]=numFeature
       
        #print ('Number of Genes annotated =', genenb) # Can have this value by grep -c 'gene'feature_output.count
        ## Count the number of reads for the chromosome Y and autosomes.
        #ChrY_count = autosomal_count = 0
        #for key in MYdict:
            #if key=='Chr_Y_B' or key=='Chr_Y_A':
                #ChrY_count+=int(MYdict[key])
            #else:
                #autosomal_count+=int(MYdict[key])
        
        #print('number of reads on Y Chr =', ChrY_count)
        #print('number of reads on autosomes =', autosomal_count)
        #print('Total of mapped reads=', ChrY_count + autosomal_count)
        
        #return MYdict


    
def match_X_scaffold (list_X_specificGene, input_file):
    " Match genes on X chromosomes file found in the gtf file and extract the gene_id of the corresponding by storig it in the dict geneMatch "
    sdval=[]
    with open (input_file, 'r') as f, open (list_X_specificGene) as xGene:
        
        next(xGene)
                
        geneXlist=[]
        for line in xGene:
            line=line.strip()
            if line:
                gene_X=line.split()[0]
                geneXlist.append(gene_X)
        next(f)
        MYdict={}
        geneMatch={}
    
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
            sdval.append(int(coverage))
            
            
            
            if chr_name in MYdict:
                MYdict[chr_name]+=coverage
            else:
                MYdict[chr_name]=coverage
             
             
            #Count the coverage of reads matched with eacht X genes in the list
            if Gene_id in geneXlist:
                if Gene_id not in geneMatch:
                    geneMatch[Gene_id]=coverage
                else:
                    geneMatch[Gene_id]+=coverage
                
                
           # To count the number of features (exon), grouped into Meta-(gene)
            #if chr_name in MYdict and int(matchLen)>=100:
                #MYdict[chr_name]+=exonNum
            #else:
                #MYdict[chr_name]=exonNum
       
        
        ChrY_count = tmp_autosomal_count = ChrX_count = 0
        Y_reads=X_reads=A_reads=[]
        
        
        #Count total number of reads on Y and autosome only
        
        for key in MYdict:
            if key=='Chr_Y_B' or key=='Chr_Y_A':
                Y_reads.append(MYdict[key])
                ChrY_count+=int(MYdict[key])
                print(MYdict[key])
                Y_reads.append(int(MYdict[key]))
                
            else:
                tmp_autosomal_count+=int(MYdict[key])
                
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
        
        readLen=sorted(MYdict.items(), key=operator.itemgetter(1))
        
        minY=sorted(Y_reads)
        
        print(minY)
        print(Y_reads)
        
        
        sd=statistics.stdev(sdval)
        moy=sum(sdval)/len(sdval)
      
       
       
# Plot dist of coverage
        
        
       
        #sdval=np.array(sdval)
        #axes=range(len(sdval))
        #fig, ax = plt.subplots()
        #ax.plot(axes, sdval)
        
        plt.plot(range(len(sdval)), (sdval))
        plt.grid(True)
        plt.title(u"Distribution of scaffold coverage")
        plt.xlabel('Number of genes annotated', fontsize=12)
        plt.ylabel('Number of read mapped', fontsize=12)
        
        #plt.show()
        
        plt.savefig(count_file+'.disttribution.pdf')
                
        
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
