
#!/usr/bin/env python
# -*- coding: utf-8 -*-






###++++++++++++++++++++++++++++++++++++++++

## This program takes as input: .fasta and a table .csv  

###++++++++++++++++++++++++++++++++++++++++

import sys, os, re
##"import" statement very IMPORTant
# Import statements are executed in two steps: (1) find a module, and initialize it if necessary; (2) define a name or names in the local namespace (of the scope where the "import" statement occurs). The statement comes in two forms differing on whether it uses the "from" keyword. The first form (without "from") repeats these steps for each identifier in the list. The form with "from" performs step (1) once, and then performs step (2) repeatedly. To understand how step (1) occurs, one must first understand how Python handles hierarchical naming of modules. To help organize modules and provide a hierarchy in naming, Python has a concept of packages. A package can contain other packages and modules while modules cannot contain other modules or packages. From a file system perspective, packages are directories and modules are files. Once the name of the module is known (unless otherwise specified, the term "module" will refer to both packages and modules), searching for the module or package can begin. The first place checked is "sys.modules", the cache of all modules that have been imported previously. If the module is found there then it is used in step (2) of import. 


## For more details of any other statement in python, visit the main website of PYTHON (https://www.python.org/doc/)

# An other way is to open python from your terminal by typing:

# $ python
# $ help()
# $ import (or any other statementin python).
import csv
import time
# So after importing all module that I will use in my scipt like "sys", "os", "re", etc.... I create a function by the command (def), see next #####comment### 

def fasta_header(fasta_in):  
    # So here is my function called "fasta_hearder", function can take argurments as input (here: "fasta_in" or output file, or anything, it depend on what you want put into it
    print 'Data processing \n Please wait...' 
    dict_proteome = {} # create a dictionnary called dict_proteome to store all sequences that I will return at the end of the function
    with open (fasta_in) as f: # Here Erica, I open the input and run a for command to differenciate the fasta_headers and the sequences in the fasta file
        for line in f: # And usally when I open a for loop, Magie take place :p. "for" means "FOR EACH LINE IN MY FILE, DO something...
            line = line.strip() # Here I delete blank space by stripping them  
            if line.startswith('>'): # Here you identify your header by a if statement
                head= line[1:].split()[0] # extract the position of the key for your dictionnary
                dict_proteome[head]= "" # create the key in your dict with empty space !
            else: # hahahaah WHAT ELSE :p
                #sequence= line
                dict_proteome[head]+=line # If you line is not a header, so it's a sequence so, and add it here as a value in the dict
            #print(dict_proteome)
         
         #Once the iteration done and all head ans seq in my dict, I 
            
    return dict_proteome # At  I retuen my dict (dict_proteome), you remember that we create one above, so now we fiil it in, we want to have it back, so we use return and we stored {header as key and sequence as value}.


def table_xl(table): # Here I create an other fonction but the input is not the same, here we have a csv file, so the way to extract element is diffrent to fasta file

    with open (table, 'r') as t: #open it and read it !
        reader = csv.reader(t) # Here happens the Magie of the importing module as CSV, because you have arguments which are already built in as csv.reader, it's why there are called "built-in function" and blah blah blah 
        r1=next(reader)
        for elem in r1:
            l = elem.split('.')[0]
            sp_name = l.split('_')[-1]
            print(sp_name)
        #t.next()
        cluster_dict={} # Oups a new dictionnary here, now you should keep how to create in pyhton :p
        nb=0 # Here I initiate a line counter
        for line in t: # for each line in the ...? I extract specie name, gene and Alg
            line=line.strip() 
            Specie = line.split(',')[0]
            Gene = line.split(',')[1]
            Alg = line.split(',')[2]
            nb+=1
            
            ###THESE MANY "for loops" could be easy avoid if you checked your data carefully and told that there were not uniforms. Next time, just bring a test file, if fasta dont bring nothing juste write down CLEARLY, what do you want ?we have a '*' in many header

            Elytraria = line.split(',')[3] # So basically what I did with each for loop is just to take different elements (3 - 10) in the list line and to remane it if it's empty :( So if you read until here and you wanna read the next comments.
            for name in Elytraria.split(';'):
                if name != '*':
                    new_Ely=name+'|'+'Elytraria'

            Mendoncia = line.split(',')[4]
            for name in Mendoncia.split(';'):
                if name != '*':
                    new_Men=name+'|'+'Mendoncia'

            Thunbergia = line.split(',')[5]
            for name in Thunbergia.split(';'):
                if name != '*':
                    new_Thu=name+'|'+'Thunbergia'

            Avicennia = line.split(',')[6]
            for name in Avicennia.split(';'):
                if name != '*':
                    new_Avi=name+'|'+'Avicennia'


            Andrographis = line.split(',')[7]
            for name in Andrographis.split(';'):
                if name != '*':
                    new_And=name+'|'+'Andrographis'

            Pachystachys = line.split(',')[8]
            for name in Pachystachys.split(';'):
                if name != '*':
                    new_Pac=name+'|'+'Pachystachys'

            Acanthus = line.split(',')[9]
            for name in Acanthus.split(';'):
                if name != '*':
                    new_Aca=name+'|'+'Acanthus'

            Aphelandra = line.split(',')[10]
            for name in Aphelandra.split(';'):
                if name != '*':
                    new_Aph=name+'|'+'Aphelandra'

            cluster= str(new_Ely) + "\n" + str(new_Men) + "\n"+ str(new_Thu) + '\n'+str(new_Avi)+'\n' + str(new_And)+'\n' +str(new_Pac)+'\n'+str(new_Aca)+'\n'+str(new_Aph)+'\n'
            cluster_dict[nb] = cluster
            # if nb == 30:
            #     print(cluster)
    return cluster_dict

def writing_sequence(dict_proteome, cluster_dict, file_out):

    new_dico={}
    num = 0
    for number in cluster_dict:
        num +=1
        cluster = cluster_dict[number]
        with open(file_out + str(num)+'.fasta', 'w') as outf:
            for full_name in cluster.split('\n'):
                name = (full_name.split('|')[0] )
                if name in dict_proteome:
                    seq= dict_proteome[name]
                    outf.write('> '+ str(full_name) +'\n')
                    outf.write(seq +'\n')

file_in =sys.argv[1] # fasta file "with "*" in some header 'sequence name'
table_in =sys.argv[2] # cvs file
file_out= '/home/issa900/Documents/prog/results/cluster_'
x= fasta_header(file_in)
y = table_xl(table_in)
writing_sequence(x, y, file_out)
print ('Time consuming', time.time())


## Run maftt auto
# for i in `ls results`; do mafft --maxiterate 1000 --globalpair $i > aln_dir/$i.aln

## Split fasta file by size (from Erica) 
# https://biowize.wordpress.com/2012/01/20/splitting-fasta-files-by-size/
