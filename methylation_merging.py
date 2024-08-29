# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 22:01:14 2020

@author: mneme
"""


import os
import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--fileList", help="List of methylation files to combine. Should have full paths.")
parser.add_argument("-c", "--cytosine_positions", help="Bed file of positions of cytosines of interest (CA, CG, etc.).")
parser.add_argument("-of", "--out_file", help="Output file of merged methylation")
parser.add_argument("-os", "--out_samples", help="Output file showing methylation of samples.")
args = parser.parse_args()

filenames=[]
samples = []
with open(args.fileList, "r") as fileList:
    for line in fileList:
        data = line.strip()
        filenames.append(data)
        samples.append(data.split("/")[-1])
        
files = [open(file,"r") for file in filenames]
read_go = [True for n in filenames]
read_store = ["" for n in filenames]



with open(args.out_file,"w") as out_file:
    with open(args.out_samples,"w") as out_samples:
        for line in open(args.cytosine_positions, "r"):
            fields = line.strip().split('\t')
            refchro = fields[0]
            refstart = int(fields[1])
            all_vals = {}
            for sample in samples:
                all_vals[sample] = [0,0]
            for n,file in enumerate(files):
                file_loop = True
                file_vals = [0,0]
                while file_loop:
                    if read_go[n]:
                        methline = file.readline()
                    else:
                        methline = read_store[n]
                    if not methline:
                        read_go[n] = True
                        file_loop = False
                        continue
                    meth_fields = methline.strip().split('\t')
                    meth_chro = meth_fields[0]
                    if len(meth_fields) == 7 :
                        if (meth_chro == refchro) & (meth_fields[4].isnumeric()) & (meth_fields[5].isnumeric()) :
                            meth_start = int(meth_fields[1])
                            if meth_start == refstart:
                                file_vals = [int(meth_fields[4]),int(meth_fields[5])]
                                file_loop = False
                                read_store[n] = ""
                                continue
                            elif meth_start < refstart:
                                print(filenames[n])
                                print(line)
                                print(methline)
                                read_go[n] = True
                                read_store[n] = ""
                                continue
                            elif meth_start > refstart:
                                read_go[n] = False
                                read_store[n] = methline
                                file_loop = False
                                continue
                        elif (meth_chro < refchro) or not ((meth_fields[4].isnumeric()) & (meth_fields[5].isnumeric())) :
                            print(filenames[n])
                            print(line)
                            print(methline)
                            read_go[n] = True
                            read_store[n] = ""
                            continue
                        elif (meth_chro > refchro) & (meth_fields[4].isnumeric()) & (meth_fields[5].isnumeric()):
                            read_go[n] = False
                            read_store[n] = methline
                            file_loop = False
                            continue
                    elif len(meth_fields) != 7 :
                        print(filenames[n])
                        print(line)
                        print(methline)
                        read_go[n] = True
                        read_store[n] = ""
                        continue
                all_vals[samples[n]][0] += file_vals[0]
                all_vals[samples[n]][1] += file_vals[1]
            #out_samples.write("\t".join(fields[0:3]+[fields[4]]+[fields[3]]))
            #for sample in samples:
                #vals = all_vals[sample]
                #out_samples.write("\t"+str(vals[0])+"\t"+str(vals[1])+"\n")
                #out_bioreps.write("\n")
            meth = 0
            cov = 0
            for sample in samples:
                vals = all_vals[sample]
                meth += vals[0]
                cov += vals[1]
            try:
                out_file.write("\t".join(fields[0:3]+[fields[4]]+[str(meth)]+[str(cov)]+[fields[3]]+[str(float(meth)/float(cov))])+"\n")
                #out_file.write("\t"+str(float(meth)/float(cov))+"\n")
            except ZeroDivisionError:
                pass
                                 
