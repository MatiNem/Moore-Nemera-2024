#!/usr/bin/env python

#Script to create an aggregate window annotation file from a given annotation (extends annotation by f, moves into the region of interest by r, splits everything up by n and m)
from optparse import OptionParser, OptionGroup
import glob
import os
import subprocess
import re
import tempfile


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-a", "--anno", type="string", dest="annoname",help="Name for the annotation file to intersect with methylation", metavar="INFILE")
    parser.add_option("-n", "--nbins", type="int", dest="n",help="Number of bins in each flanking region", metavar="INFILE")
    parser.add_option("-f", "--flank", type="int", dest="flank",help="Distance to examine outside the region", metavar="INFILE")
    parser.add_option("-o", "--outfile", type="string", dest="out",help="Distance to move into the annotation region")
    (options, args) = parser.parse_args()

    cwd = os.getcwd()
    if options.annoname[0] == "/":
        anno_file = options.annoname
    else: anno_file = cwd+"/"+options.annoname
    base_name = anno_file.split(".")[0]
    windowfile = open(options.out,"w")
    nrange = int(options.flank) / int(options.n)
    for line in open(options.annoname,"r"):
        fields = line.split("\t")
        chro =  fields[0]
        start = (float(fields[1])+float(fields[2]))/2
        end = (float(fields[1])+float(fields[2]))/2
        if fields[-1][-1] == "\n":
            fields[-1] = fields[-1][0:-1]
        for field in fields:
            if field == "-" or field == "+":
                strand = field
        anno = fields[3]
        for n in range(0,int(options.n)*2):
            m = n*nrange-int(options.flank)
            num = n-int(options.n)
            wstart = start + m
            wend = wstart + nrange-1
            strand = "+"
            if strand == "-":
                num = -num-1
                wstart = end + m
                wend = wstart + nrange-1
            if wstart < 0:
                wstart = 0
                wend = max(wend,1)
            windowfile.write("\t".join([chro,str(int(wstart)),str(int(wend)),anno,"f_"+str(int(num)),strand]))
            windowfile.write("\n")
    windowfile.close()