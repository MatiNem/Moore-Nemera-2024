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
    parser.add_option("-m", "--mbins", type="int", dest="m",help="1/2 Number of bins in between annotation regions", metavar="INFILE")
    parser.add_option("-f", "--flank", type="int", dest="flank",help="Distance to examine outside the region", metavar="INFILE")
    parser.add_option("-r", "--region", type="int", dest="region",default=0,help="Distance to move into the annotation region")
    parser.add_option("-o", "--outfile", type="string", dest="out",help="Distance to move into the annotation region")
    parser.add_option("-e", "--regionside", type="string", dest="rside",default="both",help="Side of the annotation to advance into (1, up, l, left all indicate upstream side, 2, down, r, right all indicate downstream side, 3, both, indicate both sides")
    parser.add_option("-g", "--meta", type="int", dest="meta",default=0,help="Option to make the pipeline plot meta-genes, with a number indicating how many bins the meta-gene will have. Mutually exclusive with -s")
    (options, args) = parser.parse_args()

    cwd = os.getcwd()
    if options.annoname[0] == "/":
        anno_file = options.annoname
    else: anno_file = cwd+"/"+options.annoname
    base_name = anno_file.split(".")[0]
    windowfile = open(options.out,"w")
    try:
        nrange = int(options.flank) / int(options.n)
    except ZeroDivisionError:
        nrange=0
    for line in open(options.annoname,"r"):
        fields = line.split("\t")
        chro =  fields[0]
        start = int(fields[1])
        end = int(fields[2])
        if end - start > options.region:
            anno = False
            if fields[-1][-1] == "\n":
                fields[-1] = fields[-1][0:-1]
            for field in fields[3:]:
                if field == "-" or field == "+":
                    strand = field
                if field[0].isalpha():
                    anno = field
            if not anno:
                anno = fields[3]
            anno = fields[3]
            for n in range(0,int(options.n)):
                m = n*nrange
                num1 = -(n+1)
                num2 = n+1
                side1 = False
                side2 = False
                if options.rside == "left":
                    side1 = True
                    side2 = False
                if options.rside == "right":
                    side1 = False
                    side2 = True
                if strand == "-":
                    num1 = n+1
                    num2 = -(n+1)
                    side1 = not side1
                    side2 = not side2
                side1 = True
                side2 = True
                wend1 = start - m
                wstart1 = wend1 - nrange+1
                wstart2 = end + m
                wend2 = wstart2 + nrange-1
                if wstart1 > 0 and side1:
                    windowfile.write("\t".join([chro,str(int(wstart1)),str(int(wend1)),anno,"f_"+str(int(num1)),strand]))
                    windowfile.write("\n")
                if wstart2 > 0 and side2:
                    windowfile.write("\t".join([chro,str(int(wstart2)),str(int(wend2)),anno,"f_"+str(int(num2)),strand]))
                    windowfile.write("\n")
            if options.m > 0 and options.region > 0:
                mrange = int(options.region)/ int(options.m)
                for n in range(0,int(options.m)):
                    m = n*mrange
                    wstart1 = start + m+1
                    wend1 = wstart1 + mrange-1
                    wend2 = end - m-1
                    wstart2 = wend2 - mrange+1
                    num1 = n+1
                    num2 = 2*int(options.m)-n
                    side1 = False
                    side2 = False
                    if options.rside == "left":
                        side1 = True
                        side2 = False
                    if options.rside == "right":
                        side1 = False
                        side2 = True
                    if strand == "-":
                        num1 = 2*int(options.m)-n
                        num2 = n+1
                        side1 = not side1
                        side2 = not side2
                    if options.rside == "both":
                        side1 = True
                        side2 = True
                    if wstart1 > 0 and side1 and wend1 < end:
                        windowfile.write("\t".join([chro,str(int(wstart1)),str(int(wend1)),anno,"r_"+str(int(num1)),strand]))
                        windowfile.write("\n")
                    if wstart2 > 0 and side2 and wstart2 > start:
                        windowfile.write("\t".join([chro,str(int(wstart2)),str(int(wend2)),anno,"r_"+str(int(num2)),strand]))
                        windowfile.write("\n")
            if options.meta:
                if options.rside == "both":
                    grange = float((end - start - 2*options.region))/float(int(options.meta))
                else:
                    grange = float((end - start - options.region))/float(int(options.meta))
                if grange >= 1:
                    for n in range(0,options.meta):
                        m = n*grange
                        num = n+1
                        if options.rside == "right":
                            mstart = start
                            mend = end
                        else:
                            mstart = start + options.region + 1
                            mend = end - options.region - 1
                        if strand == "+":
                            wstart = mstart + m
                            wend = wstart + grange-1
                        elif strand == "-":
                            wend = mend - m
                            wstart = wend - grange+1

                        windowfile.write("\t".join([chro,str(int(wstart)),str(int(wend)),anno,"m_"+str(int(num)),strand]))
                        windowfile.write("\n")
    windowfile.close()