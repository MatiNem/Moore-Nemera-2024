#!/usr/bin/env python
#Change calculation: average methylation level of windows such that each window is weighted equally
#User-friendly wrapper script that takes in a bed or CGmap file, runs bedtools intersect
from optparse import OptionParser, OptionGroup
import glob
import os
import subprocess
import re
import tempfile


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="Name of the Methylation data of interest", metavar="INFILE")
    parser.add_option("-a", "--anno", type="string", dest="annoname",help="Name for the annotation file to intersect with methylation", metavar="INFILE")
    parser.add_option("-o", "--outdir", type="string", dest="outdir",help="Name of the directory for the outputs", metavar="OUTFILE")
    parser.add_option("-j", "--job", type="string", dest="jobname",help="Name for the resulting job file that is submitted", metavar="OUTFILE")
    parser.add_option("-n", "--nbins", type="int", dest="n",help="Number of bins in flanking regions", metavar="INFILE")
    parser.add_option("-m", "--mbins", type="int", dest="m",default=0,help="Number of bins in between annotation regions", metavar="INFILE")
    parser.add_option("-f", "--flank", type="int", dest="flank",help="Distance to examine outside the region", metavar="INFILE")
    parser.add_option("-r", "--region", type="int", dest="region",default=0,help="Distance to move into the annotation region")
    parser.add_option("-e", "--regionside", type="string", dest="rside",default="both",help="Side of the annotation to advance into (left indicates upstream side, right indicates downstream side, both indicates both sides")
    parser.add_option("-s", "--sided", action="store_true", dest="sided",default=False,help="Option to make the pipeline one-sided, centered on the first entry if + stranded, second if -")
    parser.add_option("-g", "--meta", type="int", dest="meta",default=0,help="Option to make the pipeline plot meta-genes, with a number indicating how many bins the meta-gene will have. Mutually exclusive with -s. Compatible with regions, meta-gene will plot the the area of the gene between the -r regions")
    parser.add_option("-t", "--type", type="string", dest="type",default="m",help="Determines type of aggregation and calculation. 'c' for Chip will make aggregate_windows.py use log2(c5/c6), 'm' for Methylation will use c5/c6, 'h' for hmC will use straight means,'t' will use c5/c6 * coverage")
    parser.add_option("-w", "--weighted", action="store_true", dest="weight",default=False,help="Option to make the pipeline scale the signal file by the number of bases the signal covers")
    parser.add_option("-c", "--centered", action="store_true", dest="center",default=False,help="Option to make the pipeline use the centers of regions. Incompatible with -s, -m, and -r")
    parser.add_option("-p", "--chip_input", type="string", dest="chip_input",help="Name of the Input data of interest", default = "",metavar="INFILE")
    (options, args) = parser.parse_args()

    #if options.sided and options.meta:
    #    raise ValueError('cannot have both metagenes and single-sided plotting')

    cwd = os.getcwd()+"/"
    if options.infilename[0] == "/":
        data_file = options.infilename
    else: data_file = cwd+options.infilename

    if options.annoname[0] == "/":
        anno_file = options.annoname
    else: anno_file = cwd+options.annoname
    #base_name = data_file.split("/")
    #base_name = base_name[len(base_name)-1].split(".")[0]
    #base_name = anno_file.split(".")[0]
    #base_name = base_name.split("/")[-1]
    #base_name = options.outdir+"/"+base_name
    base_name = options.jobname

    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)
    os.chdir(options.outdir)

    windowname = base_name+"_windows.bed"

    script_dir = os.path.dirname(os.path.realpath(__file__))
    #subprocess.call(["echo",options.jobname])
    #subprocess.call(["echo",os.getcwd()])

    #Write pipeline for each run, Convert map to bed, if necessary, make window file, bedtools map, convert to .out file.
    if not options.center:
        subprocess.call(["python",script_dir+"/make_agg_windows.py","-a",anno_file,"-n",str(options.n),"-m",str(options.m),"-f",str(options.flank),"-r",str(options.region),"-o",windowname,"-g",str(options.meta),"-e",str(options.rside)])
    elif options.center:
        subprocess.call(["python",script_dir+"/make_agg_windows_center.py","-a",anno_file,"-n",str(options.n),"-f",str(options.flank),"-o",windowname])

    sort_bed = open(windowname+"_sort", "w")
    subprocess.call(["sort","-k","1,1","-k2,2n",windowname],stdout=sort_bed)
    sort_bed.close()
    mapped_bed = open(options.jobname+".bed", "w")
    if options.type == "m":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","5,6","-o","sum"],stdout=mapped_bed)
    elif options.type == "t":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","5,6","-o","sum"],stdout=mapped_bed)
        cov_bed = open(options.jobname+"_cov.bed", "w")
        subprocess.call(["bedtools","coverage","-a",options.jobname+".bed","-b",data_file,"-sorted","-counts"],stdout=cov_bed)
    elif options.type == "h":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","5,6","-o","mean"],stdout=mapped_bed)
    elif options.type == "th":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","5,6","-o","sum"],stdout=mapped_bed)
    elif options.type == "d":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","10","-o","sum"],stdout=mapped_bed)
    elif options.type == "c":
        chip_hist = subprocess.Popen(("bedtools","coverage","-a",windowname+"_sort","-b",data_file,"-sorted","-hist"), stdout=subprocess.PIPE)
        if options.chip_input:
            chip_cov = subprocess.Popen(('python', script_dir+'/group_hist_cov_2.py',"-i","stdin","-o","stdin_cov","-n",data_file), stdin=chip_hist.stdout,stdout=subprocess.PIPE)
            input_hist = subprocess.Popen(("bedtools","coverage","-a","stdin","-b",options.chip_input,"-sorted","-hist"),stdin=chip_cov.stdout,stdout=subprocess.PIPE)
            subprocess.call(('python', script_dir+'/group_hist_cov_2.py',"-i","stdin","-o","stdin_cov","-n",options.chip_input), stdin=input_hist.stdout,stdout=mapped_bed)
            chip_cov.wait()
            input_hist.wait()
        else:
            chip_cov = subprocess.call(('python', script_dir+'/group_hist_cov_2.py',"-i","stdin","-o","stdin_cov","-n",data_file), stdin=chip_hist.stdout,stdout=mapped_bed)
        chip_hist.wait()
    elif options.type == "cbg":
        chip_hist = subprocess.Popen(("bedtools","intersect","-wao","-sorted","-a",windowname+"_sort","-b",data_file), stdout=subprocess.PIPE)
        if options.chip_input:
            chip_cov = subprocess.Popen(('python', script_dir+'/group_hist_bg.py',"-i","stdin","-o","stdin_cov"), stdin=chip_hist.stdout,stdout=subprocess.PIPE)
            input_hist = subprocess.Popen(("bedtools","coverage","-a","stdin","-b",options.chip_input,"-sorted","-hist"),stdin=chip_cov.stdout,stdout=subprocess.PIPE)
            subprocess.call(('python', script_dir+'/group_hist_bg.py',"-i","stdin","-o","stdin_cov"), stdin=input_hist.stdout,stdout=mapped_bed)
            chip_cov.wait()
            input_hist.wait()
        else:
            chip_cov = subprocess.call(('python', script_dir+'/group_hist_bg.py',"-i","stdin","-o","stdin_cov"), stdin=chip_hist.stdout,stdout=mapped_bed)
        chip_hist.wait()
    elif options.type == "bg":
        bg_inter = subprocess.Popen(("bedtools","intersect","-a",windowname+"_sort","-b",data_file,"-sorted","-wao"), stdout=subprocess.PIPE)
        chip_cov = subprocess.call(('python', script_dir+'/group_hist_bg.py',"-i","stdin","-o","stdin_cov"), stdin=bg_inter.stdout,stdout=mapped_bed)
    elif options.type == "g":
        subprocess.call(["bedtools","coverage","-a",windowname+"_sort","-b",data_file,"-counts","-sorted"],stdout=mapped_bed)
    elif options.type == "mbg":
        subprocess.call(["bedtools","map","-a",windowname+"_sort","-b",data_file,"-c","4","-o","mean"],stdout=mapped_bed)
    mapped_bed.close()
    if options.type == "m":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t",options.type])
    elif options.type == "t":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+"_cov.bed","-t","m","-o",options.jobname+".out"])
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+"_cov.bed","-t","t"])
    elif options.type == "h":
        subprocess.call(["python",script_dir+"/aggregate_windows_hmC.py","-i",options.jobname+".bed","-t","h"])
        subprocess.call(["python",script_dir+"/aggregate_windows_hmC.py","-i",options.jobname+".bed","-t","m"])
    elif options.type == "th":
        subprocess.call(["python",script_dir+"/aggregate_windows_hmC.py","-i",options.jobname+".bed","-t","th"])
        subprocess.call(["python",script_dir+"/aggregate_windows_hmC.py","-i",options.jobname+".bed","-t","tm"])
    elif options.type == "c":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t","c"])
    elif options.type == "d":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t",options.type])
    elif options.type == "mbg":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t","g"])
    elif options.type == "cbg":
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t","g"])
    else:
        subprocess.call(["python",script_dir+"/aggregate_windows_3.py","-i",options.jobname+".bed","-t",options.type])