#########################################################################
# File Name: MutPrimerDesign.py
# > Author: CaoYinghao
# > Mail: yhcao@ibmc.pumc.edu.cn
# Created Time: Sat 16 Nov 2019 11:10:43 AM CST
#########################################################################
#! /usr/bin/python

import glob
import sys
import re
from Bio import SeqIO

class FastaLoc(object):
    ###start : based 1
    ###end : based 1, included
    def __init__(self,fastaFile,segments):
        self.fastafile = fastaFile
        self.segments = segments
        self.results = {}
        #self.chrom = chrom
        #self.start = start
        #self.end = end

    def getSeq(self):
        results = {}
        for gene in SeqIO.parse(self.fastafile,'fasta'):
            for segment in  self.segments:
                chrom,mrnaid,start,end,strand,tmploc = segment
                if gene.id == chrom:
                    seq1 = gene.seq[int(start)-1:int(end)]
                    if strand == "-":
                        results[mrnaid + "|" + chrom + "|" + start + "_" + end+ "|center:" + str(tmploc)] = seq1.reverse_complement()
                    else:
                        results[mrnaid + "|" + chrom + "|" + start + "_" + end+"|center:"+str(tmploc)] = seq1
                    if len(results) == len(self.segments):break
        self.results = results
        return results
        #print ("Reverse:" + seq1.reverse_complement())


if __name__ == "__main__":
    if(len(sys.argv) < 2):
        print ("fastaLoc.py fastaFile chr start end")
        exit()
    name = sys.argv[1]
    chrom = sys.argv[2]
    loc1 = loc2 = ""
    #outfile = ""
    if len(sys.argv) > 3:
        loc1 = sys.argv[3]
        #outfile = ("".join(name.split('.')[:-1]) + "_" + loc1 + ".fasta","w")
    if len(sys.argv) > 4:
        loc2 = sys.argv[4]
        #outfile = ("".join(name.split('.')[:-1]) + "_" + loc1 + "_" + loc2 + ".fasta", "w")
    for gene in SeqIO.parse(name,'fasta'):
        if gene.id == chrom:
            if loc1 != "" and loc2 != "":
                print (">"+gene.id + "_" + loc1 + "_" + loc2)
                seq1 = gene.seq[int(loc1)-1:int(loc2)]
                print (seq1)
            else:
                print (">"+gene.id)
                print(str(gene.seq))
            #print ("Reverse:" + seq1.reverse_complement())
            exit()

