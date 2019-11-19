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
import os
import argparse
from Bio import SeqIO
from gffLoader import GffLoader
from fastaLoc import FastaLoc
import primer3
import collections

import gzip
import pandas as pd
from pickle import dump,load

class VariationLoc(object):
    def __init__(self,fastaFile,genenames,gffFile,outfile,length=200):
        self.fastaFile = fastaFile
        self.genenames = genenames
        self.gffFile = gffFile
        self.results = {}
        self.outfile = outfile
        self.length = length
        #loc = 1417

    def getFlankSeq(self,hgncs):
        loader = GffLoader(self.gffFile)
        loader.load()

        lists = []
        for genename in self.genenames:
            print(genename)
            genename,ablocstr,ilen = genename.split(":")
            if genename not in hgncs:
                print("Warning: Gene (",genename,") not in hgnc ids,ignored")
                exit()
                continue
            frontloc = 0
            afterloc = 0
            if ablocstr.find("-") > -1:
                abloc,frontloc = ablocstr.split("-")
            elif ablocstr.find("+") > -1:
                abloc,afterloc = ablocstr.split("+")
            else:
                abloc = ablocstr

            gene = loader.getGeneByName(hgncs[genename].split("|")[0])
            if gene:
                gene.construct()
                #gene.printSummary()
                newlocs = gene.findLoc(int(abloc))
                n = 0
                for l in newlocs:
                    chrom,mrnaid,loc,strand = l
                    tname = mrnaid.split(".")[0]
                    if tname != hgncs[genename].split("|")[-1]:
                        n += 1
                        continue
                    if gene.strand == "-":
                        loc = loc + int(frontloc) - int(afterloc)
                    else:
                        loc = loc - int(frontloc) + int(afterloc)
                    #if chrom.find("chr") == -1 and not chrom.startswith("N"):
                    #    chrom  = "chr" + chrom 
                    lists.append([chrom,genename + "|"+ mrnaid,str(int(loc)-self.length),str(int(loc)+self.length),strand,str(loc) + "|" + ablocstr])
                if n == len(newlocs):
                    print("Waring: Transcript (" + hgncs[genename] + ") is not found in Gene (" + genename + ") in gff file, ignored")
        #print(lists)
        fl = FastaLoc(self.fastaFile,lists)
        results = fl.getSeq()
        newresults = collections.defaultdict(dict)
        for r in sorted(results.keys()):
            tag = 0
            seq = results[r]
            gname = r.split("|")[0]
            tname = r.split("|")[1].split(".")[0]
            gloc = r.split("|")[-1]
            for inputgname in self.genenames:
                genename,loc,length = inputgname.split(":")
                #if gname == genename and gloc == loc and gbase == base:
                if gname == genename and gloc == loc:
                    gbase = str(seq[self.length:self.length + int(length)])
                    newresults[r+ "|Refbase:"+gbase+"|"+length]['name'] = inputgname
                    newresults[r+ "|Refbase:"+gbase+"|"+length]['out'] = results[r]
                    tag = 1
            if tag == 0:
                print("Warning:",gname,"with",gloc,"can't be found,ignored")
        self.results = newresults
        if len(newresults) == 0:
            print(len(newresults),"amplicons found; You may check the source of gff and genome file you downloaded. NCBI's annotation is needed. Like this:")
            print(" (genome: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz)")
            print(" (gff: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz)")
        else:
            print(len(newresults),"amplicons found")
        self.writeSeq(self.outfile)
        #return newresults

    #def readHGNC(self):
    #    hgncs = {}
    #    #with open("/data1/cancer/genename/hgnc_complete_set.txt","r") as infile:
    #    with open("/home/yhcao/Drive/other/GeneTest/script/hgnc_complete_set.txt","r") as infile:
    #        for line in infile:
    #            parts =line.split("\t")
    #            if parts[3] != "protein-coding gene":continue
    #            hgncs[parts[1]] = parts[23]
    #            #print(parts[1],parts[23])
    #    return hgncs

    def writeSeq(self,outfile):
        out = open(outfile,"w")
        #print(len(self.results),"results found")
        for r in sorted(self.results.keys()):
            name = self.results[r]['name']
            result = self.results[r]['out']
            out.write(">" + name + "|"+ r +"\n")
            out.write(str(result)+"\n")
        out.close()

    def printSeq(self):
        for r in self.results.keys():
            print(">"+r+"|"+"\n"+str(results[r]))
            print("." * self.length + "*"+"." * self.length)

class PrimerDesign(object):
    def __init__(self):
        pass

    @staticmethod
    def readFasta(fastaFile):
        genes = {}
        for gene in SeqIO.parse(fastaFile,"fasta"):
            genes[gene.id] = gene
        return genes

    @staticmethod
    def designprimer(gid,input_seq,num,args):
        results = {}
        primer = (primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': gid,
                'SEQUENCE_TEMPLATE': input_seq,
                #'SEQUENCE_EXCLUDED_REGION': [0, 0],
                'SEQUENCE_INCLUDED_REGION': [int(int(args.maxLength)/2)-200, int(int(args.maxLength)/2) + 200],
                'SEQUENCE_TARGET': [int(int(args.maxLength)/2)-80, 150],
                #[120,110],
                #SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
                #SEQUENCE_EXCLUDED_REGION
            },
            {
                'PRIMER_TASK': 'generic',
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_NUM_RETURN': int(args.primerNum),
                'PRIMER_OPT_SIZE': int(args.primer_opt_size),
                'PRIMER_MIN_SIZE': int(args.primer_min_size),
                'PRIMER_MAX_SIZE': int(args.primer_max_size),
                'PRIMER_OPT_TM': int(args.primer_opt_tm),
                'PRIMER_MIN_TM': int(args.primer_min_tm),
                'PRIMER_MAX_TM': int(args.primer_max_tm),
                'PRIMER_MIN_GC': int(args.primer_min_gc),
                'PRIMER_MAX_GC': int(args.primer_max_gc),
                'PRIMER_MAX_POLY_X': int(args.primer_max_poly_num),
                'PRIMER_INTERNAL_OPT_TM':int(args.primer_opt_tm) + 5,
                'PRIMER_INTERNAL_MAX_TM':int(args.primer_max_tm) + 5,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 12,
                'PRIMER_MAX_SELF_END': 8,
                'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                'PRIMER_PAIR_MAX_COMPL_END': 8,
                'PRIMER_PRODUCT_SIZE_RANGE': [[len(input_seq)-int(int(args.maxLength) / 2),len(input_seq)]],
                }
            ))

        
        start = ['PRIMER_LEFT','PRIMER_RIGHT','PRIMER_INTERNAL']
        end = ['SEQUENCE','TM','GC_PERCENT','SELF_ANY_TH','SELF_END_TH','HAIRPIN_TH','END_STABILITY']
        pairs = ['PRIMER_PAIR']
        paire = ['COMPL_ANY_TH','COMPL_END_TH','PRODUCT_SIZE']
        gname = gid.split("|")[0]
        out = ""
        out1 = ""
        print('PRIMER_LEFT_EXPLAIN:',primer['PRIMER_LEFT_EXPLAIN'])
        print('PRIMER_RIGHT_EXPLAIN:',primer['PRIMER_RIGHT_EXPLAIN'])
        print('PRIMER_INTERNAL_EXPLAIN:',primer['PRIMER_INTERNAL_EXPLAIN'])
        print('PRIMER_PAIR_EXPLAIN:',primer['PRIMER_PAIR_EXPLAIN'])
        for i in range(primer['PRIMER_PAIR_NUM_RETURNED']):
            lpstart,lplen = primer[start[0] + "_" + str(i)]
            rpstart,rplen = primer[start[1] + "_" + str(i)]
            ipstart,iplen = primer[start[2] + "_" + str(i)]
            size = primer[pairs[0] + "_" + str(i) + "_" + paire[-1]]
            comanyth = primer[pairs[0] + "_" + str(i) + "_" + paire[0]]
            comendth = primer[pairs[0] + "_" + str(i) + "_" + paire[1]]
            lpseq = primer[start[0] + "_" + str(i) + "_" + end[0]]
            rpseq = primer[start[1] + "_" + str(i) + "_" + end[0]]
            ipseq = primer[start[2] + "_" + str(i) + "_" + end[0]]
            lptm = round(primer[start[0] + "_" + str(i) + "_" + end[1]],2)
            rptm = round(primer[start[1] + "_" + str(i) + "_" + end[1]],2)
            iptm = round(primer[start[2] + "_" + str(i) + "_" + end[1]],2)
            lpgc = round(primer[start[0] + "_" + str(i) + "_" + end[2]],2)
            rpgc = round(primer[start[1] + "_" + str(i) + "_" + end[2]],2)
            ipgc = round(primer[start[2] + "_" + str(i) + "_" + end[2]],2)
            lselanyth = round(primer[start[0] + "_" + str(i) + "_" + end[3]],2)
            rselanyth = round(primer[start[1] + "_" + str(i) + "_" + end[3]],2)
            iselanyth = round(primer[start[2] + "_" + str(i) + "_" + end[3]],2)
            lselendth = round(primer[start[0] + "_" + str(i) + "_" + end[4]],2)
            rselendth = round(primer[start[1] + "_" + str(i) + "_" + end[4]],2)
            iselendth = round(primer[start[2] + "_" + str(i) + "_" + end[4]],2)
            lselhairth = round(primer[start[0] + "_" + str(i) + "_" + end[5]],2)
            rselhairth = round(primer[start[1] + "_" + str(i) + "_" + end[5]],2)
            iselhairth = round(primer[start[2] + "_" + str(i) + "_" + end[5]],2)
            lselstable = round(primer[start[0] + "_" + str(i) + "_" + end[6]],2)
            rselstable = round(primer[start[1] + "_" + str(i) + "_" + end[6]],2)
            out += "Amplicon" + str(num) + "_" + str(i) + ":\n " + gid + " " + str(i) + "\n"
            out += " ".join([" amplicon length:",str(size)]) + "\n"
            out += " ".join([" primer informations:"]) + "\n"
            out += " ".join(["        primer pair:","AnyTH:"+str(comanyth)+"|EndTH:"+str(comendth)]) + "\n"
            out += " ".join(["     left primer(>):",lpseq +"|Loc:"+str(lpstart)+"|Len:"+str(lplen)+"|TM:"+str(lptm)+"|GC:"+str(lpgc) + "|AnyTH:" + str(lselanyth)+"|EndTH:" + str(lselendth) + "|HairpinTH:" + str(lselhairth) + "|Stability:" + str(lselstable)]) + "\n"
            out += " ".join(["    right primer(<):",rpseq+"|Loc:"+str(rpstart)+"|Len:"+str(rplen)+"|TM:"+str(rptm)+"|GC:"+str(rpgc) + "|AnyTH:" + str(rselanyth)+"|EndTH:" + str(rselendth) + "|HairpinTH:" + str(rselhairth) + "|Stability:" + str(rselstable)]) + "\n"
            out += " ".join(["  internal oligo(|):",ipseq+"|Loc:"+str(ipstart)+"|Len:"+str(iplen)+"|TM:"+str(iptm)+"|GC:"+str(ipgc) + "|AnyTH:" + str(iselanyth)+"|EndTH:" + str(iselendth) + "|HairpinTH:" + str(iselhairth)]) + "\n"
            out += input_seq + "\n"

            iseqs = list("-"*len(input_seq))
            iseqs[:lpstart] = "." * lpstart
            iseqs[rpstart:] = "." * (int(args.maxLength)-rpstart)
            iseqs[lpstart:lpstart+lplen] = ">" * lplen
            iseqs[rpstart-rplen:rpstart] = "<" * rplen
            iseqs[ipstart:ipstart+iplen] = "|" * iplen
            iseqs[int(int(args.maxLength)/2)] = "*" * int(gid.split("|")[-1])
            newseq = "".join(iseqs)

            #lpart = ipstart - lplen - lpstart
            #rpart = size - lplen - ipart - rplen -1
            #print("."*lpstart + ">"*lplen + (size - lplen -rplen)*"-" + "<"*rplen + "."*(401-rpstart))
            out += newseq + "\n" * 2
            #out += " ".join(["."*lpstart + ">"*lplen + lpart*"-" + "*"+rpart*"-" + "<"*rplen + "."*(int(args.maxLength) - rpstart)]) + "\n"*2

            out1 += "{0:25s}{1:<30}{2:>30}{3:>40}".format(gname + "_" + str(i),lpseq,rpseq,ipseq) + "\n"


        return out,out1
    @staticmethod
    def design(genes,args):

        out1 = open(args.outFile+".detail","w")
        out2 = open(args.outFile+".primer","w")
        out2.write("{0:25s}{1:<30}{2:>30}{3:>40}".format("#VariationID","#LeftPrimer","#RightPrimer","#InternalOligo")+"\n")
        n = 0
        for gid in sorted(genes.keys()):
            gene = genes[gid]
            gname = gid.split("|")[0]
            seq = str(gene.seq)
            seqout,primerout = PrimerDesign.designprimer(gid,seq,n,args)
            out1.write(seqout)
            out2.write(primerout)
            if seqout != "":
                n += 1
        print(str(n) + " primerdesign jobs done!")
        out1.close()
        out2.close()


class Process(object):
    def __init__(self):
        self.get_parser()
        pass

    def get_parser(self):
        desc = """Program: MutPrimerDesign
  Version: 0.1
  Author : Yinghao Cao
  Email  : <yhcao@ibms.pumc.edu.cn>
        """

        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc)
        igroup = parser.add_mutually_exclusive_group(required=True)
        igroup.add_argument('-i','--input',help="Input string for mutation primer design, (format:  genename:cds_location:length, separated by comma. For example:  TPP1:1417,CLN6:794:4,WAS:223:5).")
        igroup.add_argument('-f','--file',help="Input file for mutation primer design, each line stands for one mutation.")
        igroup.add_argument('-p','--primerFile',action = "store_true", help="Amplicon sequence file for mutation primer design, ignore sequence extraction.")
        parser.add_argument('-g','--genome', required = True, help="Human genome sequence file (from NCBI)." )
        parser.add_argument('-G','--gffFile', required = True, help="Human genome annotation file (from NCBI)." )
        parser.add_argument('-k','--nameDB', required = True, help="Name convertion database file ." )
        parser.add_argument('-o','--outFile', required = True, help="Output file for amplicon sequences." )
        parser.add_argument("-m",'--maxLength', default = 400, help="Amplicon maximum length. (Default:400bp)")
        parser.add_argument("-n",'--primerNum', default = 1, help="Number of primer pairs for each amplicon. (Default:1)")
        parser.add_argument("-s",'--geneSource', default = 'gene_symbol', help="Type of genenames. (Default:gene_symbol,[ensembl_id],[hgnc_id])")
        parser.add_argument('--primer_opt_size', default = 20, help="Primer optimize size. (Default:20)")
        parser.add_argument('--primer_min_size', default = 18, help="Primer minimum size. (Default:18)")
        parser.add_argument('--primer_max_size', default = 27, help="Primer maximum size. (Default:27)")
        parser.add_argument('--primer_opt_tm', default = 60, help="Primer optimize temperature. (Default:60)")
        parser.add_argument('--primer_min_tm', default = 57, help="Primer minimum temperature. (Default:57)")
        parser.add_argument('--primer_max_tm', default = 63, help="Primer maximum temperature. (Default:63)")
        parser.add_argument('--primer_min_gc', default = 20, help="Primer minimum GC. (Default:20)")
        parser.add_argument('--primer_max_gc', default = 80, help="Primer maximum GC. (Default:80)")
        parser.add_argument('--primer_max_poly_num', default = 5, help="Primer maximum ploy base number. (Default:5)")
        self.parser = parser


    def dump(self):
        f = pd.read_csv("hgnc_complete_set.txt",sep="\t")
        cols = ['hgnc_id','symbol','ensembl_gene_id','refseq_accession']
        nf = f[cols]
        nf['hgnc_new_id'] = nf['hgnc_id'].apply(lambda x:x.split(":")[-1])
        nf['new_refseq_accession'] = nf['symbol'] + "|" + nf['refseq_accession']
        symbol = dict(zip(nf['symbol'],nf['new_refseq_accession']))
        hgnc = dict(zip(nf['hgnc_new_id'],nf['new_refseq_accession']))
        ensembl = dict(zip(nf['ensembl_gene_id'],nf['new_refseq_accession']))
        alldic = {"gene_symbol":symbol,"hgnc_id":hgnc,"ensembl_id":ensembl}
        out = gzip.open("namedb","wb")
        dump(alldic,out)
        pass

    def load(self,mytype):
        handler = gzip.open("namedb","rb")
        names = load(handler)
        #print(names.keys())
        return names[str(mytype)]



    def design(self,args,genes):
        genes = PrimerDesign.readFasta(args.outFile)
        PrimerDesign.design(genes,args)


    def run_cmd(self):
        args = self.parser.parse_args()
        if args.geneSource not in ['gene_symbol','ensembl_id','hgnc_id']:
            self.parser.print_help()
            exit()


        #print(args)
        fastaFile = args.genome
        gffFile = args.gffFile

        genenames = []
        if args.input:
            genenames = args.input.split(";")
        elif args.file:
            with open(args.file,"r")as infile:
                for i in infile:
                    if i.strip() != "":
                        genenames.append(i.strip())
        if args.primerFile:
            if not os.path.exists(args.outFile):
                print("Amplicon sequence " + args.outFile + " not found, please extract sequences using -i or -f commands first." )
            else:
                self.design(args,genenames)
        else:
            print("GeneNames:",genenames)
            names = self.load(args.geneSource)
            varloc = VariationLoc(fastaFile,genenames,gffFile,args.outFile,int(int(args.maxLength)/2))
            varloc.getFlankSeq(names)
            self.design(args,genenames)



        pass


if __name__ == "__main__":
    p = Process()
    p.run_cmd()
    exit()
