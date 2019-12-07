#########################################################################
# File Name: MutPrimerDesign.py
# > Author: CaoYinghao
# > Mail: yhcao@ibmc.pumc.edu.cn
# Created Time: Sat 16 Nov 2019 11:10:43 AM CST
#########################################################################
#! /usr/bin/python

### based: 1
import glob
import sys
import re
import gzip
import collections

class GFF(object):
    def __init__(self):
        pass

class GModel(object):
    GTYPE = ['exon','five_prime_UTR','three_prime_UTR','CDS','mRNA','gene']
    def __init__(self,gchr,start,end,strand,gtype,parent,gid,gname):
        self.chr = gchr
        self.start = start
        self.end = end
        self.strand = strand
        self.type = gtype
        self.parent = parent
        self.id = gid
        self.name = gname

    def __repr__(self):
        return " ".join(["chr:",self.chr,"start:",self.start,"end:",self.end,
            "strand:",self.strand,"type:",self.type,"parent:",self.parent,"id:",self.id,"name:",self.name])

class GMRNA(GModel):
    def __init__(self,*l):
        super(GMRNA,self).__init__(*l)
        self.gmodels = collections.defaultdict(list)
        self.starts = []
        self.abstarts = []
        self.lens = []

    def addModel(self,gmodel):
        self.gmodels[gmodel.type].append(gmodel)

    def constructCDSRegion(self):
        self.starts = []
        self.abstarts = []
        self.lens = []
        if "CDS" in self.gmodels:
            cdss = self.gmodels["CDS"]
            start = 1
            if self.strand == "+":
                for cds in sorted(cdss,key = lambda x:int(x.start)):
                    self.starts.append(int(cds.start))
                    self.abstarts.append(start)
                    lens = int(cds.end) - int(cds.start) + 1
                    self.lens.append(lens)
                    start += lens
                    #print(repr(cds))
            elif self.strand == "-":
                for cds in sorted(cdss,key = lambda x:-int(x.start)):
                    self.starts.append(int(cds.end))
                    self.abstarts.append(start)
                    lens = int(cds.end) - int(cds.start) + 1
                    self.lens.append(lens)
                    start += lens
                    #print(repr(cds))
            #for s,a,l in zip(self.starts,self.abstarts,self.lens):
            #    print("start:",s,"abstart:",a,"len:",l)

        #print(self.__repr__())
        pass


class Gene(GModel):
    def __init__(self,*l):
        super(Gene,self).__init__(*l)
        self.MRNAs = []
        self.annotation = ""

    def addTranscript(self,gmrna):
        self.MRNAs.append(gmrna)

    def construct(self):
        for mrna in self.MRNAs:
            mrna.constructCDSRegion()

    def findLoc(self,loc):
        miss = 0
        results = []
        for mrna in self.MRNAs:
            abstart = mrna.abstarts[0]
            abend = mrna.abstarts[-1] + mrna.lens[-1]
            if loc < abstart or loc > abend:
                miss += 1
                continue
            else:
                for i,(s,l) in enumerate(zip(mrna.abstarts,mrna.lens)):
                    e  = s + l - 1
                    if s <= loc <= e:
                        if self.strand == "-":
                            results.append([mrna.chr,mrna.name,mrna.starts[i] - (loc - s),mrna.strand])
                        elif self.strand == "+":
                            results.append([mrna.chr,mrna.name,mrna.starts[i] + loc - s,mrna.strand])

        if miss == len(self.MRNAs):
            print("Loc:",loc,"is not in transcripts!")
        return results

        pass

    def findABLoc(self,loc):
        pass

    def printSummary(self):
        print(self.__repr__())

class GffLoader(object):
    def __init__(self,gfile):
        self.gfile = gfile
        self.genes = {}
        self.genenames = {}

    def readAtt(self,line):
        parts = line.split(";")
        values = {}
        for p in parts:
            key,value = p.split("=")
            if key == 'ID':
                values['id'] = value
            elif key == "Name":
                values['name'] = value
            elif key == "description":
                values['anno'] = value
            elif key == "Parent":
                values['parent'] = value
        return values

    def load(self):
        genes = {}
        mrnas = {}
        #models = collections.defaultdict(dict)
        models = []
        with open(self.gfile,"r") as infile:
            sys.stdout.write("read gff file" + "\n")
            sys.stdout.flush()
            for line in infile:
                if line.strip().startswith("#"):continue
                parts = line.strip().split("\t")
                gchr = parts[0]
                start = parts[3]
                end = parts[4]
                strand = parts[6]
                values = self.readAtt(parts[-1])
                gid = values.get('id',"")
                gname = values.get('name','')
                gparent = values.get('parent','')
                if parts[2] == "gene":
                    #gene = Gene(gchr=gchr,start=start,end=end,strand=strand,gtype='gene',parent=gparent,gid=gid,gname=gname)
                    gene = Gene(gchr,start,end,strand,parts[2],gparent,gid,gname)
                    genes[gene.id] = gene
                    #print(gene.id)
                    #yield gene
                elif parts[2] == "mRNA" or parts[2] == "transcript":
                    mrna = GMRNA(gchr,start,end,strand,parts[2],gparent,gid,gname)
                    mrnas[mrna.id] = mrna
                elif parts[2] in GModel.GTYPE:
                    model = GModel(gchr,start,end,strand,parts[2],gparent,gid,gname)
                    if model.parent in mrnas:
                        mrnas[model.parent].addModel(model)
                    else:
                        models.append(model)
            sys.stdout.write("merge left models in gff file" + "\n")
            sys.stdout.flush()
            abmodelnum = 0
            for model in models:
                if model.parent in mrnas:
                    mrnas[model.parent].addModel(model)
                else:
                    abmodelnum += 1
            print("Warning::Abnormal models num:",abmodelnum,"ignored")

            sys.stdout.write("merge mrnas in gff file" + "\n")
            sys.stdout.flush()
            for mrnaname,mrna in mrnas.items():
                if mrna.parent in genes:
                    gene = genes[mrna.parent]
                    gene.addTranscript(mrna)
                    self.genenames[gene.name] = gene.id
            self.genes = genes
    def getGenerators(self):
        for name,gene in self.genes.items():
            yield gene
    def getGeneIDLists(self):
        return self.genes

    def getGeneByID(self,gid):
        if gid in self.genes:
            return self.genes[gid]
        else:
            return None

    def getGeneNames(self):
        return self.genenames

    def getGeneByName(self,name):
        if name in self.genenames:
            return self.genes[self.genenames[name]]
        else:
            return None

if __name__ == "__main__":
    #loader = GffLoader("/data1/cancer/genename/ensembl/hg38/Homo_sapiens.GRCh38.87.gff3")
    #loader = GffLoader("/data1/cancer/genename/ensembl/grch37/Homo_sapiens.GRCh37.87.gff3")
    #loader = GffLoader("t.gff3")
    #"/data1/cancer/GATK/hg19/hg19.fasta"
    loader.load()
    ##for g in loader.getGenerators():
    ##    print(g.id)
    gene = loader.getGeneByName("TPP1")
    if gene:
        gene.construct()
        gene.printSummary()
        loc = 1417
        newlocs = gene.findLoc(loc)
        for l in newlocs:
            print(loc,l)

    


