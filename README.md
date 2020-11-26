# MutPrimerDesign: Primer design for human mutations

DNA mutations located in the coding region of a gene are closely related to the function of the gene. When the mutation site of the human gene protein coding region is known, how to design primers on the genome to validate the mutation becomes an important issue. In this study, MutPrimerDesign was developed as a primer design program by using Python language. Through analyzing human genome sequence and gene annotation information, MutPrimerDesign converts the gene coding region coordinates into genomic coordinates and calls the Python package interface of Primer3, which can automatically complete the primer and probe sequence design for the gene mutation in batches. MutPrimerDesign is easy to use, can recognize gene names in various databases, and can modify the general parameters of primers, thus achieving rapid adjustment of primers.
MutPrimerDesign is maintained by Yinghao Cao [yhcao@ibms.pumc.edu.cn].
## Download and Installation
```
git clone https://github.com/bioinfo-ibms-pumc/MutPrimerDesign.git

Human genome file:
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip -c GCF_000001405.39_GRCh38.p13_genomic.fna.gz > genome.fa
Human gene annotation file:
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
gunzip -c GCF_000001405.39_GRCh38.p13_genomic.gff.gz > gene.gff

```

## Dependency package installation:
```
pip3 install biopython pandas primer3-py

```
## Command Lines

```  
MutPrimerDesign.py [-h] (-i INPUT | -f FILE | -p) -g GENOME -G GFFFILE
                          -k NAMEDB -o OUTFILE [-m MAXLENGTH] [-n PRIMERNUM]
                          [-s GENESOURCE] [--primer_opt_size PRIMER_OPT_SIZE]
                          [--primer_min_size PRIMER_MIN_SIZE]
                          [--primer_max_size PRIMER_MAX_SIZE]
                          [--primer_opt_tm PRIMER_OPT_TM]
                          [--primer_min_tm PRIMER_MIN_TM]
                          [--primer_max_tm PRIMER_MAX_TM]
                          [--primer_min_gc PRIMER_MIN_GC]
                          [--primer_max_gc PRIMER_MAX_GC]
                          [--primer_max_poly_num PRIMER_MAX_POLY_NUM]

        

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input string for mutation primer design, (format:
                        genename:cds_location:length, separated by comma. For
                        example: KRAS:34:2,MET:1124:1,EGFR:2573:1).
  -f FILE, --file FILE  Input file for mutation primer design, each line
                        stands for one mutation.
  -p, --primerFile      Amplicon sequence file for mutation primer design,
                        ignore sequence extraction.
  -g GENOME, --genome GENOME
                        Human genome sequence file (from NCBI).
  -G GFFFILE, --gffFile GFFFILE
                        Human genome annotation file (from NCBI).
  -k NAMEDB, --nameDB NAMEDB
                        Name convertion database file .
  -o OUTFILE, --outFile OUTFILE
                        Output file for amplicon sequences.
  -m MAXLENGTH, --maxLength MAXLENGTH
                        Amplicon maximum length. (Default:400bp)
  -n PRIMERNUM, --primerNum PRIMERNUM
                        Number of primer pairs for each amplicon. (Default:1)
  -s GENESOURCE, --geneSource GENESOURCE
                        Type of genenames.
                        (Default:gene_symbol,[ensembl_id],[hgnc_id])
  --primer_opt_size PRIMER_OPT_SIZE
                        Primer optimize size. (Default:20)
  --primer_min_size PRIMER_MIN_SIZE
                        Primer minimum size. (Default:18)
  --primer_max_size PRIMER_MAX_SIZE
                        Primer maximum size. (Default:27)
  --primer_opt_tm PRIMER_OPT_TM
                        Primer optimize temperature. (Default:60)
  --primer_min_tm PRIMER_MIN_TM
                        Primer minimum temperature. (Default:57)
  --primer_max_tm PRIMER_MAX_TM
                        Primer maximum temperature. (Default:63)
  --primer_min_gc PRIMER_MIN_GC
                        Primer minimum GC. (Default:20)
  --primer_max_gc PRIMER_MAX_GC
                        Primer maximum GC. (Default:80)
  --primer_max_poly_num PRIMER_MAX_POLY_NUM
                        Primer maximum ploy base number. (Default:5)
  
```
## Examples
```
1. Design 5 primers for 2 bases of mutation in gene TP53, with which 
                        located in 345 base pair in the cds sequence.
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -i "TP53:345:2" --primerNum 5
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -i "ENSG00000141510:345:2" --primerNum 5 -s ensembl_id
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -i "11998:345:2" --primerNum 5 -s hgnc_id
                        
2. Design primers for amplicons of three genes, including TPP1, CLN6 and WAS.
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -i "TPP1:1417:1,CLN6:794:4,WAS:223:5"
                        
3. Design primers for amplicons of genes from file.
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -f "samples.file"
4. Redesign primers for amplicons of genes which saved in temp.fa file, 
                        ignoring sequence extraction, with modified 
                        minimum acceptable melting temperature 50C.
   python3 MutPrimerDesign.py -g genome.fa -G gene.gff -k namedb -o temp.fa \
                        -p --primer_min_tm 50

```

If you use SCSA for your research, please kindly cite the following paper:

CAO Yinghao,PENG Gongxin. MutPrimerDesign: Design primers for human gene mutations located in coding 
                          sequence region[J].Chinese Journal of Bioinformatics,2020,18(3):169-175.
