# MutPrimerDesign: Primer design for human mutations

MutPrimerDesign is maintained by Yinghao Cao [yhcao@ibms.pumc.edu.cn].
## Download and Installation
```
git clone https://github.com/bioinfo-ibms-pumc/MutPrimerDesign.git
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
                        example: TPP1:1417:1,CLN6:794:4,WAS:223:5).
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
1. 
```
