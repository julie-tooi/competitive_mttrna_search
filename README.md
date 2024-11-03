# Competitive mt-tRNA annotation

mt-tRNA annotation based on cmsearch results

## Scripts description

### Start with annotation

_mitoannotation.py_ small pipeline for the mitogenomes 
annotation based on using cmsearch with custom CM set.

Important parameters: --max (turn off filtering), -g (global search mode)


**Input data**

```
  -i INPUT, --input INPUT
                        Path to mitogenome in fasta format
  -o OUTPUT, --output OUTPUT
                        Path to cmsearch output (folder with the results can be created if needed)
  -c CMS, --cms CMS     Path to CMs set
  -t THREADS, --threads THREADS
                        Threads number for cmsearch
```

**Output data**

Directory with cmsearch output files for all the CMs

### Collect the annotation results

_collect_annotation_results.py_ collect the cmsearch results from one folder
to one file (for one folder - annotation only for one mitogenome expected)

**Input data**

```
  -i INPUT, --input INPUT
                        Path to folder with the annotation results
  -o OUTPUT, --output OUTPUT
                        Path to table with results (.tsv)
```

**Output data**

Table in _.tsv_ format with annotation information, such as tRNA 
sequence, structure, genome coordinates, bit-score and e-value of hit


### Competitive locus-specific filtering

_filter_results.py_ 

**Input data**

```
  -i INPUT, --input INPUT
                        Path to collected annotation (table in .tsv format)
  -o OUTPUT, --output OUTPUT
                        Path to updated and filtered annotation
  -c, --codon           Require codon table: strict codon filter
  -t TRN, --trn TRN     Path to trna-codon table (NCBI format)
  -f FILTER_THRESHOLD, --filter_threshold FILTER_THRESHOLD
                        Bit-score filter: set the filter threshold
```

**Output data**

Table 1 in _.tsv_ format with extended annotation information. Info
about anticodon loop added.
Table 2 in _.tsv_ format: filtered information with locus grouping 
and filtration


***
## Requirements

python ≥ 3.6

infernal ≥ 1.1


## Usage

CMs for mt-tRNAs can be found, for example, [here](https://zenodo.org/records/2672835)

Usage sample:

```commandline
python mitoannotation.py -i test_files/NC_010214.1.fasta -o test_files/NC_010214.1 -c test_files/cms -t 6

python collect_annotation_results.py -i test_files/NC_010214.1 -o test_files/NC_010214.1_annotation_results.tsv

python filter_results.py -i test_files/NC_010214.1_annotation_results.tsv -o test_files/NC_010214.1_annotation_results_filtered.tsv -t test_files/transl_table2.txt -f 20 -c

```


***


## Author

Julie Ozerova [julie@bioinf.uni-leipzig.de](julie@bioinf.uni-leipzig.de)
