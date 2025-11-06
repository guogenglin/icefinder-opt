# ICEfinder-opt
ICEfinder-opt is an optimized and enhanced local version of ICEfinder2, a tool for identifying integrative and conjugative elements (ICEs) in bacterial genomes.

The updated version includes multiple bug fixes, improved accuracy in boundary prediction, and additional support for batch processing and recognizing draft genomes.

# Preface
This repository provides an improved local version of ICEfinder2, originally developed by Meng Wang and Hong-Yu Ou. The present release includes bug fixes, revised output formatting, and several additional functions to enhance usability and detection accuracy.

For more information about the original tool, please visit:
🔗 https://tool2-mml.sjtu.edu.cn/ICEberg3/ICEfinder.php

For detailed inquiries regarding ICEfinder2, please contact the original authors.

However, if you have any questions or feedback regarding this modified version, you are also welcome to contact me.

# Major Improvements

This improved version introduces several key modifications compared to the original ICEfinder2, focusing on stability, functionality, and broader genome compatibility.

The main improvements are summarized below:

## Bug fixes and code optimization
Several bugs that previously caused runtime errors have been fixed, allowing the program to run smoothly on modern systems.

In addition, parts of the code were refactored to improve readability, maintainability, and logical consistency.

## Adjustment of matching length for tRNA search
In the original ICEfinder2, the tRNA search range was too short to accurately identify ICEs in many cases.

For example, when tested on a known sequence, ICEfinder1 successfully detected an ICE, while the web-based ICEfinder2 failed to do so due to the tRNA matching window being limited to only five genes.

In this improved version, the search range has been extended to 20 genes (for now), making the detection results more accurate and robust.

## Batch processing support
The local version of ICEfinder2 could only process one sample at a time.

The new version introduces batch mode, allowing users to analyze multiple genomes simultaneously, greatly improving efficiency for large-scale datasets.

## Support for draft and complex genomes
The original ICEfinder2 only handled chromosome-level assemblies and failed to process complete genomes with plasmids or draft genomes.

In this version, the code has been redesigned to support draft genomes and multi-contig assemblies, ensuring compatibility with a wider range of genomic data.

This improvement also fixes a critical issue where ICEs crossing the origin of the genome could not be detected in the original version.

## Improved output system
The previous version produced a JSON-based gene map by default.

In this revision, this visualization function is optional (not default), considering that batch users are usually more interested in summary statistics than individual maps.

Two new output tables — a summary table and a detailed table — are added to provide comprehensive information on detected MGEs (Mobile Genetic Elements) and their genomic distribution, facilitating downstream statistical analyses.

# Tips and Notes

## Input format recommendation
Using a GENBANK file as input is much faster than using FASTA.

When using FASTA, the file needs to be annotated to generate GFF and FAA files (e.g., via Prokka), which may take several minutes.

In the future, using Prodigal alone may further accelerate this process.

## Scope of this version
The original ICEfinder2 includes two functionalities: predicting MGEs in single bacterial genomes and in metagenomic datasets.

In this version, the focus is on single-genome analysis.

The original metagenome scripts are still included, but they have only minimal modifications to maintain the original functionality.

## Virulence and antimicrobial resistance gene prediction
The original tool used the BLAST databases provided by the authors.

In this version, ABRicate is used instead, making the analysis more flexible and easier to update with other databases.

## Unused databases
Some of the databases provided by the original authors were not used from the beginning, but they are still included for completeness.

## Limitations with draft genomes
ICE predictions in draft genomes may be biased. Experimental validation or additional sequencing is recommended to confirm results.

This is because ICE detection is pattern-based, not sequence-alignment-based.

If an ICE is split across multiple contigs (common for large ICEs, tens of kb in length, often composed of integrases and IS elements), predictions in the middle of ICEs may be less accurate.

To indicate confidence, genes essential for ICE structure are marked as:
`H` = high probability
`L` = low probability
Note that this marking system cannot fully eliminate bias.

The ICE predicted that located in 2 and more contigs will be marked a `r` as the prefix, for example, aHPS7_rICE1.

## Visualization limitations
Genetic maps of ICEs that are split across contigs will not be drawn, but antimicrobial resistance genes within these ICEs will still be predicted.

## Output changes
The FAA format output has been removed; only FASTA files are retained.

## Code quality and future improvements
You may notice some duplicate code and complex, clunky logic in the scripts.

Many areas could be improved for efficiency and readability.

Contributions to make the code faster or cleaner are welcome!

# Dependencies
`BLAST+`, `prodigal`, `prokka`, `macsyfinder`, `hmmsearch`, `vmatch`, `abricate`, `kraken`, `defense-finder`, `biopython`, `python` >= 3.8

# Usage
```
icefinder-opt [-h] -i INPUT [INPUT ...] [-t TYPE] [-o OUTPUT] [-c THREADS] [-j] [-v]
Input and Output:
  -i, --input             FASTA/Genbank format file, Genbank format file accepted only for single genome.
  -t, --type              Genome Type: Single/Metagenome (default: single)
  -o, --output            Output dir (default: ICEfinder_result)
Parameters:
  -c, --threads           Threads to use for BLAST searches (default: 4)
  -j, --json              output the json based genetic map (default: False)
  -h, --help              show this help message and exit
  -v, --version           Show version number and exit
```

## Quick start
``` Python
python icefinder-opt.py -i *.fasta 
```

# Output
The output file for every input genome will include a fasta file for every predicted MGE.

example fasta format file for MGEs:
```fasta
>04-14752_NZ_JAWDGC010000017.1_71342-106075
TTAGTCGCGGTTGGTGGATGCTTGTTCCAGTTCATCCAGAGTACGGAAGAGAAACATTTC
CCGGTTTTTCGCATCAAACGGAGCGGGAGTACTGTTTCTCATCACTGTCCTGATGCTGCC
AAACAGCCTGACCAGCGATTCTTTCGTTTTATCATCCGGCATTACATACATGGCATAGTG
CCAGCGTCCGGTGTCGCTGGCGGCCAAATGACTGTTAATAATGCTGATGTAGCGGGCGCG
GGTTTTCAGCGAGCGCTCTGTTTCCACCGCAACGATGGCACCGCTGCTTAAGGTAATAAT
GCCGTCCGGTCGGTGCCGGACACCCGGATATCGCGCCATAAACGCCCCCCGGTCGCCGTT
TAGCCAGCCTGTACCGCCTTTTTCCTCCAGGGCCAGCCTGACCCGCTGGTTTAACAGACG
ATGCTCCAGCGTCCAGTGCCTGAGTTTTCCCGGCTCAAAATAAGCGGGAAAAACAGGATC
ATCCGCGTTAACCACCCTTGCCAGTCCCGGCATAGTTATTCCCCACAATGATTTTTTCCC
GGTCAGTACCGGATATTCATGTTTACTCAACAAGCCCAGCGCAACCGCTTTATTCAGCAG
GTTATACATCGCATGGTTATGTTTTCCTTTATAGCCGACAACTTTTTTTAATGTCTGGAA
```
tab-separated output file contain details of every predicted MGE:
Column | Example | Description
-------|---------|------------
Isolate | `aHPS7` | The name of the inputfile
MGE | `aHPS7_ICE1` | The name of the predicted MGE
Location | `NZ_CP049090.1:1353271..1421762` | The location of the predicted MGE
Length | `68492` | Length of the predicted MGE
GC | `37.86` | GC content of the predicted MGE (%)
Relaxase_Type | `MOBH` | The classification of the relaxase enzyme
Systems | `typeG` | The classification of the conjugation or secretion machinery
oriT | `-` | AThe origin of transfer of the predicted MGE
attL | `1353271..1353285` | The left attachment site
attR | `1421748..1421762` | the right attachment site
att_seq | `CGGATTTTGAATCCG` | The sequence of the attachment site
tRNA | `tRNA-Leu,NZ_CP049090.1:1421704..1421790[-]` | The host tRNA gene targeted by the ICE for integration
AMR | `sul2;aph(3'')-Ib;aph(6)-Id;aph(3')-Ia;blaROB-1;tet(B)` | The antimicrobial resistance gene located in the MGE
middle | `-` | Genome contigs predicted to be located in the middle of an ICE.
middle_probability | `-` | Confidence for the prediction that a given middle region truly belongs to the ICE.

for example:
```
Isolate	MGE	Location	Length	GC	Relaxase_Type	Systems	oriT	attL	attR	att_seq	tRNA	AMR	middle	middle_probably
1082_KOXY	1082_KOXY_ICE1	NZ_JWDP01000035.1:(3757..24983)	21227	57.24	MOBP1	typeI	-	-	-	-	-	-	-	-
1082_KOXY	1082_KOXY_ICE2	NZ_JWDP01000172.1:(16845..43731)	26887	49.51	MOBH	typeG	-	16845..16859	43717..43731	AACCGTAGAAATACG	tRNA-Leu(caa),NZ_JWDP01000172.1:16804..16888[+]	-	-	-
1082_KOXY	1082_KOXY_IME2	NZ_JWDP01000026.1:(2219..10245)	8027	39.77	MOBQ	-	-	-	-	-	-	-	-	-
1084_KOXY	1084_KOXY_ICE1	NZ_JWDN01000042.1:(3757..24983)	21227	57.24	MOBP1	typeI	-	-	-	-	-	-	-	-
1084_KOXY	1084_KOXY_ICE2	NZ_JWDN01000183.1:(16816..39690)	22875	48.16	MOBH	typeG	-	16816..16830	39676..39690	AAATCGGTAGACGCA	tRNA-Leu(caa),NZ_JWDN01000183.1:16804..16888[+]	-	-	-
1085_KOXY	1085_KOXY_ICE1	NZ_JWDM01000043.1:(16801..47867)	31067	49.00	MOBH	typeG	-	16801..16815	47853..47867	TGGCGAAATCGGTAG	tRNA-Leu(caa),NZ_JWDM01000043.1:16794..16878[+]	-	-	-
1085_KOXY	1085_KOXY_IME1	NZ_JWDM01000060.1:(5435..13461)	8027	39.77	MOBQ	-	-	-	-	-	-	-	-	-
1085_KOXY	1085_KOXY_ICE2	NZ_JWDM01000157.1:(3747..24973)	21227	57.24	MOBP1	typeI	-	-	-	-	-	-	-	-
109680-17	109680-17_ICE1	NZ_VNMN01000028.1:(23715..55638)	31924	53.00	MOBF	typeF	-	-	-	-	-	-	-	-
```

tab-separated output file contain summary of every predicted MGE:
Column | Example | Description
-------|---------|------------
Isolate | `aHPS7` | The name of the inputfile
ICE(rICE_included) | `2` | The number of predicted ICE (include rICE)
ICE(rICE_excluded) | `1` | The number of predicted ICE (not include rICE)
IME | `1` | The number of predicted IME
AICE | `0` | The number of predicted AICE
rICE | `1` | The number of predicted rICE

for example:
```
Isolate	ICE(rICE_included)	ICE(rICE_excluded)	IME	AICE	rICE
1082_KOXY	2	2	1	0	0
1084_KOXY	2	2	0	0	0
1085_KOXY	2	2	1	0	0
109680-17	1	1	0	0	0
```

optional json file and genetic map:
gene json:
```json
{
        "gene": "G5S31_RS06505",
        "pos": "1353443..1354310 [-], 868",
        "prod": "tyrosine-type recombinase/integrase, Phage_integrase",
        "featu": "Integrase"
    },
    {
        "gene": "G5S31_RS06510",
        "pos": "1354337..1356263 [-], 1927",
        "prod": "MobH family relaxase, MOBH",
        "featu": "Relaxase"
    },
    {
        "gene": "G5S31_RS06515",
        "pos": "1357033..1357933 [-], 901",
        "prod": "tyrosine-type recombinase/integrase, Phage_integrase",
        "featu": "Integrase"
    }
```
MGE information json:
```json
{
    "Type": "T4SS-type ICE",
    "Location (nt)": "1353271..1421762",
    "Length (bp)": "68492",
    "GC Content (%)": "37.86",
    "oriT seq": "-",
    "DRs": "attL:1353271..1353285(CGGATTTTGAATCCG)  attR:1421748..1421762(CGGATTTTGAATCCG)",
    "Relaxase Type": "MOBH",
    "Mating pair formation systems": "typeG",
    "Close to tRNA": "tRNA-Leu (1421704..1421790) [-]"
}
```
summary json:
```json
    {
        "region": "Region1",
        "location": "1353271..1421762",
        "length": "68492",
        "gc": "37.86",
        "type": "T4SS-type ICE",
        "detail": "aHPS7_ICE1"
    }

```
and a genetic map:
<img width="1505" height="751" alt="屏幕截图 2025-11-06 161732" src="https://github.com/user-attachments/assets/0c658726-9dd8-4e15-bea4-bf9dedc1a2cc" />

# Citation
If you use or publish results obtained with this tool, please prioritize citing the paper of the original authors:

M. Wang, G Liu, M Liu, C Tai, Z. Deng, J Song, H.Y. Ou. ICEberg 3.0: functional categorization and analysis of the integrative and conjugative elements in bacteria. Nucleic Acids Research.2024 Jan 5; 52(D1):D732-D737. - [doi: 10.1093/nar/gkad935](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767825/)

If you are also willing to cite this GitHub repository, it would be greatly appreciated: 

* Genglin G, *icefinder-opt*, **Github** `https://github.com/guogenglin/icefinder-opt`
