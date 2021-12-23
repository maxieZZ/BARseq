# BAR-Seq #

Software for the analysis of Homology-Directed Repair barcoded cells.

---
### Download and Installation ###

BAR-Seq is based on [TagDust](http://tagdust.sourceforge.net/) for the extraction of barcodes from NGS reads.
A copy of the TagDust software is present in the `ext-progs` folder.

To install TagDust follow the instructions provided [here](http://tagdust.sourceforge.net/#install).
The BAR-Seq pipeline was tested with TagDust v2.33.

### Requirements ###

BAR-Seq is based on python 3 (version 3.6 or higher) and requires the following packages:

* numpy (v1.18)
* editdistance (v0.5.3)
* networkx (v2.4)
* pandas (v1.0.5)
* matplotlib (v.3.2.2)
* logomaker (v0.8)
* scipy (v1.5.0)

---
### Run ###

The BAR-Seq pipeline was developed on a Linux 64bit server, but works also on Mac OS X systems.

Usage:

```
python runBARseq.py -i IN_FILES -tagdust-dir TAGDUST_DIR -tagdust-opt TAGDUST_OPT
                    [-l IN_LABELS] [-o OUT_DIR] [-f FILTER_STRUCT] [-iupac IUPAC_STRUCT]
                    [-e ED_THR] [-s SATURATION] [-c MIN_COUNT]
```

where:

```
  -i IN_FILES, --input-files IN_FILES
                        Comma separated list of files. (default: None)
  -tagdust-dir TAGDUST_DIR, --tagdust-dir TAGDUST_DIR
                        TagDust program directory. (default: None)
  -tagdust-opt TAGDUST_OPT, --tagdust-opt TAGDUST_OPT
                        TagDust run options (e.g. "-1 P:CAGGG -2 R:N -3 S:TAGGGACAGGA -4 P:AAGCCTCCTCCT"). (default: None)
  -l IN_LABELS, --file-labels IN_LABELS
                        Comma separated list of file labels. (default: )
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory (must exist). (default: .)
  -f FILTER_STRUCT, --filter-struct FILTER_STRUCT
                        Structural filter: 'nofilter', 'percfilter', or 'fixedstruct'. (default: percfilter)
  -iupac IUPAC_STRUCT, --iupac-struct IUPAC_STRUCT
                        IUPAC structure used only with 'fixedstruct' filter. (default: N)
  -e ED_THR, --edit-distance ED_THR
                        Edit distance threshold for graph-based merging. (default: 3)
  -s SATURATION, --saturation SATURATION
                        Saturation threshold. (default: 90)
  -c MIN_COUNT, --min-count MIN_COUNT
                        Minimum count threshold for an input barcode. (default: 3)
```

### Output ###

For each input file BAR-Seq produces the following outputs:

* `<SAMPLE_NAME>.barcode.mc<MIN_COUNT>.tsv` list of barcodes having count at least `<MIN_COUNT>`
* `<SAMPLE_NAME>.barcode.mc<MIN_COUNT>.sat<SATURATION>.tsv` list of barcodes having count at least `<MIN_COUNT>` cut to saturation `<SATURATION>`
* `<SAMPLE_NAME>.barcode.saturation.png` plot of counts and saturation level of barcodes
* `<SAMPLE_NAME>.barcode.structure.png` logo plot of the barcode compositions

where `<SAMPLE_NAME>` is the label corresponding to the sample or (if not provided) it is extracted from the input FASTQ file.

If more than one input file is provided, BAR-Seq computes also the sharing results:

* `Barcode.mc<MIN_COUNT>_fullMatrix.tsv` matrix containing all the barcodes having count at least `<MIN_COUNT>` (rows) in all the provided samples (columns)
* `Barcode.mc<MIN_COUNT>.sat<SATURATION>_fullMatrix.tsv` matrix containing all the barcodes having count at least `<MIN_COUNT>` cut to saturation `<SATURATION>` (rows) in all the provided samples (columns)
* `Barcode.mc<MIN_COUNT>.sat<SATURATION>_heatmap.png` heatmap plot of all the barcodes having count at least `<MIN_COUNT>` cut to saturation `<SATURATION>` (rows) in all the provided samples (columns)

---
The repository provides a procedure to assess the library complexity of barcodes.

Usage:

```
python runLibComplexity.py -i IN_FILES [-l IN_LABELS] [-o OUT_DIR] [-p CONFIDENCE] [-c MIN_COUNT]
```

where:

```
  -i IN_FILES, --input-files IN_FILES
                        Comma separated list of TSV files. (default: None)
  -l IN_LABELS, --file-labels IN_LABELS
                        Comma separated list of file labels. (default: )
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory (must exist). (default: .)
  -p CONFIDENCE, --perc-confidence CONFIDENCE
                        Confidence level used in the analysis. (default: 95)
  -c MIN_COUNT, --min-count MIN_COUNT
                        Minimum count threshold for an input barcode. (default: 1)
```

### Output ###

For each input file the runLibComplexity procedure produces the following outputs:

* `<SAMPLE_NAME>.barcode.mc<MIN_COUNT>_LibComplexity.png` scatter plot reporting for each barcode count (x-axis) the corresponding abundance (y-axis) 
* `<SAMPLE_NAME>.barcode.mc<MIN_COUNT>_LibComplexity.tsv` a two-column table containing for a sequence of increasing numbers of potential unique cells to track, the corresponding probability of tracking those cells (at confidence `<CONFIDENCE>`) with the provided library of barcodes having count at least `<MIN_COUNT>`, plus some additional indexes of library complexity (ENS, Equitability, Chao1 Index, and Richness).

where `<SAMPLE_NAME>` is the label corresponding to the sample or (if not provided) it is extracted from the input FASTQ file.

---
### Examples ###

In the `example` folder there are 3 input files provided as examples (`Sample1`, `Sample2`, and `Sample3`) for the BAR-Seq pipeline.

To run the BAR-Seq pipeline on `Sample1`, type:

```
python runBARseq.py -i example/Sample1.fastq.gz
                    -l Sample1
                    -tagdust-dir <TAGDUST_DIR>
                    -tagdust-opt "-1 P:CAGGGGATGCGGTGGGCTCTATGG -2 R:N -3 S:TAGGGACAGGA -4 P:TTGGTGACAGAAAAGCCCCATCCTTAGGCCTCCTCCTTCCTAGT"
                    -o example
```

where `<TAGDUST_DIR>` should be replaced with the path of installation of the TagDust software.

The runtime required for the computation on a standard server/workstation on `Sample1` is ~3 min.

The `example` folder contains also a file for the library complexity assessment (`LibExample.tsv`) which can be run as:

```
python runBARseq.py -i example/LibExample.tsv
                    -l LibExample
                    -o example
```

---
### Contacts ###
* Beretta Stefano
* Ivan Merelli