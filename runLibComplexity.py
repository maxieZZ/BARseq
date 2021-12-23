#!/usr/bin/env python

import sys
import os
import argparse
import logging
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Sampling (with replacement)
def sample_wr(data, n_sample = 1000):
    idx = np.random.randint(0, len(data), n_sample)
    return data[idx]

# Number of events counted only once
def count_unique(data):
    U, C = np.unique(data, return_counts = True)
    return np.sum(C == 1)

# Compute Library complexity
def computeLibComplexity(out_dir, file_path, sample_name, conf_level, min_count, n_trials = 1000):
    barcodes = {}
    with open(file_path, "r") as tf:
        tf.readline()
        for line in tf.readlines():
            t = line.split()
            barcodes[t[0]] = int(t[1])
    v = np.array(list(barcodes.values()))
    v = v[v >= min_count]
    # Plot graph
    U, C = np.unique(v, return_counts=True)
    plt.scatter(U, C/np.sum(C), c = "blue")
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim([min(C/np.sum(C)) / 1.5, 1])
    plt.xlim([1, max(U) * 1.5])
    plt.title("Library Complexity")
    plt.xlabel("Barcode Count (log scale)")
    plt.ylabel("Normalized Count Abundance (log scale)")
    plt.savefig(out_dir + "/" + sample_name + "_LibComplexity.png")
    # Compute Indexes
    SSC = (np.sum(C) - 1) / np.sum(C)
    a1 = np.sum(v == min_count)
    a2 = np.sum(v == (min_count + 1))
    chao1 = np.round((len(v) + (a1 * (a1 - 1)) / ((a2 + 1) * 2)) * SSC) # Chao1 Index
    eL = np.exp(scipy.stats.entropy(v)) # ENS
    eq = (scipy.stats.entropy(v) / np.log(len(v))) # Equitability
    rich = len(v) / chao1
    # ---------    
    vpos = np.zeros(np.sum(v), dtype=int)
    s = 0
    for x in range(len(v)):
        for y in range(v[x]):
            vpos[s + y] = x
        s = s + v[x]
    un_cells = {10**x : 0 for x in range(0, 7)}
    for cells in un_cells.keys():
        sym_p = np.array([count_unique(sample_wr(vpos, cells)) / cells for x in range(n_trials)])
        un_cells[cells] = np.percentile(sym_p, int(conf_level))
    ENS = {"name" : "ENS", "value" : int(round(eL, 0))}
    Equit = {"name" : "Equitability", "value" : round(eq, 3)}
    Chao1 = {"name" : "Chao1_Index", "value" : int(round(chao1, 0))}
    Richness = {"name" : "Richness", "value" : round(rich, 5)}
    return {"idxs" : [ENS, Equit, Chao1, Richness] , "un_cells" : un_cells}

# Save results to TSV file
def saveRes(res, out_dir, sample_name):
    with open(out_dir + "/" + sample_name + "_LibComplexity.tsv", "w") as out:
        for idx in res["idxs"]:
            out.write("# " + idx["name"] + "=" + str(idx["value"]) + "\n")
        out.write("NumUniqueCells" + "\t" + "ProbTargetUniqueCells" + "\n")
        for n_cells in res["un_cells"]:
            out.write(str(n_cells) + "\t" + str(res["un_cells"][n_cells]) + "\n")
    return 0

def main(raw_args = None):
    parser = argparse.ArgumentParser(prog = "runLibComplexity",
                                     description = "Run analysis on Barcode Lib. Complexity.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input-files', help = "Comma separated list of TSV files.",
                        required = True, dest = 'in_files')
    parser.add_argument('-l', '--file-labels', help = "Comma separated list of file labels.",
                        dest = 'in_labels', default = "")
    parser.add_argument('-o', '--out-dir', help = "Output directory (must exist).",
                        dest = 'out_dir', default = ".")
    parser.add_argument('-p', '--perc-confidence', help = "Confidence level used in the analysis.",
                        dest = 'confidence', default = 95, type = int)
    parser.add_argument('-c', '--min-count', help = "Minimum count threshold for an input barcode.",
                        dest = 'min_count', default = 1, type = int)
    parser.add_argument('-v', '--verbose', help = "Increase verbosity.",
                        action = 'count', default = 0)
    args = parser.parse_args(raw_args)

    if args.verbose == 0:
        log_level = logging.INFO
    elif args.verbose == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level = log_level,
                        format = '[%(asctime)s] %(levelname)-8s %(message)s',
                        datefmt = "%Y-%m-%d %H:%M:%S")

    logging.info("LibComplexity Started")
    # Check output directory
    if not os.path.exists(args.out_dir):
        logging.error("Output directory %s not found!", args.out_dir)
        return 1
    # Check input files and labels
    input_files = args.in_files.split(",")
    input_labels = []
    if args.in_labels != "":
        input_labels = args.in_labels.split(",")
        if len(input_files) != len(input_labels):
            logging.error("Labels do not match input files!")
            return 1
    else:
        for in_fq_file in input_files:
            # Infer sample label
            sample_lab = in_fq_file.split("/")[-1]
            sample_lab = sample_lab.replace(".tsv", "")
            input_labels.append(sample_lab)
    # Check Confidence Level
    if args.confidence > 100 or args.confidence < 0:
        logging.error("Wrong confidence level!")
        logging.error("Only values between 0 and 100 allowed.")
        return 1
    # Check Min Count
    if args.min_count < 1:
        logging.error("Wrong min. count threshold!")
        logging.error("Only values greater than 0 allowed.")
        return 1
    # Run Lib Complexity on all the samples
    for in_tsv in range(0, len(input_files)):
        # Check input TSV file
        if not os.path.exists(input_files[in_tsv]):
            logging.error("Input TSV file %s not found!", input_files[in_tsv])
            return 1
        logging.info("Computing Lib Complexity on %s", input_labels[in_tsv])
        res = computeLibComplexity(args.out_dir, input_files[in_tsv], input_labels[in_tsv], args.confidence, args.min_count)
        saveRes(res, args.out_dir, input_labels[in_tsv])
    logging.info("LibComplexity finished")
    return 0

if __name__ == "__main__":
    main()
