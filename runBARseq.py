#!/usr/bin/env python

import sys
import os
import argparse
import logging
import pandas as pd
import pylab as plt
import numpy as np
import editdistance
import networkx as nx
import math
import gzip
from textwrap import wrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker as lm

from multiprocessing.pool import Pool
from functools import partial

#########################
### Utility Functions ###
#########################

## Parse FASTQ entry
def processFQ(lines = None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

## Match barcode sequence based on template
def match_bc(seq, template):
  count = 0
  degenerate = {'A' : 'A', 'C' : 'C', 'G' : 'G', 'T' : 'T',
                'W' : 'AT', 'S' : 'CG', 'M' : 'AC', 'K' : 'GT', 'R' : 'AG', 'Y' : 'CT',
                'B' : 'CGT', 'D' : 'AGT', 'H' : 'ACT', 'V' : 'ACG', 'N':'ACGT'}
  for x in range(len(seq)):
    if not seq[x] in degenerate[template[x]]:
      count += 1
  return count

## Functions to retrived matrix cell position from 1-d array
def calc_row_idx(k, n):
    return int(math.ceil((1/2.) * (- (-8*k + 4 *n**2 -4*n - 7)**0.5 + 2*n -1) - 1))

def elem_in_i_rows(i, n):
    return i * (n - 1 - i) + (i*(i + 1))/2

def calc_col_idx(k, i, n):
    return int(n - elem_in_i_rows(i + 1, n) + k)

def condensed_to_square(k, n):
    i = calc_row_idx(k, n)
    j = calc_col_idx(k, i, n)
    return i, j

#######################################
## Save BAR-seq outputs for a sample ##
## Prepare TSV files and plots       ##
#######################################
def save_barseq_outputs(npa, bcode_merge, bcode_ids, bcode_len, sample_lab, out_dir, ed_thr, min_count, saturation, sum_counts):
    full_fasta_name = out_dir + "/" + sample_lab + ".barcode.mc" + str(min_count) + ".fa.gz"
    with gzip.open(full_fasta_name, "wt") as fasta_out:
        for barcode in bcode_merge.keys():
            for id in bcode_ids[barcode]:
                fasta_out.write(">" + id + "\n" + barcode + "\n")
    full_out_name = out_dir + "/" + sample_lab + ".barcode.mc" + str(min_count) + ".tsv"
    full_out = open(full_out_name, "w")
    full_out.write("\t".join([sample_lab + "-" + x for x in ["Barcode", "Count", "SaturationPerc"]]) + "\n")
    sel_out_name = out_dir + "/" + sample_lab + ".barcode.mc" + str(min_count) + ".sat" + str(saturation) + ".tsv"
    sel_out = open(sel_out_name, "w")
    sel_out.write("\t".join([sample_lab + "-" + x for x in ["Barcode", "Count", "SaturationPerc"]]) + "\n")
    num_selected = 0
    whole_bcode_compos = [{'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0} for k in range(bcode_len)]
    num_selected = 0
    sum_counts_cumul = 0
    ss = 0
    bcode_counts = []
    bcode_saturations = []
    stop_selection = False
    for npa_elem in npa:
        sum_counts_cumul += npa_elem['count']
        ss = (sum_counts_cumul / sum_counts) * 100
        logging.debug("%s\t%d (%d)\t%.2f", npa_elem['barcode'], npa_elem['count'], sum_counts_cumul, ss)
        full_out.write("\"" + npa_elem['barcode'] + "\"\t" + str(npa_elem['count']) + "\t" + "{:.2f}".format(round(ss, 2)) + "\n")
        bcode_counts.append(npa_elem['count'])
        bcode_saturations.append(ss)
        if not stop_selection:
            sel_out.write("\"" + npa_elem['barcode'] + "\"\t" + str(npa_elem['count']) + "\t" + "{:.2f}".format(round(ss, 2)) + "\n")
            num_selected += 1
        if ss > saturation:
            stop_selection = True
        for p in range(bcode_len):
            if npa_elem['barcode'][p] in whole_bcode_compos[p]:
                whole_bcode_compos[p][npa_elem['barcode'][p]] += npa_elem['count']
    logging.info("Barcodes selected at saturation %d: %d", saturation, num_selected)
    # Prepare plots
    logo_df = pd.DataFrame(whole_bcode_compos)
    logo = lm.Logo(logo_df, font_name = 'DejaVu Sans')
    plt.savefig(out_dir + "/" + sample_lab + ".barcode.structure.png")
    #
    ranks = np.arange(1, len(npa) + 1, 1)
    hline = [saturation for x in range(len(npa))]
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('Rank')
    ax1.set_ylabel('Counts', color = color)
    ax1.plot(ranks, bcode_counts, color = color)
    ax1.tick_params(axis = 'y', labelcolor = color)
    plt.axvline(x = num_selected, linestyle = "--", color = color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel('Saturation (%)', color = color)  # we already handled the x-label with ax1
    ax2.plot(ranks, bcode_saturations, color = color)
    ax2.plot(ranks, hline, "--", color = color)
    ax2.tick_params(axis = 'y', labelcolor = color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(out_dir + "/" + sample_lab + ".barcode.saturation.png")
    return 0


#######################################################
## Filter barcodes based on (inferred) length        ##
## and on different filters: nofilter, percfilter,   ##
## and fixedstruct (which requires a IUPAC structure ##
#######################################################
def filter_barcodes(bcode_fq, bcode_lens, filter_struct, iupac_struct):
    bcode_len =  max(bcode_lens, key = bcode_lens.get)
    logging.info("Barcode length: %d", bcode_len)
    if filter_struct == "fixedstruct" and len(iupac_struct) < bcode_len:
        iupac_struct = iupac_struct + ''.join(["N" for x in range(bcode_len - len(iupac_struct))])
        logging.warning("IUPAC code %s is shorter than barcode length!")
        logging.warning("Adding N chracters: %s", iupac_struct)
    # Check barcode lengths and structures
    whole_bcode_compos = [{'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'sum' : 0} for k in range(bcode_len)]
    for bcode_seq in bcode_fq.keys():
        if len(bcode_seq) != bcode_len:
            continue
        for p in range(bcode_len):
            if bcode_seq[p] not in ['A', 'C', 'G', 'T']:
                continue
            whole_bcode_compos[p][bcode_seq[p]] += bcode_fq[bcode_seq]
            whole_bcode_compos[p]['sum'] += bcode_fq[bcode_seq]
    # Filter barcodes by length and structure
    bcode_filter = {}
    bcode_content = {}
    for bcode_seq in bcode_fq.keys():
        if len(bcode_seq) != bcode_len:
            continue
        # No structural filter
        if filter_struct == "nofilter":
            bcode_filter[bcode_seq] = bcode_fq[bcode_seq]
            bcode_content[bcode_seq] = [bcode_seq.count('A'), bcode_seq.count('C'),
                                        bcode_seq.count('T'), bcode_seq.count('T')]
        # Filter based on nucleotide percentage
        elif filter_struct == "percfilter":
            wrong_pos = 0
            for p in range(bcode_len):
                if bcode_seq[p] not in ['A', 'C', 'G', 'T']:
                    continue
                if (whole_bcode_compos[p][bcode_seq[p]] / whole_bcode_compos[p]['sum']) * 100 < 1:
                    wrong_pos += 1
            if wrong_pos == 0:
                bcode_filter[bcode_seq] = bcode_fq[bcode_seq]
                bcode_content[bcode_seq] = [bcode_seq.count('A'), bcode_seq.count('C'),
                                            bcode_seq.count('T'), bcode_seq.count('T')]
        # Filter based on fixed structure (IUPAC)
        elif filter_struct == "fixedstruct":
            if match_bc(bcode_seq, iupac_struct) == 0:
                bcode_filter[bcode_seq] = bcode_fq[bcode_seq]
                bcode_content[bcode_seq] = [bcode_seq.count('A'), bcode_seq.count('C'),
                                            bcode_seq.count('T'), bcode_seq.count('T')]
        # Should never reach this point
        else:
            logginng.error("Wrong structural filter")
            return 1
    logging.info("Barcodes after filter: %d", len(bcode_filter))
    return bcode_filter, bcode_len

#########################################################
## Compute saturation values for the computed barcodes ##
#########################################################
def compute_saturation(bcode_merge, bcode_len, min_count):
    dtype = [('barcode', np.unicode_, bcode_len), ('count', int)]
    vals = []
    tot_count = 0
    for barcode in bcode_merge.keys():
        count = bcode_merge[barcode]
        tot_count += count
        # Filter out barcodes having count less than min_count
        if int(count) >= min_count:
            vals.append((barcode, count))
    # Create numpy array
    npa = np.array(vals, dtype = dtype)
    logging.info("Sum of all barcode counts %d", tot_count)
    # Sort barcodes based on counts
    npa = np.sort(npa, order = 'count')[::-1]
    return npa

###############################################
## Compute edit distances among all the BARs ##
###############################################

# Compute edit distances among sequences (in a range of positions) from a linearized matrix
def computeED(seq_array, range_pos):
    num_seqs = len(seq_array)
    dm = np.zeros((num_seqs * (num_seqs - 1)) // 2, dtype = np.int8)
    for i in range(range_pos[0], range_pos[1]):
        x,y = condensed_to_square(i, num_seqs)
        dm[i] = editdistance.eval(seq_array[x], seq_array[y])
    return dm

# Split a range of values into 'num_split' sub-intevals
def splitRange(num_elem, num_split):
    int_size = int(num_elem / num_split)
    bsize = [int_size for i in range(0, num_split)]
    surplus = num_elem % num_split
    i = 0
    for ss in range(surplus, 0, -1):
        bsize[i] += 1
        i = (i+1) % num_split
    split_array = []
    k = 0
    for i in range(0, num_split):
        split_array.append([k, k + bsize[i]])
        k += bsize[i]
    return split_array

# Compute edit distances among all sequences
def compute_edit_dist(bc1, num_threads):
    logging.info("Computing edit distances among barcodes")
    m = len(bc1)
    dm = np.zeros((m * (m - 1)) // 2, dtype = np.int8)
    if num_threads == 1:
        mx = 0
        for bx in range(0, m - 1):
            for by in range(bx + 1, m):
                dm[mx] = editdistance.eval(bc1[bx], bc1[by])
                mx += 1
    else:
        split_array = splitRange(len(dm), num_threads)
        with Pool(num_threads) as pool:
            func = partial(computeED, bc1)
            for res in pool.map(func, split_array):
                dm += res
    return dm

##########################################
## Merge barcodes based on ego-networks ##
##########################################
def merge_barcodes(g, files, threshold, ids):
    seqs = np.array([g.nodes[x]['seq'] for x in g.nodes()])
    merged = {}
    if len(seqs) > 1:
        k_nodes = np.array(g.nodes())
        counts = np.array([g.nodes[x]['counts'] for x in g.nodes()])
        counts = pd.DataFrame(counts.reshape((len(counts), 1)), index = seqs, columns = files)
        max_s = pd.concat([pd.Series(k_nodes, index = seqs), np.max(counts, axis = 1)], axis = 1)
        max_s.columns = ['node', 'max_count']
        ordered_nodes = np.max(counts, axis = 1).sort_values(ascending = False)
        while 1:
            repr_seq = ordered_nodes.index[0]
            node = max_s.loc[ordered_nodes.index[0], 'node']
            subr = nx.ego_graph(g, node, threshold, distance = 'weight').nodes()
            included_nodes = max_s.index[[x in subr for x in max_s.node]]
            merged[repr_seq] = np.sum(counts.loc[max_s.index[[x in subr for x in max_s.node]]], axis=0)[0]
            excluded_nodes = [x for x in ordered_nodes.index if not x in included_nodes]
            for x in included_nodes:
                if x != repr_seq:
                    ids[repr_seq] += ids[x]
                    ids[x] = list()
            if len(excluded_nodes) == 0:
                break
            ordered_nodes = ordered_nodes.loc[excluded_nodes]
            counts = counts.loc[excluded_nodes]
            max_s = max_s.loc[excluded_nodes]
        return(merged)
    else:
        return({seqs[0] : g.nodes[list(g.nodes())[0]]['counts']})

##########################################
## Extract barcode sequences from input ##
## reads using TagDust software         ##
##########################################
def extract_barcodes(in_fq, out_dir, sample_lab, td_prog, td_opt, num_threads=1):
    # Check previous TagDust runs
    if os.path.exists(out_dir + "/" + sample_lab + ".fq"):
        os.remove(out_dir + "/" + sample_lab + ".fq")
    if os.path.exists(out_dir + "/" + sample_lab + "_un.fq"):
        os.remove(out_dir + "/" + sample_lab + "_un.fq")
    if os.path.exists(out_dir + "/" + sample_lab + "_logfile.txt"):
        os.remove(out_dir + "/" + sample_lab + "_logfile.txt")
    td_cmd = td_prog + " " + td_opt + " -t " + str(num_threads) + " -o " + out_dir + "/" + sample_lab + " " + in_fq
    exit;
    stream = os.popen(td_cmd)
    output = stream.read().rstrip()
    
    # Parse TagDust results
    td_res = out_dir + "/" + sample_lab + ".fq"
    if not os.path.exists(td_res):
        logging.error("TagDust results %s not found!", td_res)
        return 1
    bcode_fq = {} # Store BARs sequences as keys and count the different found ones
    bcode_lens = {} # Store BARs lengths as keys and count the different (possible) found ones
    bcode_ids = {} # Store BARs sequences as keysStore and FASTQ id of each BAR
    with open(td_res, 'r') as in_fq:
        lines = []
        for line in in_fq:
            lines.append(line.rstrip())
            if len(lines) == 4:
                record = processFQ(lines)
                bcode_seq = record['sequence']
                bcode_name = record['name']
                if len(bcode_seq) not in bcode_lens:
                    bcode_lens[len(bcode_seq)] = 0
                bcode_lens[len(bcode_seq)] += 1
                if bcode_seq not in bcode_fq:
                    bcode_fq[bcode_seq] = 0
                    bcode_ids[bcode_seq] = []
                bcode_fq[bcode_seq] += 1
                bcode_ids[bcode_seq].append(bcode_name)
                lines = []
    logging.info("Number of input read ids: %d", sum(len(bcode_ids[id]) for id in bcode_ids))
    logging.info("Read %d distinct barcodes from %s", len(bcode_fq), td_res)
    return bcode_fq, bcode_lens, bcode_ids


#################################################
## Perform BAR-seq analysis on a single sample ##
#################################################
def barseq(in_fq, sample_lab, out_dir, td_prog, td_opt, filter_struct, iupac_struct,
           ed_thr, min_count, saturation, num_threads=1, log_level=logging.INFO):
    logging.root.handlers = []
    logging.basicConfig(level = log_level,
                        format = '[%(asctime)s] %(levelname)-8s %(message)s',
                        datefmt = "%Y-%m-%d %H:%M:%S",
                        handlers = [logging.FileHandler(out_dir + "/" + sample_lab + "_barseq.log"),
                                    logging.StreamHandler()])
    logging.info("## Running BAR-seq on sample %s", sample_lab)

    # Extract barcodes
    bcode_seq,bcode_lens,bcode_ids = extract_barcodes(in_fq, out_dir, sample_lab, td_prog, td_opt, num_threads)

    # Filter barcodes by length and structure
    bcode_seq,bcode_len = filter_barcodes(bcode_seq, bcode_lens, filter_struct, iupac_struct)
    for bc in set(bcode_ids) - set(bcode_seq):
        del bcode_ids[bc]
    logging.debug("Number of read ids after filter: %d", sum(len(bcode_ids[id]) for id in bcode_seq))

    # Graph-based merging: compute edit distances
    bc1 = np.array(list(bcode_seq.keys()))
    dm = compute_edit_dist(bc1, num_threads)
    if sum(dm) == 0:
        logging.error("Error in computing edit distances among BARs.")
        return 1

    # Graph-based merging: build graph
    logging.info("Building graph")
    G = nx.Graph()
    l_bc = len(bc1)
    G.add_nodes_from(np.arange(l_bc))
    for x in G.nodes():
        G.nodes[x]['seq'] = bc1[x]
        G.nodes[x]['counts'] = bcode_seq[bc1[x]]
    e_mask = dm <= ed_thr
    e_idx = [condensed_to_square(x, l_bc) for x in np.where(e_mask)[0]]
    e_w = dm[e_mask]
    G.add_edges_from(e_idx)
    for x, e in enumerate(e_idx):
        G[e[0]][e[1]]['weight'] = e_w[x]

    logging.info("Merging barcodes with edit distance lower than %d", ed_thr)
    bcode_seq = {}
    for cc in nx.connected_components(G):
        # slow but...
        g = nx.subgraph(G, cc)
        bcode_seq.update(merge_barcodes(g, [sample_lab], ed_thr, bcode_ids))
    logging.info("Barcodes after merge: %d", len(bcode_seq))

    # Compute saturation
    npa = compute_saturation(bcode_seq, bcode_len, min_count)

    # Sum barcode counts
    sum_counts = float(np.sum(npa[:]['count']))
    logging.info("Barcodes having count greater than %d: %d", min_count, len(npa))
    logging.info("Sum of barcodes having count greater than %d: %.0f", min_count, sum_counts)
    
    # Print results
    save_barseq_outputs(npa, bcode_seq, bcode_ids, bcode_len, sample_lab, out_dir, ed_thr, min_count, saturation, sum_counts)
    logging.info("## BAR-seq on sample %s finished!", sample_lab)
    return 0


############
##  Main  ##
############
def main(raw_args = None):
    parser = argparse.ArgumentParser(prog = "runBARseq",
                                     description = "Run BAR-seq analysis pipeline.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input-files', help = "Comma separated list of files.",
                        required = True, dest = 'in_files')
    parser.add_argument('-tagdust-dir', '--tagdust-dir', help = "TagDust program directory.",
                        required = True, dest = 'tagdust_dir')
    parser.add_argument('-tagdust-opt', '--tagdust-opt', help = "TagDust run options (e.g. \"-1 P:CAGGG -2 R:N -3 S:TAGGGACAGGA -4 P:AAGCCTCCTCCT\").",
                        required = True, dest = 'tagdust_opt')
    parser.add_argument('-l', '--file-labels', help = "Comma separated list of file labels.",
                        dest = 'in_labels', default = "")
    parser.add_argument('-o', '--out-dir', help = "Output directory (must exist).",
                        dest = 'out_dir', default = ".")
    parser.add_argument('-f', '--filter-struct', help = "Structural filter: 'nofilter', 'percfilter', or 'fixedstruct'.",
                        dest = 'filter_struct', default = "percfilter")
    parser.add_argument('-iupac', '--iupac-struct', help = "IUPAC structure used only with 'fixedstruct' filter.",
                        dest = 'iupac_struct', default = "N")
    parser.add_argument('-e', '--edit-distance', help = "Edit distance threshold for graph-based merging.",
                        dest = 'ed_thr', type = int, default = 3)
    parser.add_argument('-s', '--saturation', help = "Saturation threshold.",
                        dest = 'saturation', default = 90, type = int)
    parser.add_argument('-c', '--min-count', help = "Minimum count threshold for an input barcode.",
                        dest = 'min_count', default = 3, type = int)
    parser.add_argument('-t', '--threads', help = "Num. of threads used in sequence comparisons.",
                        dest = 'threads', default = 1, type = int)
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

    logging.info("BAR-seq Started")
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
            sample_lab = sample_lab.replace(".gz", "")
            sample_lab = sample_lab.replace(".fq", "")
            sample_lab = sample_lab.replace(".fastq", "")
            input_labels.append(sample_lab)
    # Check filter
    if args.filter_struct not in ["nofilter", "percfilter", "fixedstruct"]:
        logging.error("Wrong structural filter: %s", args.filter_struct)
        logging.error("Allowed filters: 'nofilter', 'percfilter', or 'fixedstruct'")
        return 1
    if args.filter_struct == "fixedstruct" and args.iupac_struct == "":
        logging.warning("Empty IUPAC code reuqired for 'fixedstruct' filter")
        logging.warning("If you want to specify it use the -iupac option")
    # TagDust
    if not os.path.exists(args.tagdust_dir):
        logging.error("TagDust dir %s not found!", args.tagdust_dir)
        return 1
    td_prog = args.tagdust_dir + "/tagdust"
    if not os.path.exists(td_prog):
        td_prog = args.tagdust_dir + "/src/tagdust"
        if not  os.path.exists(td_prog):
            logging.error("TagDust program not found in %s", args.tagdust_dir)
            return 1
    # Run BAR-seq on all the samples
    for in_fq in range(0, len(input_files)):
        # Check input FQ file
        if not os.path.exists(input_files[in_fq]):
            logging.error("Input FASTQ file %s not found!", input_files[in_fq])
            return 1
        if barseq(input_files[in_fq], input_labels[in_fq],
                  args.out_dir, td_prog, args.tagdust_opt, args.filter_struct,
                  args.iupac_struct, args.ed_thr, args.min_count, args.saturation, args.threads, log_level) != 0:
            logging.error("Error in running BAR-seq on sample: %s", input_files[in_fq])
            return 1
    logging.root.handlers = []
    logging.basicConfig(level = log_level,
                        format = '[%(asctime)s] %(levelname)-8s %(message)s',
                        datefmt = "%Y-%m-%d %H:%M:%S",
                        handlers = [logging.StreamHandler()]
                        )
    # Compute Sharing
    if len(input_files) > 1:
        logging.info("-----------------------")
        logging.info("Compute Barcode Sharing")
        full_df = pd.DataFrame()
        for lab in input_labels:
            td_res = args.out_dir + "/" + lab + ".barcode.mc" + str(args.min_count) + ".tsv"
            if not os.path.exists(td_res):
                continue
            df = pd.read_csv(td_res, sep = '\t')
            df.columns = ["Barcode", lab, "Saturation"]
            df = df[["Barcode", lab]]
            if len(full_df.index) == 0:
                full_df = df
            else:
                full_df = pd.merge(full_df, df, on = 'Barcode', how = 'outer')
        full_df = full_df.fillna(0)
        full_df = full_df.sort_values(by = ['Barcode'])
        full_df.to_csv(args.out_dir + "/Barcode.mc" + str(args.min_count) + "_fullMatrix.tsv", index = False, header = True, sep = '\t', float_format='%.0f')

        # Selected res
        sel_df = pd.DataFrame()
        for lab in input_labels:
            td_sel_res = args.out_dir + "/" + lab + ".barcode.mc" + str(args.min_count) + ".sat" + str(args.saturation) + ".tsv"
            if not os.path.exists(td_sel_res):
                continue
            df = pd.read_csv(td_sel_res, sep = '\t')
            df.columns = ["Barcode", lab, "Saturation"]
            df = df[["Barcode", lab]]
            if len(sel_df.index) == 0:
                sel_df = df
            else:
                sel_df = pd.merge(sel_df, df, on = 'Barcode', how = 'outer')
        sel_df = sel_df.fillna(0)
        sel_df_hm = sel_df
        sel_df = sel_df.sort_values(by = ['Barcode'])
        sel_df.to_csv(args.out_dir + "/Barcode.mc" + str(args.min_count) + ".sat" + str(args.saturation) + "_fullMatrix.tsv", index = False, header = True, sep = '\t', float_format='%.0f')
        # Plot Heatmap
        sel_df_hm = sel_df_hm.set_index('Barcode')
        sel_df_hm = np.log(sel_df_hm + 1)
        fig, ax = plt.subplots(figsize=(10, 7))
        plt.pcolor(sel_df_hm, cmap = "Reds")
        #plt.imshow(sel_df_hm.values, cmap = "Reds")
        if len(sel_df_hm.index) < 100:
            plt.yticks(np.arange(0.5, len(sel_df_hm.index), 1), sel_df_hm.index)
        else:
            plt.yticks([])
        plt.xticks(np.arange(0.5, len(sel_df_hm.columns), 1), [ '\n'.join(wrap(l, round(80/len(sel_df_hm.columns)))) for l in sel_df_hm.columns])
        cbar = plt.colorbar()
        cbar.ax.set_title('ln(Count)')
        plt.tight_layout()
        plt.savefig(args.out_dir + "/Barcode.mc" + str(args.min_count) + ".sat" + str(args.saturation) + "_heatmap.png")
    logging.info("BAR-seq finished")
    return 0

if __name__ == "__main__":
    main()
