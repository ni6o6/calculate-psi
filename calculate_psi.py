#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 13:11:31 2024

@author: Naoko Iida

Count junction reads supporting each splicing isoform from STAR SJ.out.tab
files and calculate the ratio (PSI) of each isoform. Supports N >= 2 isoforms.

Usage:
    python calculate_psi.py \
        -input_file sample_list.txt \
        -sj_file sj_list.txt \
        -output_prefix psi_result_

Input:
    -input_dir_list : tab-separated, no header
        <sample name>  <path to SJ.out.tab.gz (bgzip + tabix indexed)>
    -sj_file        : tab-separated, no header, 1 + 2N columns
        <gene>  <junctions of iso1> ... <junctions of isoN>  <name1> ... <nameN>
        junction format: chr:start-end (multiple junctions joined by ",")
        e.g. FGFR2  chr10:121515320-121517318,chr10:121517464-121519978  chr10:121515320-121518681,chr10:121518830-121519978  2c  2b

Output:
    <output_prefix><gene>.txt (one file per gene)
        Gene, Sample, Isoforms, <name>_count x N, SUM(isoforms), <name>_ratio x N
        ratio = count / SUM(isoforms), 3 decimal places

Requirement
pysam
"""



    
import os
import pysam


def calculate_psi(input_list, output_prefix, sj_file):  
    
    def parse_junction(sj):
        chrom, pos = sj.strip().split(":")
        start, end = pos.split("-")
        return chrom, int(start), int(end)
    
    def fmt_ratio(r):
        if r == 0 or r ==1:
            return(str(float(r)))
        else:
            v = round(r, 3)
            if v >= 1.0:
                v = 0.999
            elif  v <= 0.0:
                v = 0.001
            return f"{v:.3f}"
    
    def parse_sj_line(line):    
        """Parse one line of sj_file -> (gene, [[junctions of iso1], [junctions of iso2], [name of junction1]|...])"""
        cols = [c.strip() for c in line.rstrip("\n").split("\t")]
        cols = [c for c in cols if c != ""]
        
        if len(cols) < 5 or (len(cols) - 1) % 2 != 0:
            raise ValueError("sj_file line must have 1 + 2N columns (N >= 2): " + line.rstrip("\n").split("\t")[0])
        
        n_iso = (len(cols) - 1) // 2
        gene = cols[0]
        junction_sets = []
        names = []
        for i in range(n_iso):
            junction_sets.append([parse_junction(j) for j in cols[1 + i].split(",") if j.strip() != ""])
            names.append(cols[1 + n_iso + i])
        return gene, names, junction_sets
        
    file_list = []
    with open(input_list, 'r') as hin:
        for line in hin:
            cols = line.rstrip("\n").split("\t")
            file_list.append((cols[0], cols[1]))
            
    sj_entries = []
    with open(sj_file, "r") as sin:
        for line in sin:
            sj_entries.append(parse_sj_line(line))    
    
    for gene, names, junction_sets in sj_entries:
        output_file = output_prefix + gene + ".txt"
        with open(output_file, "w") as hout:
            header = (["Gene", "Sample", "Isoforms"]
                      + [n + "_count" for n in names]
                      + ["SUM(isoforms)"]
                      + [n + "_ratio" for n in names])
            print("\t".join(header), file=hout)
            
            for sample, sj_path in file_list:
                print("Processing " + gene + ":"+ sample)
                
                if not os.path.isfile(sj_path): 
                    print("No file " + sj_path)
                    continue
                    
                junc_tb = pysam.TabixFile(sj_path)
                iso_counts = []
                for js in junction_sets:
                    read_count = 0
                    for j in js:
                        chrom, start, end = j
                        sj_fetch = junc_tb.fetch(region = chrom + ':' + str(start-3) + '-' + str(end+3))
                        if sj_fetch is not None:
                            for record_line in sj_fetch:
                                record = record_line.split('\t')
                                rj_start = int(record[1])
                                rj_end = int(record[2])
                                adj_s = int(start) - rj_start
                                adj_e = int(end) - rj_end
                                if adj_s == adj_e:
                                    read_count += int(record[6])                                  
                    iso_counts.append(read_count)
    
                sum_iso = sum(iso_counts)
                
                if sum_iso == 0:
                    ratios = [0] * len(iso_counts)
                else:
                    ratios = [c / sum_iso for c in iso_counts]
                    
                ratios_fmt = [fmt_ratio(r) for r in ratios]
                
                rec = ([gene, sample, "|".join(names)] + [str(c) for c in iso_counts] + [str(sum_iso)] + [str(r) for r in ratios_fmt])

                print("\t".join(rec), file=hout)
            
    


if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("-input_dir_list", metavar = "input_dir_list", default = None, type = str,
                            help = "list of sample name and path to SJ.out.tab (tabix indexed)")
    parser.add_argument("-output_prefix", metavar = "output_prefix", default = "psi_result_", type = str,
                            help = "output file prefix")
    parser.add_argument("-sj_file", metavar = "sj_file", default = None, type = str,
                            help = "sj list of isoforms (Gene, junctions x N, name x N)")

    o = parser.parse_args()
    calculate_psi(o.input_dir_list, o.output_prefix, o.sj_file)
