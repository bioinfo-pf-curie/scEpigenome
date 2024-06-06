#!/usr/bin/env python
# coding: utf-8

"""
Copyleft 2024 Institut Curie
Author(s): Nicolas Servant
Contact: nicolas.servant@curie.fr
This software is distributed without any guarantee under the terms of the CECILL License
See the LICENCE file for details

Add the barcode information in input files
/!\ Input and barcodes input files must be ordered in the same way, with the same number of reads
"""

import argparse
import sys
import os
import re
import gzip

def usage():
    """Usage function"""
    print ("Usage : python addBarcodeFlag.py")
    print ("-i/--input < file [fastq]>")
    print ("[-b/--barcode] <file with barcode info [txt]>")
    print ("[-o/--ofile] <output file [fastq]>")
    print ("[-v/--verbose] <Verbose>")
    print ("[-h/--help] <Help>")
    return

def args_parse():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--fastq1", required=True, help="Input file (.fastq.gz)")
    parser.add_argument("-r", "--fastq2", help="Input file 2 for paired-end data (.fastq.gz)")
    parser.add_argument("-b", "--barcode", required=True, help="Two columns text file wih read name and barcode (.txt)")
    parser.add_argument("-s", "--suffix", help="SUffix output file", default="_barcoded.fastq.gz")
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")

    args = parser.parse_args()
    return (args)


def add_bc_fastq_se(fastq_1, barcodes):
    """
    Add barcode in read names for single-end data
    """
    n_reads=0
    n_reads_barcoded=0

    ## outputs
    ofh = gzip.open('file_1.gz', 'wt')

    with gzip.open(fastq_1, 'rt') as fh, open(barcodes) as fhbc:
        while True:
            ## Get barcode info
            bcinfo = fhbc.readline()
            if len(bcinfo.strip()) == 0 :
                break
            else:
                bcrname = bcinfo.rstrip().split("\t")[0]
                bcname = bcinfo.rstrip().split("\t")[1]

            ## Fastq R1
            name = fh.readline().rstrip().split(" ")
            seq = fh.readline()
            ph = fh.readline()
            qual = fh.readline()

            if bcrname != name[0] :
                sys.stderr.write("Error in read names [" + bcrname + " / " + name[0] + " / " + "]\nAll files must be ordered in the same way with the same number of lines\n")
                sys.exit(-1)

            n_reads += 1
            if bcname != "None":
                n_reads_barcoded+=1
                ofh.write(name[0] + "_" + bcname + " " + ' '.join(name[1:]) + '\n')
                ofh.write(seq)
                ofh.write(ph)
                ofh.write(qual)

        fh.close()
        fh2.close()
        fhbc.close()
        ofh.close()
        ofh2.close()

        ## Logs
        sys.stderr.write("## addBarcode.py\n")
        sys.stderr.write("## Total reads = " + str(n_reads) + "\n")
        sys.stderr.write("## Reads with barcodes = " + str(n_reads_barcoded) + "\n")


def add_bc_fastq_pe(fastq_1, fastq_2, barcodes, osuffix):
    """
    Add barcode in read names for paired-end data
    """
    n_reads=0
    n_reads_barcoded=0

    ## outputs
    ofh = gzip.open(re.sub(".fastq.gz$", osuffix, fastq_1), 'wt')
    ofh2 = gzip.open(re.sub('.fastq.gz$', osuffix, fastq_2), 'wt')

    with gzip.open(fastq_1, 'rt') as fh, gzip.open(fastq_2, 'rt') as fh2, open(barcodes) as fhbc:
        while True:
            ## Get barcode info
            bcinfo = fhbc.readline()
            if len(bcinfo.strip()) == 0 :
                break
            else:
                bcrname = bcinfo.rstrip().split("\t")[0]
                bcname = bcinfo.rstrip().split("\t")[1]

            ## Fastq R1
            name = fh.readline().rstrip().split(" ")
            seq = fh.readline()
            ph = fh.readline()
            qual = fh.readline()
            
            ## Fastq R2
            name2 = fh2.readline().rstrip().split(" ")
            seq2 = fh2.readline()
            ph2 = fh2.readline()
            qual2 = fh2.readline()

            if bcrname != name[0] or bcrname != name2[0]:
                sys.stderr.write("Error in read names [" + bcrname + " / " + name[0] + " / " + name2[0] + "]\nAll files must be ordered in the same way with the same number of lines\n")
                sys.exit(-1)
            
            n_reads += 1
            if bcname != "None":
                n_reads_barcoded+=1
                ofh.write(name[0] + "_" + bcname + " " + ' '.join(name[1:]) + '\n')
                ofh.write(seq)
                ofh.write(ph)
                ofh.write(qual)

                ofh2.write(name2[0] + "_" + bcname + " " + ' '.join(name2[1:]) + "\n")
                ofh2.write(seq2)
                ofh2.write(ph2)
                ofh2.write(qual2)

        fh.close()
        fh2.close()
        fhbc.close()
        ofh.close()
        ofh2.close()

        ## Logs
        sys.stderr.write("## addBarcode.py\n")
        sys.stderr.write("## Total reads = " + str(n_reads) + "\n")
        sys.stderr.write("## Reads with barcodes = " + str(n_reads_barcoded) + "\n")

if __name__ == "__main__":

    args = args_parse()

    if args.fastq2:
        add_bc_fastq_pe(args.fastq1, args.fastq2, args.barcode, args.suffix)
    else:
        add_bc_fastq_pe(args.fastq1, args.barcode, args.suffix)
