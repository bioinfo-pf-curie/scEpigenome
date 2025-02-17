#!/usr/bin/env python
"""
Transform a BAM file to a count table based on barcode information and genomic bins
"""

import time
import argparse
import sys
import os
import re
import math
import pysam
from collections import defaultdict

def timing(function, *args):
    """                              
    Run a fonction and return the run time and the result of the function
    If the function requires arguments, those can be passed in
    """
    startTime = time.time()
    result = function(*args)
    print('%s function took %0.3f ms' % (function.func_name, (time.time() - startTime) * 1000))
    return result

def get_read_tag(read, tag):
    """
    Extract a flag from a read alignment
    """
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return None


def get_frag_len(read):
    """
    Get the observed template length of a read. For a paired-end read, this is
    normally just the TLEN field (i.e difference between the two read starts). 
    For SE reads this is the observed coverage of the genome (infered from CIGAR)
    """
    if abs(read.template_length) > 0:
        return abs(read.template_length)
    else:
        return(read.infer_read_length())


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    read_pairs = defaultdict(lambda: [0, 0])
    for read in bam.fetch(until_eof=True, region=region_string):
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.is_paired:
            qname = read.query_name
            yield read, None
        else:
            qname = read.query_name
            # catch orphan R2  
            if read.is_read1:
                read_pairs[read.query_name][0] += 1
            elif read.is_read2:
                read_pairs[read.query_name][1] += 1

            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                elif read.is_read2:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                elif read.is_read2:
                    yield read_dict[qname][0], read
                del read_dict[qname]

    # once the bam has been read check for read2 without read1 and return it as singleton 
    for name, (read1_count, read2_count) in read_pairs.items():
        if read2_count > 0 and read1_count == 0:
            yield read_dict[name][1], None
            del read_dict[name]



if __name__ == "__main__":

    # Init variables
    frag_counter = 0
    pair_counter = 0
    single_counter = 0


    # Reads args
    parser = argparse.ArgumentParser(prog='bamTofrag.py', description="Transform a BAM file to a fragment file")
    parser.add_argument('-i', '--input', help="BAM file with barcode tag, sorted by read names", required=True)
    parser.add_argument('-o', '--output', help="Output fragment file (BED)", default='stdout')
    parser.add_argument('-t', '--tag', help="Barcode Tag. Default: XB", default="XB", type=str)
    parser.add_argument('-s', '--se', help="Report single-end reads", action="store_true")
    parser.add_argument('-sz', '--seisize', help="Single-end reads insert size. By default, the inferred read length is used", default=0, type=int)
    parser.add_argument('-v', '--verbose', help="", action="store_true")
    parser.add_argument('-d', '--debug', help="", action="store_true")
 
    args = parser.parse_args()
         
    # Verbose mode
    if args.verbose:
        print("------", file=sys.stderr)
        print("## bamToFrag.py", file=sys.stderr)
        print("## input = " + args.input, file=sys.stderr)
        print("## output = " + args.output, file=sys.stderr)
        print("## tag = " + args.tag, file=sys.stderr)
        print("## report single-end = " + str(args.se), file=sys.stderr)
        print("## single-end extension = " + str(args.seisize), file=sys.stderr)

    # I/O
    samfile = pysam.AlignmentFile(args.input, 'rb')
    if args.output != "stdout":
        ofile = open(args.output, 'w')
    else:
        ofile = sys.stdout
    
    start_time = time.time()

    for read1, read2 in read_pair_generator(samfile):
        
        frag_counter += 1
        name1 = read1.query_name
        chrom1 = read1.reference_name
        start1 = read1.reference_start
        bc1 = get_read_tag(read1, args.tag)
        isize = get_frag_len(read1)

        if read2 is not None:
            pair_counter += 1
            chrom2 = read2.next_reference_name
            start2 = read2.next_reference_start
            if chrom1 == chrom2:
                s = min(start1, start2)
                out = chrom1 + '\t' + str(s) + '\t' + str(s+isize) + '\t' + bc1 + '\t1\n' 
                ofile.write(out)
            else:
                print("Warning - reads [" + name1 + "] mapped on different chromosomes", file=sys.stderr)
        elif args.se:
            single_counter += 1
            if args.seisize > 0:
                isize = args.seisize
            else:
                isize = get_frag_len(read1)
            out = chrom1 + "\t" + str(start1) + "\t" + str(start1 + isize) + "\t" + bc1 + "\t1\n"
            ofile.write(out)
        else:
            print("Warning - reads [" + name1 + "] discarded. Use '--se' to report singleton reads", file=sys.stderr)

        if (frag_counter % 1000000 == 0 and args.debug):
            stop_time = time.time()
            print("##", frag_counter, stop_time - start_time)
            break

    if args.verbose:
        print("## Processed Fragment = " + str(frag_counter), file=sys.stderr)
        print("## Reported Pairs = " + str(pair_counter), file=sys.stderr)
        print("## Reported Singletons = " + str(single_counter), file=sys.stderr)


    samfile.close()
    ofile.close()