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
import gzip

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
    normally just the TLEN field. For SE reads this is the observed coverage of
    the genome (infered from CIGAR)
    """
    if abs(read.template_length) > 0:
        return abs(read.template_length)
    else:
        return(read.infer_read_length())


if __name__ == "__main__":

    # Init variables
    frag_counter = 0
   
    # Reads args
    parser = argparse.ArgumentParser(prog='bamTofrag.py', description="Transform a BAM file to a fragment file")
    parser.add_argument('-i', '--input', help="BAM file with barcode tag, sorted by read nmaes", required=True)
    parser.add_argument('-o', '--output', help="Ouptut fragment file (BED)", required=True)
    parser.add_argument('-t', '--tag', help="Barcode Tag. Default: RG", default="RG", type=str)
    parser.add_argument('-s', '--se', help="Report single-end reads", action="store_true")
    parser.add_argument('-sz', '--seisize', help="Single-end reads insert size. By default, the infered read length is used", default=0, type=int)
    parser.add_argument('-z', '--gzip', help="Compress the output file", action="store_true")
    parser.add_argument('-v', '--verbose', help="", action="store_true")
    parser.add_argument('-d', '--debug', help="", action="store_true")
 
    args = parser.parse_args()
         
    # Verbose mode
    if args.verbose:
        print("------")
        print("## bamToFrag.py")
        print("## input =", args.input)
        print("## output =", args.output)
        print("## tag =", args.tag)
        print("## single-end =", args.se)
        print("## single-end isize =", args.seisize)
        print("## verbose =", args.verbose)
        print("-----")

    # I/O
    samfile = pysam.Samfile(args.input, "rb")
    if args.gzip:
        ofile = gzip.open(re.sub(".gz$","",args.output) + ".gz", 'wb')
    else:
        ofile = open(args.output, 'w')

    for read in samfile.fetch(until_eof=True):
        frag_counter += 1

        if read.is_paired:
            if read.is_read1:
                name1 = read.query_name
                chrom1 = read.reference_name
                start1 = read.reference_start
                bc1 = get_read_tag(read, args.tag)
                isize = get_frag_len(read)
            else:
                name2 = read.query_name
                chrom2 = read.reference_name
                start2 = read.reference_start
                bc2 = get_read_tag(read, args.tag)

                if chrom1 == chrom2 or bc1 != bc2:
                    s = min(start1, start2)
                    out = chrom1 + '\t' + str(s) + '\t' + str(s+isize) + '\t' + bc1 + '\t1\n'
                    if args.gzip:
                        ofile.write(out.encode('utf-8'))
                    else:
                        ofile.write(out)
                else:
                    if name1 != name2:
                        print >> sys.stderr, "## Error - reads not sorted by name"
                        sys.exit(-1)
                    else:
                        print >> sys.stderr, "## Reads " + read.query_name + " detected on different chromosomes - skipped"
        elif args.se:
            chrom1 = read.reference_name
            start1 = read.reference_start
            bc1 = get_read_tag(read, args.tag)
            if args.seisize > 0:
                isize = args.seisize
            else:
                isize = get_frag_len(read)
            out = chrom1 + "\t" + str(start1) + "\t" + str(s+isize) + "\t" + bc1 + "\t1\n"
            if args.gzip:
                ofile.write(out.encode('utf-8'))
            else:
                ofile.write(out)
        
        if (frag_counter % 1000000 == 0 and args.debug):
            stop_time = time.time()
            print( "##", frag_counter, stop_time-start_time )
            start_time = time.time()
            break

    samfile.close()
    ofile.close()
