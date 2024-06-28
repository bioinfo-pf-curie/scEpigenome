#!/usr/bin/env python
# coding: utf-8

"""
Mapping Research Pipeline                                                                                                                                                                                
Copyleft 2018 Institut Curie                                                                                                                                                                             
Author(s): Nicolas Servant                                                                                                                                                                               
Contact: nicolas.servant@curie.fr                                                                                                                                                                        
This software is distributed without any guarantee under the terms of the CECILL License                                                                                                                 
See the LICENCE file for details                                                                                                                                                                         

Remove duplicates for single-cell ChIP-seq data
"""

import getopt
import sys
import argparse
import os
import re
import pysam
import itertools 

def get_read_tag(read, tag):
    """
    Extract barcode from a read 
    """
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return None

def is_win_dup(read1, ref, distance):
    """ 
    True of two reads (R1/R2 are close to the reference
    """
    r1_start=None
    r2_start=None
    ref1_start=None
    ref2_start=None
    d1 = None
    d2 = None

    if read.query_name == ref.query_name:
        return False
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if ref.is_unmapped or ref.mate_is_unmapped:
        return False

    if read.is_read1:
        r1_start=read.reference_start
        r2_start=read.next_reference_start
    else:
        r2_start=read.reference_start
        r1_start=read.next_reference_start

    if ref.is_read1:
        ref1_start=ref.reference_start
        ref2_start=ref.next_reference_start
    else:
        ref2_start=ref.reference_start
        ref1_start=read.next_reference_start

    d1=abs(r1_start - ref1_start)
    d2=abs(r2_start - ref2_start)

    if d1 < distance and d2 < distance:
#        print("DISTANCE")
#        print(ref)
#        print(read)
        return True

    return False


def is_rt_dup(read1, ref):
    """
    R2 start is the same but not R1
    """
    if read.query_name == ref.query_name:
        return False
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if ref.is_unmapped or ref.mate_is_unmapped:
        return False

    if read.is_read1:
        r1_start=read.reference_start
        r2_start=read.next_reference_start
    else:
        r2_start=read.reference_start
        r1_start=read.next_reference_start

    if ref.is_read1:
        ref1_start=ref.reference_start
        ref2_start=ref.next_reference_start
    else:
        ref2_start=ref.reference_start
        ref1_start=read.next_reference_start

    if r2_start == ref2_start and r1_start != ref1_start:
        return True

    return False


if __name__ == "__main__":

    # Reads args
    parser = argparse.ArgumentParser(prog='rmDup.py', description="Remove RT and Window duplicates")
    parser.add_argument('-i', '--input', help="BAM file with barcode tag, sorted by read names", required=True)
    parser.add_argument('-o', '--output', help="Ouptut file (BAM)", required=True)
    parser.add_argument('-d', '--dist', help="Distance to consider reads as duplicates. Default: 150bp", default=150, type=int)
    parser.add_argument('-t', '--tag', help="Barcode tag", default="XB", type=str)
    parser.add_argument('-r', '--rt', help="Remove RT duplicates reads", action="store_true")
    parser.add_argument('-v', '--verbose', help="", action="store_true")
    args = parser.parse_args()                                                                                                                                                                             

    reads_counter = 0
    wdup_counter = 0
    rtdup_counter = 0
    pcrdup_counter = 0
    ref_barcode={}

    ## Read bam file
    samfile = pysam.Samfile(args.input, "rb")

    ## output
    if args.output == "-":
        outfile = pysam.AlignmentFile(args.output, "w", template=samfile)
    else:
        outfile = pysam.AlignmentFile(args.output, "wb", template=samfile)

    for read in samfile.fetch(until_eof=True):
        reads_counter += 1

        ## If the read is already marked as PCR, just ignore it
        ## It cannot be used as a reference
        if read.is_duplicate:
            pcrdup_counter += 1
            read.set_tag("XD","PCR",'Z')
            outfile.write(read)
            continue

        ## Get Barcode
        barcode = str(get_read_tag(read, args.tag))

        ## Initialize
        if barcode not in ref_barcode:
            ref_barcode[barcode] = read
        else:
            ## Compare reads to the reference (closest read)
            ref = ref_barcode[barcode]

            ## Case if read = ref
            ## Get the same tag
            if read.query_name == ref.query_name:
                ptag = get_read_tag(ref, "XD")
                if ptag is not None:
                    read.set_tag("XD",ptag,"Z")
                    read.flag += 1024
                    if ptag == "RT":
                        rtdup_counter += 1
                    elif ptag == "WIN":
                        wdup_counter += 1

            ## RT dup
            elif is_rt_dup(read, ref):
                rtdup_counter += 1
                read.flag += 1024
                read.set_tag("XD","RT",'Z')
 
            ## Window dup 
            elif is_win_dup(read, ref, args.dist):
                wdup_counter += 1
                read.flag += 1024
                read.set_tag("XD","WIN",'Z')
            else:
                ## Update reference reads only if not duplicates
                ref_barcode[barcode] = read

        ## Output
        outfile.write(read)
        
        ## Distinguish reverse and forward reads
##        if not r1.is_reverse:
            ## Set new reference
#        if barcode not in ref_perbarcode:
#            ref_perbarcode[barcode] = r1
#            outfile.write(r1)                  
#        else:
#            ## Compare with existing reference
#            ref = ref_perbarcode[barcode]
#            ref_start = get_read_start(ref)
#            r1_start = get_read_start(r1)

#            ## Is duplicates
#            if ref.tid == r1.tid and r1_start < (ref_start + dist) :
#                dup_counter += 1
#            else:
            ## update reference
#                ref_perbarcode[barcode] = r1
#                outfile.write(r1)
    
        if (reads_counter % 1000000 == 0 and args.verbose):
            print ("##", reads_counter)

samfile.close()

## Log
if args.verbose:
    dup_counter = pcrdup_counter + wdup_counter + rtdup_counter
    print ("## Number of reads: " + str(reads_counter))
    print ("## Number of pre-marked duplicates: " + str(pcrdup_counter))
    print ("## Number of window duplicates: " + str(wdup_counter))
    print ("## Number of RT duplicates: " + str(rtdup_counter))
    print ("## Number of total duplicates: " + str(dup_counter))
