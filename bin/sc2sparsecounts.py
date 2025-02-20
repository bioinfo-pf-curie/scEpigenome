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
from collections import OrderedDict
from collections import defaultdict
import pysam
import numpy as np
from scipy import sparse
import scipy.io as sio
from bx.intervals.intersection import Intersecter, Interval
import gzip
import shutil
    
def timing(function, *args):
    """                              
    Run a fonction and return the run time and the result of the function
    If the function requires arguments, those can be passed in
    """
    startTime = time.time()
    result = function(*args)
    print('%s function took %0.3f ms' % (function.func_name, (time.time() - startTime) * 1000))
    return result


def load_BED(in_file, featuresOverCoord=False, verbose=False):
    """
    Read a BED file and store the intervals in a tree
    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    BED file are half-open, meaning that a bin ]100, 200] covered the bases 101 to 200
    
    in_file = input file [character]
    verbose = verbose mode [logical]
    """
    x = {}
    if verbose:
        print("BED FILE :")
        print("## Loading BED file '", in_file, "'...")
    featureNames=[]
    nline = 0
    with open(in_file) as bed_handle:
        for line in bed_handle:
            if nline > 0 and nline % 5000==0 and verbose: 
                print("## %d features loaded ..." % nline)
            nline +=1
            bedtab = line.split("\t")
            chromosome, start, end = bedtab[:3]
            if len(bedtab)>3 & featuresOverCoord==True:
                name = bedtab[4]

            # BED files are zero-based, half-open as Intervals objects
            start = int(start) 
            end = int(end)
            if featuresOverCoord==True:
                featureNames.append(name.strip())
            else:
                featureNames.append(chromosome + ":" + str(start) + "-" + str(end))
            
            if chromosome in x:
                tree = x[chromosome]
                tree.add_interval(Interval(start, end, value={'pos' : nline - 1}))
            else:
                tree = Intersecter()
                tree.add_interval(Interval(start, end, value={'pos' : nline - 1}))
                x[chromosome] = tree
    bed_handle.close()
    return (x, featureNames)


def get_barcode_number_from_rg(sam):
    """
    Read the BAM header, counts the number of @RG
    """
    barcode_number = None
    if 'RG' in sam.header:
        items = sam.header['RG']
        barcode_number = len(items)
    return barcode_number


def get_barcode_number_from_header(sam):
    """
    Read the BAM header, extract the CO tag with barcode number that 
    is added during the flag step
    """
    barcode_number = None
    if 'CO' in sam.header:
        COitems = sam.header['CO']
        for com in COitems:
            if re.search("Barcodes", com):        
                barcode_number = int(com.split(":")[1])
    return barcode_number


def get_read_tag(read, tag):
    """
    Extract a flag from a read alignment
    """
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return None

def get_read_start(read):
    """
    Return the 5' end of the read
    Same as reference_start / reference_end in recent pysam version
    """
    if read.is_reverse:
        pos = read.pos + read.alen
    else:
        pos = read.pos
    return int(pos)


def get_chromosome_size_from_header(sam):
    """
    Extract chromosome size from header. 
    Chromosome are then reorder chr1-22,X,Y
    Any other contigs are removed
    """
    chrID = [ 'chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y'] ]
    chromSize = OrderedDict()
    chromSizeOrdered = OrderedDict()
    SQitems = sam.header['SQ']
    for chrom in SQitems:
        if chrom['SN'] in chrID:
            chromSize[chrom['SN']] = chrom['LN']
    for chrom in chrID:
        if chrom in chromSize:
            chromSizeOrdered[chrom]=chromSize[chrom]
    return(chromSizeOrdered)


def get_chromosome_bins(chroms, bsize):
    """
    Get number of genomic bins per chromosome
    Is used to tranform genomic coordinates, into a bin number in 
    the count matrix
    """
    chrombins = OrderedDict()
    for chrname in chroms:
        x = float(chroms[chrname]) / float(bsize)
        chrombins[chrname] = int(math.ceil(x))
    return(chrombins)


def get_features_idx(intervals, chrom, read, useWholeRead=True, verbose=False):
    """
    Intersect a given read with the set of intervals
    intervals = the fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]
    useWholeRead = True/False, use either the 5' end of the read or the full length
    """    
    if useWholeRead:
        lpos = read.pos
        rpos = read.pos + read.alen
    else :
        if read.is_reverse:
            rpos = get_read_start(read)
            lpos = rpos - 1
        else:
            lpos = get_read_start(read)
            rpos = lpos + 1

    if chrom in intervals:
        # Overlap with the left/right positions of the read (zero-based)
        feat = intervals[chrom].find(lpos, rpos)
        if len(feat) == 0:
            if verbose: print >> sys.stderr, "Warning - no feature found for read at", chrom, ":", read.pos, "- skipped"
            return None
        else:
            feat_idx = []
            for i in range(len(feat)): 
                feat_idx.append(feat[i].value['pos'])
            return feat_idx
    else:
        if verbose: print >> sys.stderr, "Warning - no feature found for read at", chrom, ":", read.pos, "- skipped"
        return None



def get_bin_idx (read, chrbins_cumsum, binSize, useWholeRead=True):
    """
    For a given chromosome and position, return the bin indice
    Return a zero-based index
    """
    chrname = read.reference_name
    if chrname in chrbins_cumsum.keys():
        gbins = chrbins_cumsum[chrname]
    else:
        print("Chromosome " + chrname + " not found !", file=sys.stderr)
        sys.exit(-1)

    if useWholeRead:
        lpos = read.pos
        rpos = read.pos + read.alen
        if read.is_reverse:
            ## Require for reverse reads that start exactly at the end of a bin ...
            rpos = rpos - 1

        ### Sum of previous chromosome + current bin (0-based)
        idx = gbins + int(math.floor(float(lpos) / float(binSize)))
        idx_end = gbins + int(math.floor(float(rpos) / float(binSize)))
        
        if idx != idx_end:
            return range(idx, idx_end + 1)
        else:
            return [idx]
    else:
        lpos = get_read_start(read)
        if read.is_reverse:
            lpos = lpos - 1
        idx = gbins + int(math.floor(float(lpos) / float(binSize)))
        return [idx]
    

def get_bins_coordinates(i, chromsize, chrom_idx, bsize):
    """
    Transform a bin indice into a genomic coordinate
    Need indice, chromosome size, number of indices per chromosome and binsize
    """
    cumidx = 0  
    for k in chrom_idx.keys():
        if cumidx == 0:
            chrname = k
        if cumidx + chrom_idx[k] > i:
            chrname = k
            break
        else:
            cumidx += chrom_idx[k]
    
    start = i * bsize - cumidx * bsize
    end = int(start) + int(bsize)
    if end > chromsize[chrname]:
        end = chromsize[chrname]
    return np.array([str(chrname), str(start), str(end)])


def select_mat(x, barcodes, nreads=500, verbose=False):
    """
    Select counts matrix columns
    """
    cx = x.tocsr()
    if len(barcodes) < cx.shape[1]:
        print("## Only " + str(len(barcodes)) + " detected barcodes. Reducing matrix shape.")
        cx = cx[:,:len(barcodes)]

    cols = np.array(cx.sum(axis=0))  # sum the columns
    idx = np.where(cols >= nreads)[1]
    if verbose:
        print("BARCODE SUMMARY :")
        print("## Initial number of barcodes: ", len(barcodes))
        print("## Number of barcodes with at least ", nreads , " counts falling on TSS regions.: ", len(idx))
    return idx


def saveSparseMatrix(x, colnames, odir, rownames=None, chromsize=None, chromidx=None, bsize=None, verbose=False):
    """
    Write the count table into a txt file without taking too much RAM
    Note that most of the args are used to convert bin coordinate into genomic coordinates
    """
    if verbose:
        print("## Writting output file ...")

    cx = x.tocsr()
    if not os.path.exists(odir):
        os.mkdir(odir)

    mtx_file = odir + "/matrix.mtx"
    features_file = odir + "/features.tsv.gz"
    barcodes_file = odir + "/barcodes.tsv.gz"

    # Save matrix.mtx
    sio.mmwrite(mtx_file, cx, field='integer')
    with open(mtx_file,'rb') as mtx_in:
        with gzip.open(mtx_file + '.gz', 'wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(mtx_file) 
    
    # Save barcodes.tsv
    handle = gzip.open(barcodes_file,'wb')
    colnames = np.array([colnames])
    np.savetxt(handle, colnames, '%s', delimiter="\n")
    handle.close()

    # Save features.tsv
    handle = gzip.open(features_file,'wb')  
    ## Genomic bins
    if chromsize is not None and chromidx is not None and bsize is not None:
        for i in range(cx.shape[0]):
            coord = get_bins_coordinates(i, chromsize, chromidx, bsize)
            coord = np.array([coord])
            np.savetxt(handle, coord, fmt='%s', delimiter="\t")
            if (i > 0 and i % 10000 == 0 and verbose):
                print("## Write line " + str(i))
    elif rownames is not None:
        finalNbLines=0  
        for i in range(cx.shape[0]):
            name = np.array([rownames[i]])
            np.savetxt(handle, name, '%s', delimiter="\t")
            finalNbLines+=1
            if (i > 0 and i % 10000 == 0 and verbose):
                print("## Write line " + str(i)) # je ne comprends pas cette ligne 
        print("tot lines:", finalNbLines)
    handle.close()


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    read_pair_counts = defaultdict(lambda: [0, 0])
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
                read_pair_counts[read.query_name][0] += 1
            elif read.is_read2:
                read_pair_counts[read.query_name][1] += 1

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
    for name, (read1_count, read2_count) in read_pair_counts.items():
        if read2_count > 0 and read1_count == 0:
            print(f"Read2 orphelin : {name}")
            yield read_dict[name][1], None
            del read_dict[name]


if __name__ == "__main__":

    # Init variables
    pair_counter = 0
    singleton_counter = 0
    non_overlapping_counter = 0
    overlapping_counter = 0
   
    # Reads args
    parser = argparse.ArgumentParser(prog='sc2counts.py', description="Transform a BAM file to a count table based on barcode information and genomic bins/features", 
                                     epilog="Note that --bin and --bed options are exclusive ! Counts are generated either per bin (--bin) or per genomics features (--bed)")
    parser.add_argument('-i', '--input', help="BAM file with barcode tag", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-b', '--bin', help="Size of genomic bin size", default=None, type=int)
    group.add_argument('-B', '--bed', help="BED file of genomic features", default=None, type=str)
    parser.add_argument('-o', '--output', help="Output directory where to save matrix.mtx, features.tsv and barcodes.tsv files.", default="./sparse_matrix", type=str)
    parser.add_argument('-s', '--barcodes', help="Number of barcodes in the BAM file. Default: Extracted from the BAM 'CO' field", type=int)
    parser.add_argument('-t', '--tag', help="Barcode Tag. Default: XB", default="XB", type=str)
    parser.add_argument('-f', '--filt', help="Select barcodes with at least FILT counts. Default: None", default=None, type=str)
    parser.add_argument('-w', '--useWholeRead', help="Use the whole read in the count instead of the 5' end. Default: False", default=False, action="store_true")
    parser.add_argument('-F', '--featuresOverCoord', help="When counting on BED file, write feature name (column 4 of BED) as rownames of count matrix instead of coordinates. Default: False", 
                        default=False, action="store_true")
    parser.add_argument('-v', '--verbose', help="", action="store_true")
    parser.add_argument('-d', '--debug', help="", action="store_true")
 
    args = parser.parse_args()
 
    # check args
    if args.bed is None and args.bin is None:
        print("Error: --bin or --bed must be defined")
        sys.exit(-1)
        
    # Verbose mode
    if args.verbose:
        print("------")
        print("## sc2counts.py")
        print("## input =", args.input)
        print("## output =", args.output)
        print("## tag =", args.tag)
        print("## binSize =", args.bin)
        print("## bedFile =", args.bed)
        print("## barcodeNumber =", args.barcodes)
        print("## minCounts =", args.filt)
        print("## useWholeReads =", args.useWholeRead)
        print("## verbose =", args.verbose)
        print("## featuresOverCoord =", args.featuresOverCoord)
        print("-----")

    # Read the SAM/BAM file
    if args.verbose:
        print("INPUT FILE :")
        print("## Opening SAM/BAM file '", args.input,"'...")
    samfile = pysam.AlignmentFile(args.input, 'rb')

    # Get info from header
    chromsize = get_chromosome_size_from_header(samfile)
    if len(chromsize)==0:
        print("Error : chromosome lengths not available in BAM file. Exit", file=sys.stderr)
        sys.exit(-1)
    else:
        print("## " + ','.join(key for key, val in chromsize.items()))

    # Get counts dimension
    if args.barcodes is None:
        if args.verbose:
            print("BARCODE EXTRACTION :")
        if args.tag == "RG":
            N_barcodes = get_barcode_number_from_rg(samfile)
            if args.verbose:
                print("## Barcode number from RG: " + str(N_barcodes))
        else:
            N_barcodes = get_barcode_number_from_header(samfile)
            if N_barcodes is None :
                print("Erreur : unable to find barcodes number. Exit", file=sys.stderr)
                sys.exit(-1)
            if args.verbose:
                print("## Barcode number from header: " + str(N_barcodes))
    elif args.barcodes is not None:
        N_barcodes = args.barcodes

    # Chromosome/feature bins
    if args.bin is not None:
        ## Get number of bins per chromosome
        chromsize_bins = get_chromosome_bins(chromsize, args.bin)
        ## Calculate cumsum 
        csum=np.cumsum(list(chromsize_bins.values()))
        csum=np.append(0, csum[:-1])
        #csum=csum - csum[0]
        chromsize_bins_cumsum = dict(zip(chromsize_bins.keys(), csum))
        N_bins = sum(list(chromsize_bins.values()))
        if args.verbose:
            print("## Number of bins: " + str(N_bins))
    elif args.bed is not None:
        feat_bins = load_BED(args.bed, args.featuresOverCoord, args.verbose)
        N_bins = len(feat_bins[1]) 
        if args.verbose:
            print("## Number of Features from bed file: " + str(N_bins))
 
    ## Init count table
    ## Note that lil matrix is a sparse (ie. RAM eficient) matrix
    ## design to speed incrementation
    counts = sparse.lil_matrix((N_bins, N_barcodes))
    allbarcodes = {}
    start_time = time.time()
    # Vérifier la présence dans chromsize uniquement pour les reads existants
    for read1, read2 in read_pair_generator(samfile):
        # Gérer le cas où le read est un singleton
        if read2 is None:
            singleton_counter += 1 
        else: 
            pair_counter += 1
            
        if read1.reference_name not in chromsize.keys():
            continue
        if read2 is not None:
            if read2.reference_name not in chromsize.keys():
                continue
        
        print("fragment name: ", read1.query_name)

        ## get chrom name
        r1_chrom = read1.reference_name

        ## Get barcode
        barcode = str(get_read_tag(read1, args.tag))
        
        ## Get Barcode (ie col) indices
        if barcode in allbarcodes.keys():
            j = allbarcodes[barcode]
        else:
            allbarcodes[barcode]=len(allbarcodes)
            j = len(allbarcodes)-1

        ## Get bin indice (rows) and increment count matrix with one read count 
        if args.bin is not None:
            i = get_bin_idx(read1, chromsize_bins_cumsum, args.bin, useWholeRead=args.useWholeRead)
            for ii in i: 
                counts[ii, j] += 1
        ## Get features indice (rows) and increment count matrix with one read count 
        elif args.bed is not None:
            i = get_features_idx(feat_bins[0], r1_chrom, read1, useWholeRead=args.useWholeRead)
            if i is not None:
                for ii in i: 
                    counts[ii, j] += 1
                overlapping_counter += 1
            else:
                non_overlapping_counter += 1

        if (pair_counter % 1000000 == 0 and args.debug):
            stop_time = time.time()
            print( "##", pair_counter, stop_time-start_time )
            start_time = time.time()
            break
    samfile.close()

    if args.verbose and non_overlapping_counter > 0:
        tot_frag=pair_counter+singleton_counter
        print("READS :")
        print("## Number of paired reads:", pair_counter)
        print("## Number of singletons :", singleton_counter)
        print("Total fragments:", tot_frag)
        print("## Fragment overlaping features:", overlapping_counter)
        print("## Fragment NOT overlaping features:", non_overlapping_counter)

    ## Filter count matrix
    if args.filt is not None:
        filters = map(int, args.filt.split(","))
        for filt in filters:
            sel_idx = select_mat(x=counts, barcodes=allbarcodes, nreads=filt, verbose=args.verbose)
            counts_reduced = counts[:, sel_idx]
            z=np.array(list(allbarcodes.keys()))
            allbarcodes_reduced = np.array(list(allbarcodes.keys()))[sel_idx]
            allbarcodes_reduced = allbarcodes_reduced.tolist()
            ## save Matrix
            if args.bin is not None:
                saveSparseMatrix(counts_reduced, allbarcodes_reduced, args.output, chromsize=chromsize, chromidx=chromsize_bins, bsize=args.bin, verbose=args.verbose)
            elif args.bed is not None:
                saveSparseMatrix(counts_reduced, allbarcodes_reduced, args.output, rownames=feat_bins[1], verbose=args.verbose)
    else:
        ## save Matrix
        if args.bin is not None:
            saveSparseMatrix(counts, allbarcodes, args.output, chromsize=chromsize, chromidx=chromsize_bins, bsize=args.bin, verbose=args.verbose)
        elif args.bed is not None:
            saveSparseMatrix(counts, allbarcodes, args.output, rownames=feat_bins[1], verbose=args.verbose)



