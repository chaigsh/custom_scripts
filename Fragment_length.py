import HTSeq
import logging
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os
# import seaborn as sns
import numpy as np

"""
This program processes paired-end bam file sorted by read name,
calculates the fragment length for a pair of reads aligned to the same chromosome,
and visualizes it.
"""

## *** Logging module produces run log ***
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.DEBUG)

# Formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# FileHandler
file_handler = logging.FileHandler('result.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# StreamHandler
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

def conventional_fragment_length(bamfile=None, mapq=None):
    """
    paired-end bam file sorted by read name
    only calculates single alignment
    Fragment length = abs(r1 read alignment start position - r2 read alignment start position)
    :return:  Fragment length list
    """
    logger.info("Read into the bam file...")
    bam_reader = HTSeq.BAM_Reader(bamfile)
    total_read_count = 0
    single_alignment_read_pair_count = 0
    fragment_length_list = []
    for pair in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):
        r1_read_chr, r2_read_chr, r1_read_strand, r2_read_strand, r1_read_position, r2_read_position = \
            None, None, None, None, None, None
        if len(pair) == 1:  ## single alignment
            r1_read, r2_read = pair[0]
            if r1_read:
                if r1_read.aligned:
                    if r1_read.aQual >= mapq and r1_read.flag != 1024:
                        r1_read_chr = r1_read.iv.chrom
                        r1_read_position = r1_read.iv.start_d
                        r1_read_strand = r1_read.iv.strand
            if r2_read:
                if r2_read.aligned:
                    if r2_read.aQual >= mapq and r2_read.flag != 1024:
                        r2_read_chr = r2_read.iv.chrom
                        r2_read_position = r2_read.iv.start_d
                        r2_read_strand = r2_read.iv.strand
        if r1_read_chr:
            if r1_read_chr == r2_read_chr:
                if (r1_read_strand == '+' and r2_read_strand == '-') or \
                    (r1_read_strand == '-' and r2_read_strand == '+'):
                    single_alignment_read_pair_count += 1
                    Length = r1_read_position - r2_read_position
                    fragment_length_list.append(abs(Length))
        total_read_count += 1
        if not total_read_count % 1000000:
            logger.info("%s M reads have been processed" % str(total_read_count/float(1000000)))
    return total_read_count, single_alignment_read_pair_count, fragment_length_list

def conditional_fragment_length(bamfile=None, mapq=None):
    """
    paired-end bam file sorted by read name
    only calculates single alignment
    the 5' end for r1 read and r2 read must align to genome
    only calculates single alignment
    Fragment length = abs(r1 read alignment start position - r2 read alignment start position)
    :return: Fragment length list
    """
    logger.info("Read into the bam file...")
    bam_reader = HTSeq.BAM_Reader(bamfile)
    total_read_count = 0
    single_alignment_read_pair_count = 0
    fragment_length_list = []
    for pair in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):
        r1_read_chr, r2_read_chr, r1_read_strand, r2_read_strand, r1_read_position, r2_read_position = \
            None, None, None, None, None, None
        if len(pair) == 1:  ## single alignment
            r1_read, r2_read = pair[0]
            if r1_read:
                if r1_read.aligned:
                    if r1_read.aQual >= mapq and r1_read.flag != 1024:
                        if (r1_read.iv.strand == '+' and r1_read.cigar[0].type == 'M') or (
                                        r1_read.iv.strand == '-' and r1_read.cigar[-1].type == 'M'):
                            r1_read_chr = r1_read.iv.chrom
                            r1_read_position = r1_read.iv.start_d
                            r1_read_strand = r1_read.iv.strand
            if r2_read:
                if r2_read.aligned:
                    if r2_read.aQual >= mapq and r2_read.flag != 1024:
                        if (r2_read.iv.strand == '+' and r2_read.cigar[0].type == 'M') or (
                                        r2_read.iv.strand == '-' and r2_read.cigar[-1].type == 'M'):
                            r2_read_chr = r2_read.iv.chrom
                            r2_read_position = r2_read.iv.start_d
                            r2_read_strand = r2_read.iv.strand
        if r1_read_chr:
            if r1_read_chr == r2_read_chr:
                if (r1_read_strand == '+' and r2_read_strand == '-') or \
                        (r1_read_strand == '-' and r2_read_strand == '+'):
                    single_alignment_read_pair_count += 1
                    Length = r1_read_position - r2_read_position
                    fragment_length_list.append(abs(Length))
        total_read_count += 1
        if not total_read_count % 1000000:
            logger.info("%s M reads have been processed" % str(total_read_count / float(1000000)))
    return total_read_count, single_alignment_read_pair_count, fragment_length_list

def visualization_it(data=None, title=None):
    # x= np.random.randn(100)
    # x = list(x)
    # with PdfPages('multipage_pdf.pdf') as pdf:
    plt.hist(data, color='blue', edgecolor='blue')
    plt.title(title, fontsize=10)
    plt.xlabel('Fragment Length (bp)', fontsize=8)
    plt.ylabel('Read Pair Count', fontsize=8)
    pdf.savefig()
    plt.close()

    plt.hist(data, density=True, color='blue', edgecolor='blue')
    plt.title(title, fontsize=10)
    plt.xlabel('Fragment Length (bp)', fontsize=8)
    plt.ylabel('Frequency', fontsize=8)
    pdf.savefig()
    plt.close()
    return

def cmd_parameters():
    parser = argparse.ArgumentParser(
                                     usage='python %(prog)s bam_files [options]',
                                     description=
                                     '''This script takes one or multiple paired-end alignment BAM files
                                     (supported formats: BAM), performs a simply fragment length calculation
                                     and produces histogram plots. The plots are output as a PDF file.''',
                                     epilog=
                                     '''Written by chaigsh (chaigsh@mail2.sysu.edu.cn), Lab of Bioinformatics
                                     and Epigenomics.(c) 2019. Released under the terms of the GNU General
                                     Public License v3.'''
                                     )
    parser.add_argument('bamfile', type=str, nargs='+',
                        help='Required positional parameter. Only BAM format is ok, \
                        and the BAM file must be sorted by read name. Multiple BAM files are separated by space.')
    parser.add_argument("-o", "--outfile", dest="outfile", default='results.pdf', type=str,
                        help="output filename (default is results.pdf)")
    parser.add_argument("--conditional", action="store_true", dest="conditional",
                        help="Only process r1 and r2 read pairs which the 5' end must align to genome"
                        )
    parser.add_argument("--mapq", type=int, dest="mapq", default=10,
                        help="the mapping quality that appears in the BAM file (default: 10)")
    parser.add_argument("--min", type=int, dest="min", default=0,
                        help="the min of x axis for the histogram plot (default: 0)")
    parser.add_argument("--max", type=int, dest="max", default=1000,
                        help="the max of x axis for the histogram plot (default: 1000)")
    args = parser.parse_args()
    print(args)
    logger.info(args)
    return args

if __name__ == "__main__":
    args = cmd_parameters()
    Fragment_length_list = None
    if args.bamfile:
        # pdf = PdfPages(args.outfile)
        with PdfPages(args.outfile) as pdf:
            for eachfile in args.bamfile:
                Title = os.path.basename(eachfile)
                logger.info(Title + " is processing!")
                if args.conditional:
                    Total_read_count, Single_alignment_read_pair_count, Fragment_length_list = \
                        conditional_fragment_length(bamfile=eachfile, mapq=args.mapq)
                else:
                    Total_read_count, Single_alignment_read_pair_count, Fragment_length_list = \
                        conventional_fragment_length(bamfile=eachfile, mapq=args.mapq)
                if Fragment_length_list:
                    logger.info("Produce plot for %s sample" %Title)
                    # visualization_it(data=Fragment_length_list, title=Title)
                    plt.hist(Fragment_length_list, color='blue', edgecolor='blue',
                             bins=np.arange(1001))
                    plt.xlim(args.min, args.max)
                    plt.title(Title, fontsize=10)
                    plt.xlabel('Fragment Length (bp)', fontsize=8)
                    plt.ylabel('Read Pair Count', fontsize=8)
                    pdf.savefig()
                    plt.close()

                    plt.hist(Fragment_length_list, density=True, color='blue', edgecolor='blue',
                             bins = np.arange(1001))
                    # plt.xlim(0, 1000)
                    plt.xlim(args.min, args.max)
                    plt.title(Title, fontsize=10)
                    plt.xlabel('Fragment Length (bp)', fontsize=8)
                    plt.ylabel('Frequency', fontsize=8)
                    pdf.savefig()
                    plt.close()
    else:
        logger.error("No bam files!!")