import subprocess
import sys
import pyfasta
import os
import shutil
import gzip
import argparse

# accept command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("short_amplicon_fasta_ref", help="path to FASTA file containing individual amplicon sequences to normalize against")
parser.add_argument("fastq_R1", help="path to first paired-end FASTQ file")
parser.add_argument("fastq_R2", help="path to second paired-end FASTQ file")
parser.add_argument("reads_per_amplicon", help="approximate number of reads from each amplicon to save")
args = parser.parse_args()

def downsampleAmplicons(appRoot, R1_FASTQ, R2_FASTQ, READS_PER_AMPLICON, SHORT_AMPLICON_FASTA_REF):
    # for virus sequences generated with Nick Loman ZIKV technique using multiplexed pools with many small amplicons
    # normalizes coverage for each amplicon by mapping full readset to small amplicon and retaining only subset of mapped reads
    # this is essential for normalizing coverage differences in the 30 amplicons that comprise one genome
    # original version of this tool in exp 18442
    # this version eliminates mapping of normalized reads to reference because exp 18582 workflow for mapping and variant calling is better

    # get path to input FASTQ file
    # save downsampled FASTQ file to same location
    FASTQ_BASENAME = os.path.dirname(R1_FASTQ)

    # create temporary file directories if they do not already exist
    TMP_DIR = FASTQ_BASENAME + '/tmp'
    if not os.path.exists(FASTQ_BASENAME + '/tmp'):
        os.makedirs(FASTQ_BASENAME + '/tmp')

    # create temporary file directories if they do not already exist
    LOG_DIR = FASTQ_BASENAME + '/log'
    if not os.path.exists(FASTQ_BASENAME + '/log'):
        os.makedirs(FASTQ_BASENAME + '/log')

    # remove temporary directory if it already exists
    if os.path.exists(TMP_DIR):
        shutil.rmtree(TMP_DIR, ignore_errors=False, onerror=None)

    # create temporary directory and subdirectories
    os.makedirs(TMP_DIR)
    os.makedirs(TMP_DIR + '/split_reference_fasta')
    os.makedirs(TMP_DIR + '/mapped_reads')
    os.makedirs(TMP_DIR + '/filtered_reads')

    # merge overlapping R1 and R2 reads
    # trim 22bp from 5' end of each sequence
    # examine 30461 technical replicate first

    bbmerge = appRoot + '/bin/bbmap_36.86/bbmerge.sh'

    # get sample name to use for merged reads and downsmapled reads
    SAMPLE_NAME = os.path.basename(R1_FASTQ).split('.')[0]
    MERGED_FASTQ_OUT = FASTQ_BASENAME + '/' + SAMPLE_NAME + '.merged.fastq.gz'
    LEFT_TRIM = '22'

    # run bbmerge
    # capture stderr statistics to log file

    with open(LOG_DIR + '/bbmerge.log', 'w') as f:
        subprocess.call([bbmerge,
                         'in=' + R1_FASTQ,
                         'in2=' + R2_FASTQ,
                         'forcetrimleft=' + LEFT_TRIM,
                         'out=' + MERGED_FASTQ_OUT], stderr=f)

    # split reference FASTA into single sequences
    # this reference sequence has each of the amplicons separated as individual sequences by Shelby and Katie
    # note that this is different than the preferred ZIKV PR reference for mapping in a later step
    # the use of Shelby and Katie's sequence here won't hurt anything, though

    # store reference FASTA in array
    f = pyfasta.Fasta(SHORT_AMPLICON_FASTA_REF)

    # iterate over each entry in REF_FASTA and create new FASTA file with single sequence
    for i in f:
        header = str(i)
        sequence = str(f[i])

        # create FASTA file
        file = open(TMP_DIR + '/split_reference_fasta/' + header + '.fasta', "w")
        file.write('>' + header + '\n')
        file.write(sequence + '\n')
        file.close()

        # map meregd reads to each individual amplicon FASTA

    # path to folder of individual reference files
    SPLIT_REF_DIR = TMP_DIR + '/split_reference_fasta'

    # path to bbmap
    BBMAP = appRoot + '/bin/bbmap_36.86/bbmap.sh'

    # path to reformat.sh
    REFORMAT = appRoot + '/bin/bbmap_36.86/reformat.sh'

    # path to seqtk
    SEQTK = appRoot + '/bin/seqtk'

    # run BBMAP on each reference sequence individually
    # save mapped reads for each file separately

    # initialize counter
    current_ct = 0
    amplicon_ct = len(os.listdir(SPLIT_REF_DIR))

    for fn in os.listdir(SPLIT_REF_DIR):

        # print update on data processing
        current_ct +=1
        print('--Read normalizing is ' + str("%.1f" % (current_ct / amplicon_ct * 90)) + '% complete.--')

        # run bbmap
        # default sensitivty of minid=0.76 sufficient to capture only merged reads that map to a single amplicon
        # store reference in memory with nodisk

        with open(LOG_DIR + '/bbmap.log', 'w') as j:
            subprocess.call([BBMAP,
                             'in=' + MERGED_FASTQ_OUT,
                             'nodisk=t',
                             'ref=' + SPLIT_REF_DIR + '/' + fn,
                             'outm=' + TMP_DIR + '/mapped_reads/' + str(fn) + '.bam'], stderr=j)

        # run bbmap reformat.sh to filter bam file to include only reads spanning entire amplicon
        # Shelby's amplicon file includes 44bp of primer that has been trimmed during merging
        # so set "full length" to size of amplicon minus 60bp to allow some flexibility

        # get length of reference sequence
        g = pyfasta.Fasta(SPLIT_REF_DIR + '/' + fn)

        for i in g:
            REF_LENGTH = len(str(g[i]))

        MIN_LENGTH = str(REF_LENGTH - 60)

        with open(LOG_DIR + '/filter_bam.log', 'w') as k:
            subprocess.call([REFORMAT,
                             'in=' + TMP_DIR + '/mapped_reads/' + str(fn) + '.bam',
                             'out=' + TMP_DIR + '/filtered_reads/' + str(fn) + '.filtered.fastq.gz',
                             'minlength=' + MIN_LENGTH], stderr=k)

        # downsample filtered reads so each amplicon has same number of reads
        # use seqtk to filter
        # save stdout from seqtk sample to file

        with open(FASTQ_BASENAME + '/downsampled.fastq', 'a') as h:
            subprocess.call([SEQTK,
                             'sample',
                             TMP_DIR + '/filtered_reads/' + str(fn) + '.filtered.fastq.gz',
                             READS_PER_AMPLICON], stdout=h)

    # gzip compress downsampled FASTQ

    f_in = open(FASTQ_BASENAME + '/downsampled.fastq', 'rb')
    f_out = gzip.open(FASTQ_BASENAME + '/' + SAMPLE_NAME + '.' + READS_PER_AMPLICON + '.reads.per.amplicon.fastq.gz', 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    print('--Read normalizing is 100% complete.--')

    # remove uncompressed FASTQ and TMP folder
    os.remove(FASTQ_BASENAME + '/downsampled.fastq')
    shutil.rmtree(TMP_DIR, ignore_errors=False, onerror=None)

downsampleAmplicons(appRoot=os.path.dirname(sys.argv[0]),
                   R1_FASTQ=args.fastq_R1,
                   R2_FASTQ=args.fastq_R2,
                   READS_PER_AMPLICON=args.reads_per_amplicon,
                   SHORT_AMPLICON_FASTA_REF=args.short_amplicon_fasta_ref)