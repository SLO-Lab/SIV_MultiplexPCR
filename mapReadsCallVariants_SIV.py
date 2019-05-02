import argparse
import sys
import os
import subprocess
import shutil
import fileinput

# accept command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("ref_accession", help="NCBI accession of reference sequence (note - needs to be pre-indexed with prepareGeneiousReference.py")
parser.add_argument("--fastq", nargs='+', help="path to FASTQ sequence of reads, or path to tuple of two FASTQ seuquences (XX_R1.fq.gz, XX_R2.fq.gz)")
parser.add_argument("--number_reads", help="(optional) number of reads from input FASTQ to use for mapping (default=all)")
parser.add_argument("--min_var_percent", help="(optional) exclude variants below min_var_percent from VCF reporting (e.g., 0.05 for 5%) (default=0.05)")
parser.add_argument("--qtrim", help="(optional) apply bbmap quality trimming of ends with default parameters (t/f) (default=true)")
parser.add_argument("--min_read_length", help="(optional) exclude reads from mapping shorter than this value (default=100bp)")
parser.add_argument("--debug", help="(optional) save all intermediate files for debugging (t/f) (default=f)")
args = parser.parse_args()

# if single FASTQ file, process as string
if len(args.fastq) == 1:
    fastq = args.fastq
elif len(args.fastq) == 2:
    fastq = tuple(args.fastq)

# set default values if optional arguments not provided
if args.number_reads is None:
    number_reads = '0'
else:
    number_reads = args.number_reads

if args.min_var_percent is None:
    min_var_percent = '0.01'
else:
    min_var_percent = args.min_var_percent

if args.qtrim is None:
    qtrim = 't'
else:
    qtrim = args.qtrim

if args.min_read_length is None:
    min_read_length = '100'
else:
    min_read_length = args.min_read_length

if args.debug is None:
    debug = False
elif args.debug == 'f':
    debug = False
else:
    debug = True

def downsampleFastq(appRoot, fastq, output_dir, number_reads='0', qtrim='t', minlength='100'):
    # use bbmap reformat.sh to remove low quality sequences and prune short sequences
    # can also downsample FASTQ to appropriate number of reads
    # will also output intereleaved, decompressed FASTQ regardless of input FASTQ
    # this streamlines subsequent steps by ensuring input is uncompressed FASTQ

    print('--Preprocessing FASTQ reads with bbmap reformat.sh--')

    # if single FASTQ is provided
    reformat_cmd = ''

    if len(fastq) == 1:
        reformat_cmd = [appRoot + '/bin/bbmap_36.86/reformat.sh',
                    'in=' + fastq[0],
                    'out=' + output_dir + '/reformat.fq',
                    'qtrim=' + qtrim,
                    'minlength=' + minlength,
                    'samplereadstarget=' + number_reads,
                    'sampleseed=3']

        with open(output_dir + "/log/reformat.log", "wb") as err:
            subprocess.call(reformat_cmd, stderr=err)

    # if tuple of two FASTQ is provided

    if len(fastq) == 2 :
        reformat_cmd = [appRoot + '/bin/bbmap_36.86/reformat.sh',
                    'in=' + fastq[0],
                    'in2=' + fastq[1],
                    'out=' + output_dir + '/reformat.fq',
                    'qtrim=' + qtrim,
                    'minlength=' + minlength,
                    'samplereadstarget=' + number_reads,
                    'sampleseed=3']

        with open(output_dir + "/log/reformat.log", "wb") as err:
            subprocess.call(reformat_cmd, stderr=err)

    return (reformat_cmd, output_dir + '/reformat.fq')

def mapNovoalign(appRoot, ref, fastq, output_dir):
    # map reads to novoindexed reference genome in Novoalign
    # running downsampling function converts files to interleaved FASTQ

    print('--Map reads to reference with Novoalign--')

    # set reference index to use
    novoindex = appRoot + '/ref/' + ref + '/' + ref + '.fa.nix'

    # map reads with novoalign
    # on Nov 8, 2018, SLO added the line '-r', 'ALL' so that reads that map to multiple places will actually map to both places, which is imp for LTRs
    with open(output_dir + "/mapping.sam", "wb") as out, open(output_dir + "/log/mapping.log", "wb") as err:
        mapping_cmd = [appRoot + '/bin/novoalign',
                '-d',
                novoindex,
                '-f',
                fastq,
                '-i',
                '300,50',
                '-r',
                'ALL',
                '-o',
                'SAM']

        subprocess.call(mapping_cmd, stdout=out, stderr=err)

    return mapping_cmd

def extractMappedReads(appRoot, sam, output_dir):
    """
    Extract only mapped reads from Novoalign SAM file
    Minimizes size of SAM/BAM file for subsequent steps
    :param appRoot:
    :param sam:
    :param output_dir:
    :return:
    """

    filter_sam_cmd = [appRoot + '/bin/samtools',
                'view',
                '-h',
                '-F',
                '4',
                sam]

    print('--Retain only mapped reads in SAM file--')
    with open(output_dir + "/mapping.filtered.sam", "wb") as out:
        subprocess.call(filter_sam_cmd, stdout=out)

    return filter_sam_cmd

def callVariantsVarscan(appRoot, ref, sam, output_dir, sample_name, min_percentage='0.01'):
    # sort SAM file, generate mpileup, and call variants with Varscan
    samtools_ref = appRoot + '/ref/' + ref + '/' + ref + '.fa'

    # sort SAM file
    print('--Sort SAM file and convert to BAM file--')

    with open(output_dir + "/" + sample_name + ".mapping.sorted.bam", "wb") as out:
        subprocess.call([appRoot + '/bin/samtools',
                         'sort',
                         sam,
                         '--reference',
                         samtools_ref], stdout=out)

    # generate mpileup
    print('--Generate mpileup file needed by Varscan--')
    with open(output_dir + "/mapping.sorted.pileup", "wb") as out:
        subprocess.call([appRoot + '/bin/samtools',
                         'mpileup',
                         '-f',
                         samtools_ref,
                         output_dir + "/" + sample_name + ".mapping.sorted.bam"], stdout=out)

    # run VarScan
    # report both indels and SNPs supported by more than min_percentage reads
    print('--Call variants with VarScan--')

    with open(output_dir + "/mapping.vcf", "wb") as out:
        varscan_cmd = ['java',
                         '-jar',
                         appRoot + '/bin/VarScan.v2.4.3.jar',
                         'mpileup2cns',
                         output_dir + '/mapping.sorted.pileup',
                         '--variants',
                         '--strand-filter',
                         '1',
                         '--min-var-freq',
                         min_percentage,
                         '--p-value',
                         '99e-02',
                         '--output-vcf']

        subprocess.call(varscan_cmd, stdout=out)

    return varscan_cmd

def annotateSnpeff(appRoot, ref, vcf, output_dir, sample_name):
    # run snpEff to annotate vcf file
    # only annotate variants within features (by setting ud = 0)
    print('--Annotate variant impact on coding sequence with snpEff--')

    annotated_vcf = output_dir + "/" + sample_name + ".mapping.annotated.vcf"

    with open(annotated_vcf, "wb") as out:
        snpeff_cmd = ['java',
                      '-jar',
                      appRoot + '/bin/snpEff.jar',
                      '-ud',
                      '-onlyProtein'
                      '0',
                      ref,
                      vcf]

        subprocess.call(snpeff_cmd, stdout=out)

    # replace generic Sample1 sample identifier with sample_name
    with fileinput.FileInput(annotated_vcf, inplace=True) as file:
        for line in file:
            print(line.replace('Sample1', sample_name), end='')

    # copy snpEff summary file to output location
    shutil.copyfile('snpEff_summary.html', output_dir + '/log/snpEff_summary.html')

    # remove snpEff summary file and text summary invokation location
    os.remove('snpEff_summary.html')
    os.remove('snpEff_genes.txt')

    return snpeff_cmd

def debugCleanup(output_dir):
    # if debug mode is disabled, delete these intermediate files from output folder
    os.remove(output_dir + '/mapping.sam')
    os.remove(output_dir + '/mapping.sorted.pileup')
    os.remove(output_dir + '/mapping.vcf')
    os.remove(output_dir + '/reformat.fq')
    os.remove(output_dir + '/mapping.filtered.sam')
    shutil.rmtree(output_dir + '/tmp', ignore_errors=False, onerror=None)

def copyRef(ref, appRoot, output_dir):
    ref_gbk = appRoot + '/ref/' + ref + '/' + ref + '.gbk'
    ref_gff = appRoot + '/ref/' + ref + '/' + ref + '.gff'
    ref_fasta = appRoot + '/ref/' + ref + '/' + ref + '.fa'

    if os.path.isfile(ref_gbk):
        # copy reference file in Genbank format to preserve annotations to output folder so BAM/VCF can be visualized in Geneious
        # for reference files imported from NCBI
        shutil.copyfile(ref_gbk, output_dir + '/' + ref + '.gbk')
    elif os.path.isfile(ref_gff):
        # copy reference file in GFF format and FASTA file to preserve annotations to output folder so BAM/VCF can be visualized in Geneious
        # for reference files imported from Geneious
        shutil.copyfile(ref_gff, output_dir + '/' + ref + '.gff')
        shutil.copyfile(ref_fasta, output_dir + '/' + ref + '.fa')

def mapToAnnotatedVCF(appRoot, ref, fastq, min_variant_percentage, number_reads_to_keep, qtrim, min_read_length, debug):

    # put output files in same location as input FASTQ
    # handle both single and multiple FASTQ file parameters
    if isinstance(fastq, tuple):
        assert os.path.exists(fastq[0]), 'First FASTQ file provided does not exist. Please check file path and try again.'
        assert os.path.exists(fastq[1]), 'Second FASTQ file provided does not exist. Please check file path and try again.'
        OUT_DIR = os.path.dirname(fastq[0])
        SAMPLE_NAME = os.path.basename(fastq[0]).split('.')[0]
    else:
        assert os.path.exists(fastq[0]), 'FASTQ file does not exist. Please check file path and try again.'
        OUT_DIR = os.path.dirname(fastq[0])
        SAMPLE_NAME = os.path.basename(fastq[0]).split('.')[0]

    # create temporary file directories if they do not already exist
    TMP_DIR = OUT_DIR + '/tmp'
    if not os.path.exists(OUT_DIR + '/tmp'):
        os.makedirs(OUT_DIR + '/tmp')

    # create temporary file directories if they do not already exist
    LOG_DIR = OUT_DIR + '/log'
    if not os.path.exists(OUT_DIR + '/log'):
        os.makedirs(OUT_DIR + '/log')

    # make sure reference files exist
    assert os.path.exists(appRoot + '/ref/' + ref + '/' + ref + '.fa.nix'), 'Novoalign index for reference sequence does not exist. Please check reference and try again.'

    preprocessing = downsampleFastq(appRoot, fastq, OUT_DIR, number_reads_to_keep, qtrim, min_read_length)
    mapping = mapNovoalign(appRoot, ref, preprocessing[1], OUT_DIR)
    filter = extractMappedReads(appRoot, OUT_DIR + '/mapping.sam', OUT_DIR)
    variants = callVariantsVarscan(appRoot, ref, OUT_DIR + '/mapping.filtered.sam', OUT_DIR, SAMPLE_NAME, min_variant_percentage)
    annotation = annotateSnpeff(appRoot, ref, OUT_DIR + '/mapping.vcf', OUT_DIR, SAMPLE_NAME)
    copyRef(ref, appRoot, OUT_DIR)

    # save parameters to file
    with open(OUT_DIR + '/log/' + SAMPLE_NAME + 'parameters.txt', "w") as text_file:
        text_file.write('Reference genome = ' + ref + '\n\n')
        # handle interleaved FASTQ and separate FASTQ files
        if isinstance(fastq, str):
            text_file.write('FASTQ input = ' + fastq[0] + '\n\n')
        elif isinstance(fastq, tuple):
            text_file.write('FASTQ_R1 input = ' + fastq[0] + '\n')
            text_file.write('FASTQ_R2 input = ' + fastq[1] + '\n\n')
        text_file.write('Minimum variant percentage to report = ' + min_variant_percentage + '\n\n')
        text_file.write('Debug mode = ' + str(debug) + '\n\n')
        text_file.write('Read preprocessing = ' + ' '.join(preprocessing[0]) + '\n\n')
        text_file.write('Novoalign mapping command = ' + ' '.join(mapping) + '\n\n')
        text_file.write('Remove unmapped reads = ' + ' '.join(filter) + '\n\n')
        text_file.write('Varscan variant calling command = ' + ' '.join(variants) + '\n\n')
        text_file.write('SnpEff annotation command = ' + ' '.join(annotation) + '\n\n')

    # if debug is True, keep all output files, but typically save only sorted BAM file and annotated VCF
    if debug == False:
        debugCleanup(OUT_DIR)

# execute workflow
mapToAnnotatedVCF(appRoot=os.path.dirname(sys.argv[0]),
                  ref=args.ref_accession,
                  fastq= fastq,
                  min_variant_percentage=min_var_percent,
                  number_reads_to_keep=number_reads,
                  min_read_length=min_read_length,
                  qtrim=qtrim,
                  debug=debug)