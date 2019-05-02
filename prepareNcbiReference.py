import os
import re
import sys
from Bio import SeqIO
from Bio import Entrez
import subprocess
from shutil import copyfile
import argparse

# requires BioPython
# install with ~/anaconda/bin/pip install biopython

# accept command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("accession", help="NCBI Genbank accession number of reference to process")
args = parser.parse_args()

def getGenbank(appRoot, accession):

    # download temporary Genbank files for accession
    tmp_gbk_filename = 'tmp.gbk'
    tmp_gbk_filename_cleaned = 'tmp.rename.gbk'

    Entrez.email = 'dhoconno@wisc.edu'  # Always tell NCBI who you are

    # Downloading...
    net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    out_handle = open(tmp_gbk_filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved")

    # with at least some sequences (e.g., KU501215) the NCBI sequence record contains a versioning .1 or .2suffix
    # this is treated inconsistently by some tools
    # to eliminate this issue, remove versioning suffix from temporary genbank file
    # in testing, discovered that snpEff uses chromsome name from LOCUS field of Genbank file
    # other tools use accession field for chromosome name
    # use a regular expression to change the value of the LOCUS field to match value of ACCESSION field
    # biopython expects an exact number of spaces between fields, so we need to calculate the number of spaces to add after the replaced accession number

    # there is almost certainly a more elegant way to do this with biopython

    with open(tmp_gbk_filename) as infile, open(tmp_gbk_filename_cleaned, 'w') as outfile:
        for line in infile:
            # remove version info
            line = re.sub(accession + '\.[1-9]', accession, line)

            # overwrite locus field with accession number

            # number of spaces after accession
            space_ct = 23 - len(accession)
            spacer = ' ' * space_ct

            line = re.sub('LOCUS(\s*)\S*\s*', 'LOCUS' + r'\1' + accession + spacer, line)
            outfile.write(line)

    # read Genbank file and extract sequence ID

    record = SeqIO.read(tmp_gbk_filename_cleaned, "genbank")
    seq_id = record.id
    print(seq_id)

    # create reference folder for seq_id if it does not already exist
    if not os.path.exists(appRoot + '/ref/' + seq_id):
        os.makedirs(appRoot + '/ref/' + seq_id)

    # create Genbank file in seq_id location if it does not already exist
    genbank_filename = appRoot + '/ref/' + seq_id + '/' + seq_id + '.gbk'

    if not os.path.isfile(genbank_filename):

        # create Genbank file
        SeqIO.convert(tmp_gbk_filename_cleaned, "genbank", genbank_filename, "genbank")

    # create FASTA file if it doesn't exist
    fasta_filename = appRoot + '/ref/' + seq_id + '/' + seq_id + '.fa'
    if not os.path.isfile(fasta_filename):

        # create FASTA file
        SeqIO.convert(tmp_gbk_filename_cleaned, "genbank", fasta_filename, "fasta")

    # remove temporary Genbank files
    os.remove(tmp_gbk_filename)
    os.remove(tmp_gbk_filename_cleaned)

    return (genbank_filename,fasta_filename)

def indexFasta(appRoot, fasta_filename):

    # if FASTA index doesn't exist, run samtools faidx

    faidx_filename = fasta_filename + '.fai'
    if not os.path.isfile(faidx_filename):
        subprocess.call([appRoot + '/bin/samtools',
                        'faidx',
                        fasta_filename])

def createNovoindex(appRoot, fasta_filename):

    # if Novoalign index doesn't exist, create it
    novoindex_filename = fasta_filename + '.nix'
    if not os.path.isfile(novoindex_filename):
        subprocess.call([appRoot + '/bin/novoindex',
                        novoindex_filename,
                        fasta_filename])

    return novoindex_filename

def createSnpeffDb(appRoot, genbank_filename):
    # add genome to config file and create required database

    # get sequence id and description from Genbank file
    record = SeqIO.read(genbank_filename, "genbank")
    seq_id = record.id
    seq_description = record.description

    # check whether sequence id is already in snpEff config file
    if seq_id in open(appRoot + '/bin/snpEff.config').read():
        add_sequence = False
    else:
        add_sequence = True

    # append genome definition to snpEff config if sequence is not in file
    if add_sequence is True:
        with open(appRoot + '/bin/snpEff.config', "a") as myfile:
            myfile.write('\n' + seq_id + '.genome : ' + seq_description)

    # if snpEff database doesn't exist, create it
    snpeff_filename = appRoot + '/ref/' + seq_id +'/snpEffectPredictor.bin'
    if not os.path.isfile(snpeff_filename):
        # create temporary Genbank file named genes.gbk as needed by snpeff
        tmp_genbank = appRoot + '/ref/' + seq_id +'/genes.gbk'
        copyfile(genbank_filename, tmp_genbank)

        # create snpeff database
        subprocess.call(['java',
                         '-jar',
                         appRoot + '/bin/snpEff.jar',
                         'build',
                         '-genbank',
                         seq_id])

        # remove temporary genes.gbk file
        os.remove(tmp_genbank)

def makeRefIndices (appRoot, accession):
    ncbi_fetch = getGenbank(appRoot, accession)
    indexFasta(appRoot, ncbi_fetch[1])
    createNovoindex(appRoot, ncbi_fetch[1])
    createSnpeffDb(appRoot, ncbi_fetch[0])

# prepare indices
makeRefIndices(os.path.dirname(sys.argv[0]), args.accession)