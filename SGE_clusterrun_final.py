import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re
import argparse

#########################
###### FROM NIKKI #######
#########################

# Define arguments
def get_args():
    #What this script does
    parser = argparse.ArgumentParser(description="Generate SGE cluster runs for RepeatMasker; built in RepeatMasker parameters are -xsmall [softmasks repetitive regions] -a [.align output file] -gff [generates a GFF format output] -pa [runs in parallel], please see RepeatMasker for details of these run options", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')
    #Give input genome FASTA
    parser.add_argument('-i', '--input', type=str, help='genome file in FASTA format', required=True)
    #Argument of species name
    parser.add_argument('-sp', '--species', type=str, help='Source species of query DNA FASTA', required=False)
    # Desired batch number
    parser.add_argument('-b', '--batch_count', type=int, help='Batch count', default=50)
    # Input genome directory
    parser.add_argument('-dir', '--genome_dir', type=str, help='Path to genome FASTA', required=True)
    # Argument for output directory
    parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output subdirectory', default='.')
    # Which queue to use
    parser.add_argument('-q', '--queue', type=str, help='Select the queue to run RepeatMasker in [quanah|hrothgar] with the quanah option being the general quanah omni queue, and hrothgar being the communitycluster Chewie queue', choices=['quanah', 'hrothgar'], default='quanah')
    #Argument of RepeatMasker run parameter
    parser.add_argument('-lib', type=str, help='RepeatMasker run parameter custom library "-lib [filename]" option', required=False)
    #Argument of RepeatMasker run parameter
    parser.add_argument('-xsmall', action='store_true', help='Select a RepeatMasker masking option as lowercase bases [-xsmall], default is to mask as Ns')
    #Argument of RepeatMasker run parameter
    parser.add_argument('-nolow', action='store_true', help='RepeatMasker parameter flag "-nolow" option; does not mask low complexity DNA or simple repeats')
    #Argument of RepeatMasker run parameter
    parser.add_argument('-s', '--speed', type=str, help='RepeatMasker run parameter "-q" or "-s" option; q=quick search; 5-10% less sensitive, 3-4 times faster than default; s=slow search; 0-5% more sensitive, 2.5 times slower than default', choices=['q', 's'], required=False)
    # Desired batch number
    parser.add_argument('-t', '--threads', type=int, help='Number of processors to be used per batch-***.sh file', default=12)


    args = parser.parse_args()
    GENOME = args.input
    SPECIES = args.species
    BATCH_COUNT = args.batch_count
    GENOME_DIR = args.genome_dir
    OUTDIR = args.outdir
    QUEUE = args.queue
    LIBRARY = args.lib
    XSMALL = args.xsmall
    NOLOW = args.nolow
    SPEED = args.speed
    PROC = args.threads

    return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, NOLOW, SPEED, PROC
    
GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, NOLOW, SPEED, PROC = get_args()

# Sanity checks
print("The query genome is {}.\n".format(GENOME))
print("{} batches will be made.\n".format(str(BATCH_COUNT)))
print("The genome FASTA is located in '{}'.\n".format(GENOME_DIR))
print("The output directory is '{}'.\n".format(OUTDIR))
print("The job queue is {}.\n".format(QUEUE))

if not SPECIES:
    if not LIBRARY:
        sys.exit("Must supply value for option 'species' or 'lib'!")

#FLAGS = [LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV]
FLAGS = [LIBRARY, XSMALL, NOLOW, SPEED]
if not FLAGS:
    print("Default RepeatMasker parameters used, no custom library, -inv -a -gff -pa options used.")
else:
    print("Custom parameters used:\n")
    if XSMALL:
        print("-xsmall flag used.\n")
    if NOLOW:
        print("-nolow flag used.\n")
    if LIBRARY:
        print("-lib flag used. Custom library is '{}'.\n".format(os.path.basename(LIBRARY)))
    if SPEED:
        print("-{} flag used. Search sensitivity has changed.\n".format(SPEED))

if not os.path.isdir(GENOME_DIR):
    sys.exit("The given genome directory, '{}', does not exist.".format(GENOME_DIR))

GENOME_FASTA = os.path.join(GENOME_DIR, GENOME)

if LIBRARY:
    try:
        if not os.path.getsize(LIBRARY) > 0:
            sys.exit("The library file, '{}', is empty.".format(LIBRARY))
    except OSError as e:
        sys.exit("The library file '{}' does not exist or is inaccessible.".format(LIBRARY))
'''  
try:
    if not os.path.getsize(LIBRARY) > 0:
        sys.exit("The library file, '{}', is empty.".format(LIBRARY))
except OSError as e:
    sys.exit("The library file '{}' does not exist or is inaccessible.".format(LIBRARY))
'''    
if not os.path.isdir(OUTDIR):
    sys.exit("The output directory '{}' does not exist.".format(OUTDIR))

PARTITION_DIR = os.path.join(GENOME_DIR, "RMPart")

SLOTS_PER_BATCH = 10
NUM_BATCHES = BATCH_COUNT

#check_empty(PARTITION_DIR)

if LIBRARY:
    copyfile(LIBRARY, PARTITION_DIR)
    LIB_FILE = os.path.basename(LIBRARY)
    LIBRARY = os.path.join(PARTITION_DIR, LIB_FILE)


PARTITION_DIR = os.path.abspath(PARTITION_DIR)
if not os.path.exists(PARTITION_DIR):
    try:
        os.makedirs(PARTITION_DIR)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    print("Made '{}' directory.".format(PARTITION_DIR))
else:
    if not os.listdir(PARTITION_DIR):
        print("'{}' is empty. Continuing.".format(PARTITION_DIR))
    else:
        print("'{}' is not empty. Removing contents and continuing.".format(PARTITION_DIR))
    #   os.remove(os.path.join(PARTITION_DIR, '*'))
    ## To remove anything in the RMPart folder from a previous run (all symbolic links (not expected here) and files and subdirectories) without deleting the RMPart directory itself
        for FILE in os.listdir(PARTITION_DIR):
            FILE_PATH = os.path.join(PARTITION_DIR, FILE)
            try:
                shutil.rmtree(FILE_PATH)
            except OSError:
                os.remove(FILE_PATH)


#########################
###### FROM AUSTIN ######
#########################

# Set chunk number
batch_number = BATCH_COUNT

# Create a list of all the seq records inside of the genome
records = list(SeqIO.parse(GENOME, "fasta"))

# Define the chunk size
# Sum of all lengths of all records diveded by "batch_number"
chunk_size = sum(len(i) for i in records) // batch_number + 1


# Start creating batches
def create_batch(records, chunk_size):

    # Create an object out of each record and go through them iteratively
    record_it = iter(records)

    # The "next" object in the list of "records"...
    # Basically going through each contig one at a time
    record = next(record_it)

    # Initiallize base pair counting
    current_base = 0

    # Create a dictionary for batches and initialize the batch size
    # "Batch" is defined as the output fasta file that has a collection of "chunks"
    batch = []
    batch_size = 0

    # While there are still records left in the list, keep creating new batches
    while record:

        # Loop over records untill the batch is full (i.e. reached the max chunk size), or there are no new records
        while batch_size != chunk_size and record:

            # Define the end... which sums up to the chunk size
            end = current_base + chunk_size - batch_size

            # Define the output sequence, which is the current base (beginning base), seperated by a ":" and the end base of the contig
            seq = record[current_base:end]

            # Define where to cut the contig off, which is the current base + the length of the output sequence defined above
            end_of_slice = current_base + len(seq)

            # Create the fasta headers to match that of the original SGE script
            # <original_contig_name> ":" <beginning_base> "-" <end_base>
            fasta_header = record.id + ":{}-{}".format(current_base, end_of_slice)

            # Change the seq.id to the fasta header defined above.
            seq.id = seq.name = fasta_header

            # Set a blank description for the sequence.
            # For some reason this throws off Biopython if there is nothing present in the description object.
            seq.description = ''

            # Add the sequence to the current batch
            batch.append(seq)

            # This is where we start doing the math.
            # Add the lenth of the current sequence we are iterating through to the current base.
            # When doing this, we also need to keep track of the batch_size... we want to make everything as equal as possible.
            current_base += len(seq)
            batch_size += len(seq)

            # When we have "added" all of the bases from the current sequence we are iterating through,
            # then we need to go and grab the next sequence in the list.
            if current_base >= len(record):
                record = next(record_it, None)
                current_base = 0

        # Once we have a batch with the correct size, yield the batch.
        # OR... we have run out of sequences in the genome, so stop.
        yield batch
        batch = []
        batch_size = 0

# Make the RMPart directory
#os.mkdir(PARTITION_DIR)

LIB_OR_SPECIES = ""
if LIBRARY:
    LIB_OR_SPECIES = " -lib " + LIBRARY + " "
else:
    LIB_OR_SPECIES = " -species " + SPECIES + " "

ADD_PARAMS = str(LIB_OR_SPECIES)
if NOLOW:
    ADD_PARAMS = ADD_PARAMS + "-nolow "
if XSMALL:
    ADD_PARAMS = ADD_PARAMS + "-xsmall "
if SPEED:
    ADD_PARAMS = ADD_PARAMS + "-" + str(SPEED) + " "
else: ADD_PARAMS = ADD_PARAMS + "-s "

#Set up GENOME prefix
PREFIX, EXTENSION = os.path.splitext(GENOME)

# Write out the batches as new fasta files.
for i, batch in enumerate(create_batch(records, chunk_size)):

    #Name the file and keep track of the numbering.
    filename = "{:03d}.fa".format(i)

    #Create the filepath to put the fasta file in the appropriate RMPart subdirectory
    num_dir = "{:03d}".format(i)
    path2_batch_dir = os.path.join(PARTITION_DIR, num_dir)
    os.mkdir(path2_batch_dir)
    filepath = os.path.join(PARTITION_DIR, num_dir, filename)

    # Write all the batch'e's sequences and their appropriate headers to the output fasta file.
    SeqIO.write(batch, filepath, "fasta")

    #Write the .lft file as well using the .fa file just created.
    lift_file = "{:03d}.lft".format(i)
    lift_filepath = os.path.join(PARTITION_DIR, num_dir, lift_file)

    OUTPUT = open(lift_filepath,'w')
    with open(filepath, 'r') as FIN:
        for LINE in FIN:
            if LINE.startswith('>'):
                HEADER = re.sub('>','',LINE)
                HEADER2 = re.sub('\n','',HEADER)
                PART1_HEADER = HEADER2.split(":")
                CONTIG = str(PART1_HEADER[0])
                PART2_HEADER = PART1_HEADER[1]
                SPLIT_PART2 = PART2_HEADER.split("-")
                START = int(SPLIT_PART2[0])
                END = int(SPLIT_PART2[1])
                LENGTH = END-START
                OUTPUT.write(str(START) + '\t' + str(HEADER2) + '\t' + str(LENGTH) + '\t' + str(CONTIG) + '\t' + str(END) + '\n')
    OUTPUT.close()

    #write the batch-***.sh file as well using the same system
    SH_file = "batch-{:03d}.sh".format(i)
    SGEBATCH_PATH = os.path.join(PARTITION_DIR, num_dir, SH_file)
    
    if QUEUE == 'hrothgar':
        PROJECT = 'communitycluster'
        KUE = 'Chewie'
    elif QUEUE == 'quanah':
        PROJECT = 'quanah'
        KUE = 'omni'
    else:
        sys.exit('Please choose hrothgar or quanah as the queue to use.')

#    BATCH_SH_FILE = open(SGEBATCH_PATH,'w')
    with open(SGEBATCH_PATH, 'w') as BATCH_FILE:
    #with open("RMPart/" + BATCH + "/" + QSUBFILENAME, 'w') as THISFILE:
        BATCH_FILE.write("#!/bin/sh\n")
        BATCH_FILE.write('#$ -V\n')
        BATCH_FILE.write('#$ -cwd\n')
        BATCH_FILE.write('#$ -S /bin/bash\n')
        BATCH_FILE.write("#\$ -N {}-batch-{}\n".format(PREFIX, num_dir))
        BATCH_FILE.write("#\$ -o \$JOB_NAME.o.\$JOB_ID\n")
        BATCH_FILE.write("#\$ -e \$JOB_NAME.e.\$JOB_ID\n")
        BATCH_FILE.write('#$ -q ' + KUE + '\n')
        BATCH_FILE.write('#$ -pe sm {}\n'.format(str(PROC)))
        BATCH_FILE.write('#$ -P ' + PROJECT + '\n')
        BATCH_FILE.write('\n')
        BATCH_FILE.write("cd {}\n\n".format(SGEBATCH_PATH))
        BATCH_FILE.write("{}/RepeatMasker{}-inv -a -gff -pa {} {}.fa >& run.log\n".format(REPEATMASKER, ADD_PARAMS, str(PROC - 1), num_dir))
    BATCH_FILE.close()


'''
#########################
###### FROM Dr.RAY ######
#########################



LIB_OR_SPECIES = ""
if LIBRARY:
    LIB_OR_SPECIES = " -lib " + LIBRARY + " "
else:
    LIB_OR_SPECIES = " -species " + SPECIES + " "

ADD_PARAMS = str(LIB_OR_SPECIES)
if NOLOW:
    ADD_PARAMS = ADD_PARAMS + "-nolow "
if XSMALL:
    ADD_PARAMS = ADD_PARAMS + "-xsmall "
if SPEED:
    ADD_PARAMS = ADD_PARAMS + "-" + str(SPEED) + " "
else: ADD_PARAMS = ADD_PARAMS + "-s "
        
# Build SGE file
BATCH_NUMBER = BATCH_COUNT
for BATCH in range(BATCH_COUNT):
    BATCHNUMBER = '{:03}'.format(BATCH)
    SGEBATCH_NAME = "batch-" + str(BATCH_NUMBER - 1) + ".sh"
    SGEBATCH_PATH = os.path.join(PARTITION_DIR, BATCHNUMBER)
    SGEBATCH = os.path.join(SGEBATCH_PATH, SGEBATCH_NAME)

    if QUEUE == 'hrothgar':
        PROJECT = 'communitycluster'
        KUE = 'Chewie'
    elif QUEUE == 'quanah':
        PROJECT = 'quanah'
        KUE = 'omni'
    else:
        sys.exit('Please choose hrothgar or quanah as the queue to use.')

    with open(SGEBATCH) as BATCH_FILE:
    #with open("RMPart/" + BATCH + "/" + QSUBFILENAME, 'w') as THISFILE:
        BATCH_FILE.write("#!/bin/sh\n")
        BATCH_FILE.write('#$ -V\n')
        BATCH_FILE.write('#$ -cwd\n')
        BATCH_FILE.write('#$ -S /bin/bash\n')
        BATCH_FILE.write("#\$ -N {}-batch-{}\n".format(PREFIX, BATCH_NUMBER))
        BATCH_FILE.write("#\$ -o \$JOB_NAME.o.\$JOB_ID\n")
        BATCH_FILE.write("#\$ -e \$JOB_NAME.e.\$JOB_ID\n")
        BATCH_FILE.write('#$ -q ' + KUE + '\n')
        BATCH_FILE.write('#$ -pe sm {}\n'.format(str(PROC)))
        BATCH_FILE.write('#$ -P ' + PROJECT + '\n')
        BATCH_FILE.write('\n')
        BATCH_FILE.write("cd {}\n\n".format(SGEBATCH_PATH))
        BATCH_FILE.write("{}/RepeatMasker{}-inv -a -gff -pa {} {}.fa >& run.log\n".format(REPEATMASKER, ADD_PARAMS, str(PROC - 1), BATCH_NUMBER))
    close.BATCH_FILE()

'''
#########################
###### FROM KEVIN #######
#########################


#Sanity Check
print("Creating doLift.sh file...\n")

#Create the DoLift file
DoLift_file = "doLift.sh"
DoLift_filepath = os.path.join(GENOME_DIR, DoLift_file)
OUT = open(DoLift_filepath,'w+')

#Determine which computer will be used for analyses. 
if QUEUE == 'hrothgar':
    PROJECT = 'communitycluster'
    KUE = 'Chewie'
elif QUEUE == 'quanah':
    PROJECT = 'quanah'
    KUE = 'omni'
else:
    sys.exit('Please choose hrothgar or quanah as the queue to use.')


#Begin writing to the DoLift file
OUT.write( "#! /bin/csh\n")
OUT.write( "#$ -V\n")
OUT.write( "#$ -cwd\n")
OUT.write( "#$ -N "+PREFIX+"-doLift\n")
OUT.write( "#$ -o $JOB_NAME.o.$JOB_ID\n")
OUT.write( "#$ -e $JOB_NAME.e.$JOB_ID\n")
OUT.write( "#$ -q "+KUE+"\n")
OUT.write( "#$ -pe sm 1\n")
OUT.write( "#$ -P "+PROJECT+"\n")
OUT.write("cd "+PARTITION_DIR+"/..\n")
OUT.write("\n")
OUT.write('foreach d0 ( RMPart/??? )\n')
OUT.write('  set bNum = $d0:t\n')
OUT.write('  /lustre/work/daray/software/liftUp -type=.out stdout $d0/$bNum.lft error $d0/$bNum.fa.out > $d0/$bNum.fa.liftedOut\n')
OUT.write('  if ( -e $d0/$bNum.fa.align ) then\n')
OUT.write('    /lustre/work/daray/software/liftUp -type=.align stdout $d0/$bNum.lft error $d0/$bNum.fa.align > $d0/$bNum.fa.liftedAlign\n')
OUT.write('  endif\n')
OUT.write('end\n')
OUT.write("\n")
OUT.write("/lustre/work/daray/software/liftUp "+GENOME+".out /dev/null carry RMPart/???/*.liftedOut\n")
OUT.write("/lustre/work/daray/software/liftUp "+GENOME+".align /dev/null carry RMPart/???/*.liftedAlign\n")
OUT.write("\n")
OUT.write('## Create a summary file\n')
OUT.write('# In some cases the file system delays cause the *.out file not to be available \n')
OUT.write('# Give it some time for things to settle down\n')
OUT.write('sleep 30\n')
OUT.write("/lustre/work/daray/software/RepeatMasker/util/buildSummary.pl -useAbsoluteGenomeSize -genome "+PREFIX+".2bit "+GENOME+".out > "+PREFIX+".summary\n")
OUT.write("gzip -f "+PREFIX+".summary\n")
OUT.write("\n")
OUT.write('## Generate RepeatLandscape\n')
OUT.write("perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s "+PREFIX+".divsum "+GENOME+".align\n")
OUT.write("perl /lustre/work/daray/software/RepeatMasker/util/createRepeatLandscape.pl -div "+PREFIX+".divsum -twoBit "+PREFIX+".2bit > "+PREFIX+"-landscape.html\n")
OUT.write("gzip -f "+PREFIX+".divsum\n")
OUT.write("gzip -f "+PREFIX+"-landscape.html\n")
OUT.write("\n")
OUT.write("\n")
OUT.write("##\n")
OUT.write('## NOTE NOTE NOTE: Only useful for UCSC\n')
OUT.write('## Generate data for the UCSC browser tracks\n')
OUT.write("##\n")
OUT.write("if ( -e "+GENOME+".out && -e "+GENOME+".align ) then\n")
OUT.write("  /lustre/work/daray/software/RepeatMasker/util/rmToUCSCTables.pl -out "+GENOME+".out -align "+GENOME+".align\n")
OUT.write("  gzip -f "+GENOME+".out.tsv\n")
OUT.write("  gzip -f "+GENOME+".align.tsv\n")
OUT.write("  gzip -f "+GENOME+".join.tsv\n")
OUT.write('endif \n')
OUT.write("\n")
OUT.write("##\n")
OUT.write("## NOTE NOTE NOTE: Only useful to ISB\n")
OUT.write("## Generate files for loading into RMDB\n")
OUT.write("##\n")
OUT.write("cat "+GENOME+".out | /home/rhubley/cAlign > "+GENOME+".out.c\n")
OUT.write("cat "+GENOME+".align | /home/rhubley/cAlign > "+GENOME+".align.c\n")
OUT.write("gzip -f "+GENOME+".out.c\n")
OUT.write("gzip -f "+GENOME+".align.c\n")
OUT.write("\n")
OUT.write("## Zip up final results\n")
OUT.write("gzip -f "+GENOME+".out\n")
OUT.write("gzip -f "+GENOME+".align\n")
close.OUT()
