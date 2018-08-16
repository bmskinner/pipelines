#!/usr/bin/env python 

from collections import OrderedDict
from datetime import datetime
# from joblib import Parallel, delayed
from pipeline import Sample, CommandFactory, Pipeline
import argparse
import errno
import glob
import logging
import multiprocessing
import os
import os.path
import pickle
import Queue
import re
import signal
import subprocess
import sys
import time
import pipeline


__version__ = 'v01'

logger_name = "analysis"

class SraSample(Sample):
    '''Class to hold sample data and simplify access to 
       folders
    '''
    def __init__(self, sra, breed, sex, mbases):
        Sample.__init__(self, sra)
        self._sra_acc = sra
        self._breed = breed
        self._sex = sex
        self._mbases = mbases

    def getAccession(self):
        return self._sra_acc

    def getLogFile(self, suffix):
        return args.basedir + ("logs/%s.%s.log" % (self._sra_acc, suffix))

    def getFastqFolder(self):
        return (args.basedir+"fastq/"+self._breed+"/")

    def getFastq(self):
        return (self.getFastqFolder()+self._sra_acc+".fastq.gz")

    def getFastqRead1(self):
        return (self.getFastqFolder()+self._sra_acc+"_1.fastq.gz")

    def getFastqRead2(self):
        return (self.getFastqFolder()+self._sra_acc+"_2.fastq.gz")

    def hasFastqFiles(self):
        return self.has_file(self.getFastqRead1()) and self.has_file(self.getFastqRead2())

    def getTrimmedReadFolder(self):
        return (args.basedir+"trimmed_reads/"+self._breed+"/")

    def getTrimmedRead1(self):
        return (self.getTrimmedReadFolder()+self._sra_acc+"_1_val_1.fq.gz")

    def getTrimmedRead2(self):
        return (self.getTrimmedReadFolder()+self._sra_acc+"_2_val_2.fq.gz")

    def hasTrimmedReads(self):
        return self.has_file(self.getTrimmedRead1()) and self.has_file(self.getTrimmedRead2())

    def getMappedReadFolder(self):
        return (args.basedir+"mapped_reads/"+self._breed+"/")

    def getMappedReadFile(self):
        return (self.getMappedReadFolder()+self._sra_acc+".sam")

    def hasMappedReadFile(self):
        return self.has_file(self.getMappedReadFile())

    def hasMappedReads(self):
        return self.has_file(self.getMappedReadFile())

    def getBamFile(self):
        return (self.getMappedReadFolder()+self._sra_acc+".bam")

    def hasBamFile(self):
        return self.has_file(self.getBamFile())

    def getSortedBamFile(self):
        return (self.getMappedReadFolder()+self._sra_acc+".sorted.bam")

    def hasSortedBamFile(self):
        return self.has_file(self.getSortedBamFile())

    def getBamStatsFile(self):
        return (self.getMappedReadFolder()+self._sra_acc+".bam.bas")

    def hasBamStatsFile(self):
        return self.has_file(self.getBamStatsFile())

    def getFilteredReadFolder(self):
        return (args.basedir+"filtered_reads/"+self._breed+"/")

    def getFilteredReadFile(self):
        return (self.getFilteredReadFolder()+self._sra_acc+".MT.bam")

    def hasFilteredReadFile(self):
        return self.has_file(self.getFilteredReadFile())

    def __unicode__(self):
        return ('Accession %s: %s %s %s' % (self._sra_acc, self._breed, self._sex, self._mbases))
        
    def __str__(self):
        return unicode(self).encode('utf-8')

    def __repr__(self):
        return str(self)

class BwaFactory(CommandFactory):
    '''Generate command line call to bwa
    '''
    def __init__(self, genome, n_threads=1, memory=8):
        CommandFactory.__init__(self, args.basedir, n_threads, memory, "map")
        self._genome = genome

    def make_cmd(self, s):
        ''' get the run command for the given sample

            The verbosity is set to 0 - no output. Excess verbosity (default 3) causes /tmp/slurmd.log overflow
            Setting -v from 0-2 has no effect on the logging to console - as documentation says, 'This option has not been fully supported throughout BWA'
        '''
        script_file = s.getMappedReadFolder()+s.getAccession()+".sh"
        bwa_cmd = ' bwa mem -v 1 -t %s %s "%s" "%s"  ' % (self.threads(), self._genome, s.getTrimmedRead1(), s.getTrimmedRead2()) #> "%s" removed due to test in sdout() below , s.getMappedReadFile()

        commands = []
        commands.append(bwa_cmd)
        commands.append("")

        self.write_script(commands, script_file, s)
        return("sh %s " % script_file)

    def has_output(self, s):
        return s.hasMappedReads() or s.hasBamFile()

    def has_input(self, s):
        return s.hasTrimmedReads()

    def lock_file(self, s):
        return s.getMappedReadFolder()+s.getAccession()+"."+self._log_suffix+".lck"

    def stdout(self, s):
        ''' WARNING: don't muck about with this sdout redirection.
            The following (conventional) command can drain a node by filling /tmp/slurmd.log:
            srun -c 8  --mem=40G  -e logfile.log  bwa mem -v 1 -t 8 genome.file fq1 fq2 > mapped.sam

            To prevent this, we redirect sdout via slurm to the SAM file:
            --mem=40G -o mapped.sam -e logfile.log  bwa mem -v 1 -t 8 genome.file fq1 fq2
        '''
        return s.getMappedReadFile()

    def __str__(self):
        return "read mapping"

class TrimFactory(CommandFactory):
    ''' Generate command line calls to trim galore
    '''

    def __init__(self, n_threads=1, memory=8):
        CommandFactory.__init__(self, args.basedir, n_threads, memory, "trim")

    def make_cmd(self, s):
        ''' get the trim command for the given sample
        '''
        cmd        = ' trim_galore --gzip --paired --fastqc --fastqc_args "--nogroup --extract" --output_dir "%s" "%s" "%s" ' % (s.getTrimmedReadFolder(), s.getFastqRead1(), s.getFastqRead2())
        commands   = []
        commands.append(cmd)
        commands.append("")

        script_file = s.getTrimmedReadFolder()+s.getAccession()+".sh"
        self.write_script(commands, script_file, s)
        return("sh %s " % script_file)

    def lock_file(self, s):
        return s.getTrimmedReadFolder()+s.getAccession()+"."+self._log_suffix+".lck"

    def has_output(self, s):
        return s.hasTrimmedReads() or s.hasMappedReads() or s.hasBamFile()

    def has_input(self, s):
        return s.hasFastqFiles()

    def __str__(self):
        return "trimming"

class DownloadFactory(CommandFactory):
    ''' Generate command line calls to download SRA fastq files
    '''

    def __init__(self, n_threads=1, memory=4):
        CommandFactory.__init__(self, args.basedir, n_threads, memory, "dl")

    def make_cmd(self, s):
        ''' get the download command for the given sample
        '''
        cmd        = "/usr/bin/Rscript %ssrc/downloadSRA.R %s %s '%s' " % (args.basedir, args.basedir, s.getAccession(), s._breed)
        commands   = []
        commands.append(cmd)
        commands.append("")

        script_file = s.getFastqFolder()+s.getAccession()+".sh"
        self.write_script(commands, script_file, s)
        return("sh %s " % script_file)

    def has_output(self, s):
        return s.hasFastqFiles() or s.hasTrimmedReads() or s.hasMappedReads() or s.hasBamFile()

    def has_input(self, s):
        return True # default, first step in pipeline

    def lock_file(self, s):
        return s.getFastqFolder()+s.getAccession()+"."+self._log_suffix+".lck"

    def __str__(self):
        return "downloading"

class BamCompressFactory(CommandFactory):
    ''' Generate command line calls to convert SAM files to BAM
    '''
    def __init__(self, n_threads=8, memory=40):
        CommandFactory.__init__(self, args.basedir, n_threads, memory, "bam")

    def make_cmd(self, s):
        ''' get the conversion command for the given sample
    
            view
                -b - compress to bam. 
                -h - include header
                -S - legacy identifing input as SAM.  Compatability with old versions
                -@ - the number of threads is fixed at 2; more than this has no improvement
            sort 
                -T - PREFIX for temporary files to avoid filename collision with parallel instances;
                     mapped_read_folder/sample_accession.nnnn.bam
                -@ - number of threads
                -o - output file
        '''

        filt = "samtools view -@ 2 -bSh %s " % (s.getMappedReadFile()) # compress to bam
        sort = 'samtools sort -T %s -@ %s -o %s ' % (s.getMappedReadFolder()+s.getAccession(), self.threads(), s.getBamFile())
        cmd  = "%s | %s" %(filt, sort)
        clean_sam = "rm %s" % s.getMappedReadFile()

        commands = (cmd, clean_sam)

        script_file = s.getMappedReadFolder()+s.getAccession()+".sh"
        self.write_script(commands, script_file, s)
        return("sh %s " % script_file)

    def has_output(self, s):
        return s.hasBamFile()

    def has_input(self, s):
        return s.hasMappedReads()

    def lock_file(self, s):
        return s.getMappedReadFolder()+s.getAccession()+"."+self._log_suffix+".lck"

    def stdout(self, s):
        return ""

    def __str__(self):
        return "compressing and sorting"

class ReadExtractionFactory(CommandFactory):
    ''' Generate command line calls to extract read pairs with one read in the MT
        and the other read elsewhere in the genome.
    '''
    def __init__(self, n_threads=8, memory=40):
        CommandFactory.__init__(self, args.basedir, n_threads, memory, "rex")

    def make_cmd(self, s):
        '''
        '''
        filt_sam_file     = s.getFilteredReadFolder()+s.getAccession()+".sam"
        header_file       = filt_sam_file+".header"
        mismatch_sam_file = s.getFilteredReadFolder()+s.getAccession()+".mismatch.sam"
        mt_temp_file      = s.getFilteredReadFolder()+s.getAccession()+".MT.tmp.sam"
        mt_sam_file       = s.getFilteredReadFolder()+s.getAccession()+".MT.sam"
        mt_temp_bam_file  = s.getFilteredReadFolder()+s.getAccession()+".MT.tmp.bam"
        mt_bam_file       = s.getFilteredReadFolder()+s.getAccession()+".MT.bam"
        
        filt            = "samtools view -F 14 %s > %s" % (s.getBamFile(), filt_sam_file) # fetch reads that are properly mapped, but not in a proper pair
        save_header     = "samtools view -H %s > %s" % (s.getBamFile(), header_file) # get the file header
        extract_chrs    = "awk '($3!=$7 && $7!=\"=\" && ( $3==\"MT\" || $7==\"MT\" )) ' %s > %s " % (filt_sam_file, mt_temp_file) # get the reads with mismatched chromosomes
        add_header      = "cat %s %s > %s " % (header_file, mt_temp_file, mt_sam_file) # add back the header
        convert_to_bam  = "samtools view -b %s > %s" % (mt_sam_file, mt_temp_bam_file) # convert to bam
        sort_bam        = "samtools sort -o %s %s" % (mt_bam_file, mt_temp_bam_file) # sort bam
        index_bam       = "samtools index %s" % (mt_bam_file) # index the bam file


        # Clean up the temporary files
        clean_tmp       = "rm %s" % mt_temp_file
        clean_header    = "rm %s" % header_file
        clean_sam       = "rm %s" % filt_sam_file
        clean_bam       = "rm %s" % mt_temp_bam_file
        clean_mt_sam    = "rm %s" % mt_sam_file

        commands = (filt, save_header, extract_chrs, add_header, convert_to_bam, sort_bam, index_bam, clean_tmp, clean_header, clean_sam, clean_bam, clean_mt_sam)

        script_file = s.getFilteredReadFolder()+s.getAccession()+".sh"
        self.write_script(commands, script_file, s)
        return("sh %s " % script_file)

    def has_output(self, s):
        return s.hasFilteredReadFile()

    def has_input(self, s):
        return s.hasBamFile()

    def lock_file(self, s):
        return s.getFilteredReadFolder()+s.getAccession()+"."+self._log_suffix+".lck"

    def stdout(self, s):
        ''' Override otherwise the output will go to the log.
            Using "" so srun will not divert from console.
        '''
        return ""

    def is_batch(self):
        ''' Should the command be run via srun or sbatch
        '''
        return False

    def __str__(self):
        return "extracting MT reads"

def create_logger():
    '''Configure the logger. 

    Use a log file to store all messages of level DEBUG 
    or higher. Send all messages of level INFO or higher
    to sdout.
    '''
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)


    pipeline.make_sure_path_exists(args.basedir+"logs/")
    logfile = datetime.now().strftime('analysis.%Y_%m_%d_%H_%M.log')
    fh = logging.FileHandler(args.basedir+"logs/"+logfile)
 
    formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(funcName)s\t%(levelname)s\t%(message)s')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG) # Log ALL THE THINGS to file
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO) # Log important things to console
    logger.addHandler(sh)

    # Add the log handlers to the pipeline logger
    pipe_logger = logging.getLogger("pipeline")
    pipe_logger.setLevel(logging.DEBUG)
    pipe_logger.addHandler(fh)
    pipe_logger.addHandler(sh)

    return(logger)

def read_sample_list():
    logger.info("Reading sample file")

    if( not os.path.isfile(args.sample_file)):
        logger.error("Expected sample file %s is not present" % args.sample_file)
        sys.exit(1)
    sampleList = []
    with open(args.sample_file, 'r') as f:
        for line in f:
            entry = line.strip().split("\t")
            entry = [w.replace('"', '') for w in entry] # replace all quotes
            entry = [w.replace(' ', '_') for w in entry] # replace all spaces
            if(entry[0]=="breed"): # skip header line
                continue
            sampleList.append( SraSample(entry[6], entry[0], entry[7], entry[4]) )
    logger.info("Read sample file of %d samples" % len(sampleList))
    return(sampleList)
                
def run_analysis():
    '''Create and run the analysis pipeline
    '''

    # Read sample file to get breed folder and sample ids
    samples = read_sample_list()

    pl = Pipeline(args.ninstances, args.test)
    pl.\
        add_method(instances=1, out_folder="fastq", runner=DownloadFactory(), is_skip=args.no_download).\
        add_method(out_folder="trimmed_reads/", runner=TrimFactory(n_threads=8, memory=40), is_skip=args.no_trim).\
        add_method(out_folder="mapped_reads/", runner=BwaFactory(genome="/mnt/cgs-fs3/Sequencing/Genome/Sscrofa11.1/bwa/Sus_scrofa.Sscrofa.11.1.dna.fa", n_threads=8, memory=40), is_skip=args.no_map).\
        add_method(out_folder="mapped_reads/", runner=BamCompressFactory(n_threads=8, memory=40), is_skip=args.no_convert).\
        add_method(out_folder="filtered_reads/", runner=ReadExtractionFactory(n_threads=1, memory=10), is_skip=args.no_extract).\
        run(samples)

    # add_method(out_folder="mapped_reads/", runner=BamSortFactory(n_threads=8, memory=40), is_skip=args.no_convert).\

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='run_analysis.py',description = '')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--basedir', 
        help="Base directory for analysis",
        default='./')
    parser.add_argument('--ninstances', 
        help="Number of parallel srun instances",
        default=3, type=int)
    parser.add_argument('--sample_file', 
        help="List of samples to process",
        default='./Filtered_samples.csv')
    parser.add_argument('--test', help="Run in test mode only", action='store_true')
    parser.add_argument('--no_download', help="Skip the download step", action='store_true')
    parser.add_argument('--no_trim', help="Skip the trimming step", action='store_true')
    parser.add_argument('--no_map', help="Skip the mapping step", action='store_true')
    parser.add_argument('--no_convert', help="Skip the bam conversion", action='store_true')
    parser.add_argument('--no_extract', help="Skip the MT read extraction", action='store_true')
    args   = parser.parse_args()

    logger = create_logger()
    logger.info("Invoked with arguments: "+str(args))

    logger.info("Running%sanalysis" % " test " if args.test else " ")
    run_analysis()
    logger.info("Analysis complete")
    sys.exit(0)
