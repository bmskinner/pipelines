from collections import OrderedDict
from datetime import datetime
import errno
import logging
import multiprocessing
import os
import os.path
import Queue
import signal
import subprocess
import sys
import time

class Sample:
    '''Abstract base class to hold sample data
    '''

    def __init__(self, name):
    	'''Default constructor.

    		--name - the sample name
    	'''
        self._name = name

    def name(self):
        return self._name

    def log_file(self, suffix):
        return ""

    def has_file(self, file):
        return os.path.isfile(file) and os.path.getsize(file)>0

    def __unicode__(self):
        return ('Sample %s' % (self._name))
        
    def __str__(self):
        return unicode(self).encode('utf-8')

    def __repr__(self):
        return str(self)


class CommandFactory:
    '''Abstract base class for the command factories that
       generate command line calls for Samples
    '''

    def __init__(self, base_dir="./", n_threads=1, memory=8, suffix="null"):
    	self._base_dir   = base_dir 
        self._n_threads  = n_threads # number of threads per instance of the method
        self._memory     = memory       # memory required per instance
        self._log_suffix = suffix   # suffix of the sample log file (e.g. sample_name.suffix.log)

    def write_script(self, commands, script_file, s):
        ''' Write the given list of commands to a shell script.

            Includes commands to create a lock file for the current factory and sample.
            The lock file and the script itself are deleted once execution is complete
        '''
        with open(script_file, 'w') as f:
            f.write( "#!/bin/sh\ntouch '%s'\n%s\nrm '%s'\nrm '%s'\n" % ( self.lock_file(s), "\n".join(commands), self.lock_file(s), script_file )   )
        f.close()

    def has_output(self, s):
        ''' Does the output for the sample already exist?

            This must be overridden by subclasses.
        '''
        return False

    def has_input(self, s):
        ''' Does the input for the sample already exist?

            This must be overridden by subclasses.
        '''
        return False

    def lock_file(self, s):
        ''' Define the lock file for a sample.

            This must be overridden by subclasses.
        '''
        return ""

    def has_lockfile(self, s):
        ''' Is there a lock file for the sample, indicating leftover from an
            interrupted analysis?
        '''
        return os.path.isfile(self.lock_file(s))

    def make_cmd(self, s):
        ''' get the run command for the given sample

            This must be overridden by subclasses.
        '''
        return ""

    def base_dir(self):
    	'''Get the base directory for the analysis
    	'''
    	return self._base_dir

    def threads(self):
        '''Get the thread requirement for the command
        '''
        return self._n_threads

    def memory(self):
        '''Get the memory requirement for the command
        '''
        return self._memory

    def logfile(self, s):
        ''' Get the sample log file for the command
        '''
        return self._base_dir + ("logs/%s.%s.log" % (s.name(), self._log_suffix))

    def stdout(self, s):
        '''Default output file for stdout. Use "" for console.
        '''
        return self.logfile(s)

    def stderr(self, s):
        '''Default output file for stderr. Use "" for console.
        '''
        return self.logfile(s)

    def exclusive(self):
        ''' Override this to make the job take exclusive control of the node,
            whatever the resource requirement
        '''
        return False

    def is_batch(self):
        ''' Should the command be run via srun or sbatch
        '''
        return False

    def job_name(self, s):
        return s.name()+"_"+self.__class__.__name__;


class AbstractParallelMethod():
    ''' Class to perform parallelised methods
        on samples.

        When jobs are submitted to the cluster, ensures not 
        all the nodes are filled. Maintains an internal task queue
        which keeps <ninstances> srun jobs populated.
    '''  
    def __init__(self, instances=1, test=False, method=None, out_folder=None, runner=None):
    	''' Default constructor

    		--instances  - the number of parallel instances in the pool
    		--test       - if the method should be run in test mode
    		--method     - the method to be invoked (defaults to internal method which submits tasks to a slurm cluster via srun)
    		--out_folder - the ouput folder for the analysis
    		--runner     - the CommandFactory that will generate commands to be executed for each sample
    	'''
        self._is_test = test
        self._ninstances = instances
        if(method==None):
            self._base_method = self._run
        else:
            self._base_method = method
        # self._base_method = method if method not None else self._run
        self._out_folder = out_folder
        self._runner = runner

    def create_output_folders(self, samples):
        ''' Create a sub folder for each unique breed in the samples

            --samples - the samples
        '''
        logger = logging.getLogger(__name__)
        top_level_dir = slash_terminate(self._runner.base_dir()+self._out_folder)
        logger.debug("Creating output folders for directory %s" % top_level_dir)

        if(self._is_test):
            logger.debug("Would create output folder %s" % top_level_dir)
        else:
            make_sure_path_exists(top_level_dir)
        
        breeds = set(map(lambda s: s._breed, samples)) # get distinct breeds
        for b in breeds:
            b_folder = top_level_dir+b+"/"
            if(self._is_test):
                logger.debug("Would create output folder %s" % b_folder)
            else:
                make_sure_path_exists(b_folder)

    def worker_function(self, job_queue):
        ''' A worker for the substitute multiprocessing pool

            Takes a sample from the job queue and executes the base method (default _run) using the command factory

            --job_queue - the queue to process
        '''
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        while not job_queue.empty():
            try:
                s = job_queue.get(block=False)
                self._base_method(s, self._runner)
            except Queue.Empty:
                pass

    def run(self, samples):
        ''' Run the method.

            A worker pool is constructed to handle the samples, and invoked through srun
            until all are processed.
        '''
        self.create_output_folders(samples)
        logger = logging.getLogger(__name__)
        
        # NEW METHOD
        # Given the problem with the old method shown below, here we implement a custom 
        # job queue and worker pool that can respond to keyboard interrupts. It's 
        # not perfect - anything sent to slurm before the interrupt will keep running,
        # but no new jobs will be submitted, and the script will exit.
        job_queue = multiprocessing.Queue()

        logger.info("Running %s on %d samples" % (str(self._runner), len(samples)))
        logger.debug("Invoking %d workers" % self._ninstances)

        for s in samples:
            job_queue.put(s)

        workers = []
        for i in range(self._ninstances):
            tmp = multiprocessing.Process(target=self.worker_function,
                                          args=(job_queue,))
            tmp.start()
            workers.append(tmp)

        try:
            for worker in workers:
                worker.join()
            logger.debug("Worker pool closed normally")
        except KeyboardInterrupt:
            for worker in workers:
                worker.terminate()
                worker.join()
            logger.warn("Script interrupted; manual file cleanup may be needed")
            logger.warn("Continuing jobs already submitted to slurm")
            sys.exit(1)
        logger.info("Finished %s" % str(self._runner))

        # OLD METHOD using a multiprocessing pool
        # There is an issue with using a worker pool; worker processes can handle the KeyboardInterrupt and call sys.exit, 
        # but the processes persist and still receive future tasks. 
        # Additionally, the KeyboardInterrupt is not delivered to the parent process until all jobs are completed.
        # This means Ctrl-C will only add the next worker into slurm; it will not cancel all the jobs
        # See here for possible solution: https://bryceboe.com/2010/08/26/python-multiprocessing-and-keyboardinterrupt/

        
        # _pool = multiprocessing.Pool(processes=self._ninstances) # pool must be global
        # logger.info("Running %s on %d samples" % (str(self._runner), len(samples)))
        # logger.debug("Invoking apply_async with %d processes" % self._ninstances)

        # try:
        #     for s in samples:
        #         _pool.apply_async(self._base_method, args=(s,self._runner))
        # except KeyboardInterrupt:
        #     logger.debug("Caught KeyboardInterrupt, terminating workers and quitting")
        #     _pool.terminate()
        #     logger.warn("Script interrupted; manual file cleanup may be needed")
        #     sys.exit(1)
        # else:
        #     _pool.close()
        #     _pool.join()
        #     logger.debug("Worker pool closed normally")

    def _run(self, s, runner):
        ''' Run the commands generated by a command factory for the given sample

            --s      - the Sample to run
            --runner - the CommandFactory
        '''
        logger = logging.getLogger(__name__)
        try:
            if(runner.has_output(s) and not runner.has_lockfile(s)):
                logger.debug("Output files of %s exists for %s, skipping" % (str(runner), s.name()))
                return
            if(not runner.has_input(s)):
                logger.debug("Missing input files of %s for %s, skipping" % (str(runner), s.name()))
                return
            if(runner.has_lockfile(s)):
                logger.debug("Lock file of %s exists for %s, repeating" % (str(runner), s.name()))

            start = time.time()
            sdout_str = (", directing sdout to %s" % runner.stdout(s)) if runner.stdout(s) else ""
            sderr_str = (", directing sderr to %s" % runner.stderr(s)) if runner.stderr(s) else ""
            logger.info("Running %s on sample %s%s%s" % (str(runner), s.name(), sdout_str, sderr_str) )
            if(self._is_test):
                logger.debug("Would run: %s" % runner.make_cmd(s))
            else:
                if(runner.is_batch()):
                    sbatch(cmd_list=runner.make_cmd(s), sdout=runner.stdout(s), sderr=runner.stderr(s))
                else:
                    srun(cmd=runner.make_cmd(s), ncores=runner.threads(), mem=runner.memory(), sdout=runner.stdout(s), sderr=runner.stderr(s), exclusive=runner.exclusive(), job_name=runner.job_name(s))
            end = time.time()
            logger.info("Completed %s %s in %s" % (str(runner), s.name(), format_seconds(end-start)))
            return
        except:
            e = sys.exc_info()[0]
            logger.exception("Error running %s on sample %s" % (str(runner), s.name()) )


class Pipeline:
    '''Simplify running of each method in the analyis by chaining them into a pipeline
    '''

    def __init__(self, instances=1, test=False):
    	'''Default constructor

    		--instances - the default number of parallel instances.
    		--test      - true if the pipeline should be run in test mode, false otherwise 
    	'''
        self._ninstances = instances
        self._is_test = test
        self._methods = OrderedDict() # store methods in the order they were inserted

    def add_method(self, instances=None, out_folder=None, runner=None, is_skip=False):
        '''Add a method to the pipeline

			--instances  - the number of parallel instances. Overrides the defaults from the constructor.
			--out_folder - the output folder for the analysis
			--runner     - the CommandFactory to generate commands
            --is_skip    - true if the method should be skipped, false otherwise

            returns - the pipeline
        '''
        logger = logging.getLogger(__name__)

        if(instances==None):
        	instances = self._ninstances

        method = AbstractParallelMethod(instances=instances, test=self._is_test, out_folder=out_folder, runner=runner)
        self._methods[method] = is_skip
        logger.debug("Added %s to pipeline; skip is %s" % (runner.__class__.__name__, is_skip))
        return(self)

    def run(self, samples):
        '''Run the methods in the order they were added, if the
           method is not set to skip

           --samples - the Sample collection to run the pipeline on
        '''
        logger = logging.getLogger(__name__)
        for method in self._methods.keys():
            if(self._methods[method]):
                logger.info("Skipping %s" % method._runner.__class__.__name__)
            else:
                method.run(samples)


def format_seconds(sec):
    days = sec // 86400
    hours = sec // 3600 % 24
    minutes = sec // 60 % 60
    seconds = sec % 60
    return("%dd %dh %dm %ds" % (days, hours, minutes, seconds))

def runAndCheck(cmd, msg):
    '''Run the given command and check the return value.

        If the return value of the command is not 0 (success),
        print the given error message and return
        exit code 1.

        --cmd - the command to be run
        --msg - the message to print if the command returns an error
    '''
    logger = logging.getLogger(__name__)
    logger.debug("Invoking: " + cmd)
    try:
        retcode = subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
        if retcode != 0:
            logger.debug("System call returned "+str(retcode))
            logger.error(msg+". Quitting.")
        return(retcode)
    except OSError as e:
        logger.exception("Error invoking command")
        print >>sys.stderr, "Execution exception:", e
        return(1)

def srun(cmd, ncores="1", mem="", sdout="", sderr="", exclusive=False, job_name=""):
    '''Run the given command through srun.

    If srun is available, the command will be run
    with the given cores and memory. If srun is not
    installed, the command is invoked directly.

    --cmd       - the command to be run (required)
    --ncores    - the number of cores (optional)
    --mem       - the memory in Gb (optional)
    --sdout     - the desired output file for sdout redirection (optional)
    --sderr     - the desired output file for sderr redirection (optional)
    --exclusive - if true, the exclusive flag will be set; only this job will use the node (optional)
    --job_name  - the name for the job to display in sview / squeue (optional)
    '''
    logger = logging.getLogger(__name__)
    if(srun_is_installed()):
        s = "srun -c %s " % (str(ncores))
        s = s + " --mem="+str(mem)+"G " if mem else s
        s = s + (" -o %s " % sdout) if sdout else s
        s = s + (" -e %s " % sderr) if sderr else s
        s = s + " --exclusive " if exclusive else s
        s = s + (" -J %s " % job_name) if job_name else s
        cmd = s+cmd
    else:
        logger.debug("srun not found")

    return runAndCheck(cmd, "Error in srun command") 

def sbatch(cmd_list, sdout="", sderr=""):
    sdout = sdout if sdout else "/dev/null"
    sderr = sderr if sderr else "/dev/null"
    run = "sbatch -o %s -e %s <<EOF\n#!/bin/sh\n%s\nEOF" % (sdout, sderr, "\n".join(cmd_list))
    return runAndCheck(run, "Error in sbatch command") 


def srunRscript(cmd, ncores=1, mem="", sdout="", sderr=""):
    ''' Run the given R script via srun.

    This is a wrapper to functions::srun, which simply adds
    the Rscript executable to the front of the command.

    --cmd - the command to be run
    --ncores - the number of cores
    --mem - the memory in Gb
    --log_file - the desired output file for sdout and sderr redirection
    '''
    logger = logging.getLogger(__name__)
    cmd = '/usr/bin/Rscript '+ cmd
    logger.debug("Invoking Rscript: %s " % cmd)
    srun(cmd, ncores, mem, sdout, sderr)


def testLogging(logger):
    '''Testing for the logger
    ''' 
    logger = logging.getLogger(__name__)
    logger.info("Testing logging info level")
    logger.debug("Testing logging debug level")
    logger.error("Testing logging error level")
    try:
        logger.info("Testing stack trace")
        raise RuntimeError
    except Exception, err:
        logger.exception("Expected exception")

def srun_is_installed():
    '''Tests if srun is installed. 

        Uses 'command -v' for POSIX compliance; works in sh and bash.
        When srun is not found, the result will be 1, else null.
    '''
    cmd = 'command -v srun >/dev/null 2>&1 || { echo "1" >&2; }'
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    output = p.stdout.read()
    return(output != "1\n")

def slash_terminate(s):
    '''Add a trailing '/' to the string if not present
    '''
    return(s if s.endswith('/') else s + '/')

def make_sure_path_exists(path):
    ''' Create the given path if it does not exist
    '''
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise