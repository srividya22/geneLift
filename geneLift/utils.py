
"""
Provides utility functions for handler/taskrunner level code to provide consistent
behavior across different taskrunners and reduce code repetition of similar tasks.
"""

from __future__ import (division, print_function, absolute_import,unicode_literals)
import sys,os,shutil
import subprocess
#import simplejson
import logging 
import time
import operator

def cleanup(logger=None, directory=None):
    """
    Clean up after the job.  At the moment this just means removing the working
    directory, but later could mean other things.
    """
    
    try:
        shutil.rmtree(directory, ignore_errors=True) # it would not delete if fold is not empty
        # need to iterate each entry
    except IOError as e:
        logger.error("Unable to remove working directory {0}".format(directory))
        raise

def setupWorkingDir(logger=None,directory=None):
    """
	Clean up an existing workingdir and create a new one
    """
    try:
        if os.path.exists(directory): cleanup(logger,directory)
        os.makedirs(directory)
    except IOError as e:
        logger.error("Unable to setup working dir {0}".format(directory))
        raise

def run(cmnd):
    """ Run command and report status. """
    log('Running : %s' % cmnd)
    if subprocess.call(cmnd, shell=True, executable='/bin/bash') != 0:
       raise RuntimeError('Failed : %s ' % cmnd)


def log(message):
    """ Log messages to standard output. """
    print(time.ctime() + ' --- ' + message, flush=True)
