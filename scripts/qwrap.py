#!/usr/bin/env python
# encoding: utf-8
"""
qwrap.py

Copyright (c) 2008, Michael Kuhn
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the European Molecular Biology Laboratory nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY Michael Kuhn ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Michael Kuhn BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import sys
import os
import re
import time
import itertools

from collections import defaultdict

def run_submission(command):
    """Run the command passed as parameter and sniff the list of job ids from the output."""
    
    job_list = []
    
    re_job_id = re.compile(r"Your job (\d*) \(.*\) has been submitted").match
    
    for line in os.popen(command):
        print line,
        m = re_job_id(line)
        if m:
            job_list.append(m.group(1))
        
    return job_list
        
        
def monitor_jobs(job_list):
    """Monitor the progress of the passed jobs. Exit when none of the jobs remain in the queueing system."""
    
    jobs = set(job_list)
    
    n_submitted = len(jobs)
    
    # you mave have to adapt this your queuing system: Where does 'qstat' indicate the status of the job?
    state_index = 40
    status_line = ""
    
    for spinner in itertools.cycle("-\\|/"):
        
        status_counts = defaultdict(int)
        
        for line in os.popen("qstat"):
            fields = line.strip().split()
            
            if fields[0] in jobs:
                status_counts[ line[state_index:state_index+3].strip() ] += 1
            
        if not status_counts:
            break
            
        sys.stderr.write("%s\r" % (" "*(10+len(status_line))))    
            
        status_line = " ".join( [ "%3s: %3d" % (k, v) for (k, v) in sorted(status_counts.items()) ] )
            
        sys.stderr.write("%s Status: %s\r" % (spinner, status_line))    
        
        time.sleep(5)

    print >> sys.stderr, ""
    print >> sys.stderr, "No jobs remaining."
    


def main():
    
    if len(sys.argv) < 2:
        sys.exit("Usage: qwrap.py command\nwhere 'command' is a script or program that submits several jobs")
        
    job_list = run_submission(" ".join(sys.argv[1:]))
    
    if not job_list:
        sys.exit("Error: no job submissions were detected.")
    
    print >> sys.stderr, ""
    print >> sys.stderr, "Scheduling finished. %d jobs have been detected." % len(job_list)
    
    monitor_jobs(job_list)
        


if __name__ == '__main__':
    main()
