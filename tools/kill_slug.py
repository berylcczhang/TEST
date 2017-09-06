"""
This is a little script that just kills off all slug processes on a linux-like system
"""

import subprocess
import os
import signal

pl = subprocess.Popen(['ps', '-A'], stdout=subprocess.PIPE).communicate()[0].splitlines()
slugproc = [p for p in pl if p.endswith('slug')]
slugpid = [p.split()[0] for p in slugproc]
for p in slugpid:
    os.kill(int(p), signal.SIGTERM)
