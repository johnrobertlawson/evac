"""Run this in the evac top level directory
        after committing code changes.
"""

import subprocess
import os

evacdir = os.getcwd()
# evacdir = os.path.join(os.getcwd(),'evac')
rootdir = os.path.dirname(evacdir)
htmldir = os.path.join(rootdir,'evac_docs','html')

def cmd(s,d):
    os.chdir(d)
    os.system(s)
    # p = subprocess.Popen(s.split(' '),cwd=d)
    # p.wait()
    return

# Commit changes
cmd('git add -A .',evacdir)
cmd('git commit -m "Updates docs"',evacdir)
cmd('git push',evacdir)
   
# Make the html
cmd('make html 2>&1 | tee log.log',os.path.join(evacdir,'docs'))

# Git add
cmd('git add -A .',htmldir)
cmd('git commit -m "Doc html push"',htmldir)
cmd('git push',htmldir)

