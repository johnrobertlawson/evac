"""Run this in the evac top level directory
        after committing code changes.
"""

import subprocess
import os

evacdir = os.path.join(os.getcwd(),'evac')
rootdir = os.path.dirname(evacdir)
htmldir = os.path.join(rootdir,'evac_docs','html')

def cmd(s,d):
    p = subprocess.Popen(s,cwd=d)
    p.wait()
    return
   
# Make the html
cmd('make html',os.path.join(evacdir,docs))

# Git add
cmd('git add .',html)
cmd('git commit -m "Docs update"',html)
cmd('git push origin gh-pages',html)

