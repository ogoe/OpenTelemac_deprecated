#!/usr/bin/env python
import sys
from os import path,environ,system
PWD = path.dirname(path.dirname(sys.argv[0]))
CODE = path.splitext(path.basename(sys.argv[0]))[0]
PYTEL = path.join(PWD,'pytel')
SYSTELCFG = path.join(PWD,'config')
if environ.has_key('SYSTELCFG'): SYSTELCFG = environ['SYSTELCFG']
system(path.join(PYTEL,'runcode.py')+' -f '+SYSTELCFG+' '+CODE+' '+' '.join(sys.argv[1:]))
