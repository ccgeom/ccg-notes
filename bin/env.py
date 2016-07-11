# -*- coding: utf-8 -*-

import sys
import os
import virtualenv as venv

"""
Colorful output
"""

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = "\033[1m"

def head(msg):
    print (HEADER + msg + ENDC)

def info(msg):
    print msg

def infog(msg):
    print (OKGREEN + msg + ENDC)

def infob(msg):
    print (OKBLUE + msg + ENDC)

def warn(msg):
    print (WARNING + msg + ENDC)

def err(msg):
    print (FAIL + msg + ENDC)

"""
Check python version
"""

info("checking python version...")

cur_version = sys.version_info

if cur_version < (2,7):
    err("Your Python interpreter is too old. Please consider upgrading.")
    sys.exit(-1)

if cur_version >= (3,0):
    err("Your Python interpreter is 3.x but 2.7 is required.")
    sys.exit(-1)

"""
Check virtual enviroment
"""

if not os.path.exists(".py"):
    sys.argv = ['virtualenv', '--system-site-packages', '.py']
    venv.main()