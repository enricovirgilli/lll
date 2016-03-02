#!/usr/bin/env python

"""USAGE:

OPTIONS:
    -h, -?, --help: obtain help

TODO:
"""
import sys, getopt
import XCOM

# checks options
def checkopt():
    flags={}
    sopt='h?'
    lopt=['help']
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], sopt, lopt)
    except getopt.GetoptError:
        sys.stderr.write(__doc__)
    for opt, optarg in opts:
        if opt in ('-h', '-?', '--help'):
            sys.stdout.write(__doc__)
            sys.exit()
    return opts, args, flags

opts, args, flags=checkopt()

if __name__=='__main__':
    # Import Psyco if available
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    # ...your code here...
    print "\n".join(["%s\t%s" %(e, XCOM.main(29, e/1000.)) for e in range(1, 1000, 100)])
