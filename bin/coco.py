#!/usr/bin/env python3
"""
coco.py
This script is called by the bash executable "coco" from the /bin/ directory.
It calls one of three other scripts: correct_annotation, CorrectCount or CorrectBedgraph depending on the parameters used by the user.
"""
__author__      = "Vincent Boivin, Gabrielle Deschamps-Francoeur, and Michelle Scott"
__email__ = "Michelle.Scott@Usherbrooke.ca"
__version__='0.1.0'

import sys
import os
import tests

help="CoCo: Count Corrector for embedded and multi-mapped genes.\n" \
              "Version: "+str(__version__)+"\n" \
              "Usage: CoCo <Run mode> <Run mode specific arguments>\n" \
              "\n" \
              "Run modes:\n" \
              "\tcorrect_annotation, CA, ca:\tProduce modified annotation for embedded genes.\n" \
              "\tcorrect_count, CC, cc:\t\tProduce gene expression values, taking multi-mapped reads into account.\n" \
              "\tcorrect_bedgraph, CB, cb:\tProduce paired-end bedgraphs.\n" \
              "\n" \
              "-h, --help\tshow this help message and exit\n" \
              "-v, --version\tshow CoCo version\n" \
              "\n" \
              "For more info, please refer to the README file."
if len(sys.argv)==1:
    print("error: A run mode must be selected for CoCo. (choose from 'CorrectAnnotation', 'CorrectCount', 'CorrectBedgraph')")
    print(help)
    sys.exit(1)
else:
    tests.check_dependencies()
    run_mode=sys.argv[1]
    run_mode_dict={'CA':'correct_annotation','ca':'correct_annotation','CC':'correct_count','cc':'correct_count',
                   'CB':'correct_bedgraph','cb':'correct_bedgraph'}
    if run_mode in run_mode_dict:
        run_mode=run_mode_dict[run_mode]
    if run_mode in ['correct_annotation','correct_count','correct_bedgraph']:
        print('Launching',run_mode)
        os.system('python3 '+(os.path.realpath(__file__).replace('coco.py',''))+run_mode+'.py '+(' '.join(sys.argv[2:])))
    elif run_mode=='--help' or run_mode=='-h':
        print(help)
        sys.exit()
    elif run_mode=='--version' or run_mode=='-v':
        print('CoCo version :',__version__)
        sys.exit()
    else:
        print("error: the run mode %s is invalid. (choose from 'correct_annotation', 'correct_count',"
              " 'correct_bedgraph')" %(run_mode))
        print(help)
        sys.exit(1)

