#!/usr/bin/env python3
import sys, os

version='0.1.0'
help="CoCo: Count Corrector for embedded and multi-mapped genes.\n" \
              "Version: "+str(version)+"\n" \
              "Usage: CoCo <Run mode> <Run mode specific arguments>\n" \
              "\n" \
              "Run modes:\n" \
              "\tCorrectAnnotation:\tProduce modified annotation for embedded genes.\n" \
              "\tCorrectCount:\t\tProduce gene expression values, taking multi-mapped reads into account.\n" \
              "\tCorrectBedgraph:\tProduce paired-end bedgraphs.\n" \
              "\n" \
              "-h, --help\tshow this help message and exit\n" \
              "\n" \
              "For more info, please refer to the README file."
if len(sys.argv)==1:
    print("error: a run mode must be selected for CoCo. (choose from 'CorrectAnnotation', 'CorrectCount', 'CorrectBedgraph')")
    print(help)
    sys.exit(1)
else:
    run_mode=sys.argv[1]
    if run_mode in ['CorrectAnnotation','CorrectCount','CorrectBedgraph']:
        print('Launching',run_mode)
        os.system('python3 '+(os.path.realpath(__file__).replace('CoCo.py',''))+run_mode+'.py '+(' '.join(sys.argv[2:])))
    elif run_mode=='--help' or run_mode=='-h':
        print(help)
        sys.exit()
    else:
        print("error: the run mode %s is invalid. (choose from 'CorrectAnnotation', 'CorrectCount', 'CorrectBedgraph')" %(run_mode))
        print(help)
        sys.exit(1)

