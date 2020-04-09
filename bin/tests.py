import os,sys
import subprocess
import shutil
from distutils.spawn import find_executable


def find_repair(repair_file):
    fC_path = find_executable('featureCounts')
    potential_paths = [os.path.join(fC_path, os.pardir, 'utilities', 'repair'),
                       os.path.join(fC_path, os.pardir, 'repair'),
                       '/usr/lib/subread/repair']
    found = 0
    for p in potential_paths:
        if os.path.exists(p):
            with open(repair_file, 'w') as w:
                w.write(p)
            found = 1
            break
    return found


def check_dependencies(run_mode):
    if run_mode=='correct_annotation':
        dependencies=['samtools','bedtools']
    elif run_mode == 'correct_count':
        dependencies = ['samtools', 'bedtools', 'featureCounts']
    elif run_mode == 'correct_bedgraph':
        dependencies = ['samtools', 'bedtools', 'pairedBamToBed12', 'featureCounts', 'repair']
    errors=[]
    for dependency in dependencies:
        exist=find_executable(dependency)
        if exist==None and dependency != 'repair':
            errors.append(dependency)
        elif dependency == 'repair':
            if 'featureCounts' in errors:
                pass
            else:
                repair_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'repair_path.txt')
                if os.path.exists(repair_file):
                    with open(repair_file, 'r') as f:
                        repair_path = f.readline().strip()
                        if os.path.exists(repair_path) is True:
                            found = 1
                        else:
                            found = find_repair(repair_file)

                else:
                    found = find_repair(repair_file)
                if found == 0:
                    errors.append('repair')
                    print('Subread\'s repair executable was not found in standard paths. Please write the path in :'
                          '%s' %repair_file)
    if len(errors) != 0:
        print('error: %s is/are not installed. Please read the README.md file '
              'for more information about the prerequisites of CoCo and help for their installation.' %(', '.join(errors)))
        sys.exit(1)


def check_bam(bam_file):
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    if is_binary_string(open(bam_file, 'rb').read(1)) == False:
       print('error: File is not binary. Please use a binary alignment map (.bam) format as input.')
       sys.exit(1)
    command="samtools view %s | head -n1" %(bam_file)
    view=subprocess.Popen([command], shell=True, stdout=subprocess.PIPE)
    communicate=view.communicate()
    line=str(communicate[0],'utf-8')
    if line == '':
        print('error: samtools view cannot proccess provided bam file %s , a proper bam file is required.' %(bam_file))
        sys.exit(1)
    else:
        linesplit=line.split('\t')
        if len(linesplit) < 11:
            print("first line :")
            print(line)
            print("error: entries from bam file are incomplete. the file provided should have at least 11 fields, only %d fields in the bam's first entry." %(len(linesplit)))
            sys.exit(1)


def check_gtf(gtf_file,check_gene_biotype=False):
    with open(gtf_file,'r') as f:
        for line in f:
            if line.startswith('#')==False:
                if line.split('\t')[2] == 'gene':
                    gene_test_line=line
                elif line.split('\t')[2] == 'transcript':
                    transcript_test_line=line
                elif line.split('\t')[2] == 'exon':
                    exon_test_line=line
                try:
                    gene_test_line
                    transcript_test_line
                    exon_test_line
                except:
                    pass
                else:
                    break
    #Check if any gene, transcript and exon entry was found in the gtf file.
    try:
        gene_test_line
    except:
        print('error: bad annotation format, no gene entry. Please use a gene transfer format (.gtf) as input annotation.')
        sys.exit(1)
    try:
        transcript_test_line
    except:
        print('error: bad annotation format, no transcript entry. Please use a gene transfer format (.gtf) as input annotation.')
        sys.exit(1)
    try:
        exon_test_line
    except:
        print('error: bad annotation format, no exon entry. Please use a gene transfer format (.gtf) as input annotation.')
        sys.exit(1)

    #Check if there is any gene_biotype entry
    if check_gene_biotype:
        if 'gene_biotype' not in gene_test_line:
            print('error: There is no "gene_biotype" entry for genes. Please use a gene transfer format (.gtf) annotation file obtained from Ensembl to correct this.')
            sys.exit(1)
    if 'gene_id' not in gene_test_line:
            print('error: There is no "gene_id" entry for genes. Please use a gene transfer format (.gtf) annotation file obtained from Ensembl to correct this.')
            sys.exit(1)


def check_output(output):
    if os.path.isdir(output):
        print('error: output specified is a directory, please provide a file name to be created')
        sys.exit(1)
    output_path='/'.join(output.rstrip('/').split('/')[:-1])
    cwd = os.getcwd()
    if os.path.isdir(output_path)==False:
        try:
            if os.path.isdir(output_path)== False:
                if os.path.isdir(cwd+'/'+output_path)== False:
                    sys.exit(1)
        except:
            print("error: path to output does not exist: %s. Please specify a valid output path." %(output))
            sys.exit(1)

def check_unique_count(unique_counts):
    with open(unique_counts) as f:
        line = f.readline()
        while line[0]=='#':
            line = f.readline() #skips comment line
        f.readline() #Skips header if there is one, otherwise, skips first line, which is fine
        line = f.readline()
        line=line.strip().split('\t')
        if len(line) != 7:
            print('Error: Wrong number of fields in unique_count matrix, seen %i, expected 7. Format should be: '
                  'gene_id  seqname start   end strand  length  accumulation (separated by a tabulation)'%len(line))
            sys.exit(1)
        try :
            float(line[6])
        except ValueError:
            print('Error: Last field should be a number')
            sys.exit(1)
