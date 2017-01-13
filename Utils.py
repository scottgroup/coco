import re
import os
from statistics import mean
import pandas as pd

def line_prepender(filename, line):
    f = open(filename, 'r')
    file_string=f.read()
    file_string=line+file_string
    o=open(filename,'w')
    o.write(file_string)
    o.close()

def print_CD_box_intervals(sno_seq):
    '''
    Returns all possible C or D boxes for a sequence
    :param sno_seq: sequence to investiguate.
    :return: list of re.finditer objects for C boxes and D boxes. Attributes are .span() and .match().
    '''
    C_boxes=re.finditer('[GATC][GAT]GA[TG]G[ATG]', sno_seq)
    D_boxes=re.finditer('(CTGA)', sno_seq)
    #print([ x.group()+':'+str(x.span()) for x in list(D_boxes) ])
    C_boxes=list(C_boxes)
    D_boxes=list(D_boxes)
    if len(C_boxes) >= 1:
        C_boxes = [x for x in C_boxes if x.span()[-1] < len(sno_seq)/2]
    return C_boxes, D_boxes

def print_ACA_box_intervals(sno_seq):
    '''
    Returns the more probable H and ACA boxes for a sequence.
    Finds the H box closest to the middle of the sequence
    Finds an ACA box that must be at most 5 nucleotides away from the end of the sequence.
    If no ACA box is found, no H box is given.
    If no H box is found (quite unlikely), the ACA box is given anyway.
    :param sno_seq: sequence to investiguate.
    :return: re.finditer objects for H box and ACA box. Attributes are .span() and .match().
    '''
    length=len(sno_seq)
    middle=length/2
    H_boxes=re.finditer('A[GATC]A[GATC]{2,3}?A', sno_seq)
    ACA_boxes=re.finditer('ACA', sno_seq)
    #print([ x.group()+':'+str(x.span()) for x in list(D_boxes) ])
    H_boxes=list(H_boxes)
    ACA_boxes=list(ACA_boxes)

    try:
        ACA_box = ACA_boxes[-1]
        if int(ACA_box.span()[-1]) < length - 5:    #since we want the ACA box near the very end...
            ACA_box = None
    except:
        ACA_box = None
    H_box = None
    if ACA_box != None:
        if len(H_boxes) >= 1:
            H_boxes = [x for x in H_boxes if x.span()[-1] < ACA_box.span()[0]]  #Removes H boxes that would overlap with ACA box
            if len(H_boxes) >= 1:
                H_box = min(H_boxes, key=lambda x:abs(mean((x.span()))-middle)) #Takes the H box closest to the middle.
    return H_box, ACA_box


def rainbow_dict_maker(size,min=0,type='RGB'):
    if size%4 != 0:
        size=size-1
        if size%4 != 0:
            size=size-1
            if size%4 != 0:
                size=size-1
    size_quarter=size/4
    color_jump=250/(size/4)
    position=0
    RGB=[250,0,0]
    rainbow_dict={}
    if type=='RGB':
            rainbow_dict[1]=[ str(int(round(value, 0))) for value in RGB[:] ]
    elif type=='hex':
        rainbow_dict[1+min]=('#%02x%02x%02x' % tuple([ int(round(value, 0)) for value in RGB[:] ]))
    while position < size_quarter:
        RGB[1]+=color_jump
        if type=='RGB':
            rainbow_dict[position+2+min]=[ str(int(round(value, 0))) for value in RGB[:] ]
        elif type=='hex':
            rainbow_dict[position+2+min]=('#%02x%02x%02x' % tuple([ int(round(value, 0)) for value in RGB[:] ]))
        position+=1
    while position < (size_quarter*2):
        RGB[0]-=color_jump
        if type=='RGB':
            rainbow_dict[position+2+min]=[ str(int(round(value, 0))) for value in RGB[:] ]
        elif type=='hex':
            rainbow_dict[position+2+min]=('#%02x%02x%02x' % tuple([ int(round(value, 0)) for value in RGB[:] ]))
        position+=1
    while position < (size_quarter*3):
        RGB[2]+=color_jump
        if type=='RGB':
            rainbow_dict[position+2+min]=[ str(int(round(value, 0))) for value in RGB[:] ]
        elif type=='hex':
            rainbow_dict[position+2+min]=('#%02x%02x%02x' % tuple([ int(round(value, 0)) for value in RGB[:] ]))
        position+=1
    while position < (size_quarter*4):
        RGB[1]-=color_jump
        if type=='RGB':
            rainbow_dict[position+2+min]=[ str(int(round(value, 0))) for value in RGB[:] ]
        elif type=='hex':
            rainbow_dict[position+2+min]=('#%02x%02x%02x' % tuple([ int(round(value, 0)) for value in RGB[:] ]))
        position+=1
    return rainbow_dict



def Intersect(dataf1,dataf2,output='/home/vincent/Desktop/bedtools_intersect',name='name',score=0,keep='all',r=0,join='loj',option=''):
    intersect_output=output+'.temp.intersect.bed'
    dataf_1_output=output+'.temp.1.bed'
    dataf_2_output=output+'.temp.2.bed'
    if score == 0:
        dataf1['score']=0
        dataf2['score']=0
        score='score'
    dataf1=dataf1[['chr','start','end',name,score,'strand']]
    dataf2=dataf2[['chr','start','end',name,score,'strand']]

    dataf1['start']=dataf1['start'].map(int)
    dataf1['end']=dataf1['end'].map(int)
    dataf2['start']=dataf2['start'].map(int)
    dataf2['end']=dataf2['end'].map(int)
    dataf1.to_csv(path_or_buf=dataf_1_output,
                          index=False, sep='\t', header=False)
    dataf2.to_csv(path_or_buf=dataf_2_output,
                          index=False, sep='\t', header=False)
    command= "bedtools intersect \
    -a %s \
    -b %s \
    -s " \
    %(dataf_1_output, dataf_2_output)
    command+="-"+join+" "
    if r > 0 :
        command+= "-f %s -r " \
        %(str(r))
    if option != '':
        command+= option
    command+="| sed 's/\t\.\t/\tx\t/g' \
    >  %s" \
    %(intersect_output)
    os.system(command)
    Intersect_dataf=pd.read_csv(filepath_or_buffer=intersect_output,
                       index_col=False, sep='\t', header=None,
                         names=['chr','start','end',name,score,'strand',
                                'dataf2_file','start_2','end_2',name+'_2','overlap','strand_2'])
    del Intersect_dataf['dataf2_file'], Intersect_dataf['strand_2']
    if keep=='name_only':
        Intersect_dataf=Intersect_dataf[[name,name+'_2']]
    command= "rm %s &&\
    rm %s && \
    rm %s" \
    %(dataf_1_output, dataf_2_output, intersect_output)
    os.system(command)
    return Intersect_dataf