#!/usr/bin/env python3

import sys, os
import re
import math
import glob
import random
import pickle
import tmhmm
import numpy as np
import pandas as pd
import seaborn as sns
import ruptures as rpt
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import poisson, binned_statistic, norm, lognorm

###########
# Define global tmhmm model path 
tmhmm_dir = os.path.dirname(tmhmm.__file__)
model_path = os.path.join(tmhmm_dir, 'TMHMM2.0.model')

###########

# Save and load dictionaries with pickle
def SavePickle(obj, File):
    with open(File, "wb") as myFile:
        pickle.dump(obj, myFile)
        myFile.close()

def LoadPickle(File):
    with open(File, "rb") as myFile:
        obj = pickle.load(myFile)
        myFile.close()
        return obj

###########

def alt_startcodon(ide, nt_seqs):
    
    codons = [nt_seqs[ide][i:i+3] for i in range(0, len(nt_seqs[ide]), 3)]
    return codons

def alt_startcodon_seq(seq, starts=['ATG', 'TTG', 'GTG']):
    return [1 if seq[i:i+3] in starts else 0 for i in range(0, len(seq), 3)]

def check_frame(ins, start, strand='pos'):
    f = (start-ins)%3
    if strand=='pos':
        return {0:1, 2:2, 1:3}[f]
    else:
        return {0:3, 2:2, 1:1}[f]
    
def predict_tm(aa_sequence, posterior=False):
    return tmhmm.predict(aa_sequence, header='', model_or_filelike=model_path, compute_posterior=posterior)

def tm2nttm(aa_sequence, posterior=False):
    return ''.join([i*3 for i in predict_tm(aa_sequence, posterior=posterior)])

def plot_piles(st, en, strand, all_piles):
    for experiment, piles in all_piles.items():
        pile = piles[strand]
        expression = np.log2(pile[st-1:en]+1)
        plt.plot(range(st-1, en), expression, linewidth=3, label=experiment, alpha=0.5)

def return_sample(dic, select_all, select_any=None, exclude=None):
    if select_any:
        _p = {k:v for k,v in dic.items() if all([i in k for i in select_all]) and (any([i in k for i in select_any]))}
    else:
        _p = {k:v for k,v in dic.items() if all([i in k for i in select_all])}
        
    if exclude:
        _p = {k:v for k,v in _p.items() if all([i not in k for i in exclude])}
    return _p

def minmaxstandarization(value, distribution=None):
    """
    Apply a min max standarization to the value.
    Value can be another distribution
    """
    value = np.array([float(i) for i in value])
    if hasattr (distribution, "__len__"):
        return (np.array(value)-min(distribution))/(max(distribution)-min(distribution))
    else:
        return (np.array(value)-min(value))/(max(value)-min(value))    

def windows(array, window_size=False, bins=False, metric='mean', circular=True, overlapping=True):
    """
    Computes the windows for array with window size.
    If circular, the array is considered to be circular

    If overlapping==False, this function produces a bin profile

    bins= represent the number of portions of splitted data (if 100, the profile is separated in 100
    non-overlapping subarray.
    window_size= how many bases each bin has to cover
    """
    accepted_metrics = set(['sum', 'std', 'mean', 'median', 'count', 'L', 'min', 'max', 'I', 'dens', 'R', 'RI', 'CPM', 'RPKM'])
    if type(array)!=np.ndarray:
        array = np.array(array)
    array = array.astype('float')

    if not window_size and not bins:
        sys.exit('At least one of the following arguments is required: window_size or bin')
    if metric not in accepted_metrics:
        sys.exit(metric+' not accepted. Please provide one of the following:\n- '+'\n- '.join(accepted_metrics))
    if metric=='L':
        metric='count'

    if metric in ['dens', 'RI', 'CPM', 'RPKM']:
        original_array = array.copy() # Required to work with recursion

    if circular and overlapping:
        atail = array[-int(window_size/2.0):]
        ahead = array[:int((window_size)/2.0)]
        array = np.concatenate((atail,array,ahead))

    if metric=='I' or metric=='dens':
        nz_array = array.copy()
        if overlapping:
            nz_array[nz_array==0.0]='nan'
        else:
            nz_array[nz_array>0.0]=1.0  # Binary array

    if overlapping:
        if metric=='I':
            roll = pd.Series(np.array(nz_array)).rolling(window_size, center=True)
        else:
            if metric not in ['RI', 'RPKM', 'CPM']:
                roll = pd.Series(np.array(array)).rolling(window_size, center=True)

        if metric in ['sum', 'R']:
            roll = roll.sum()
        elif metric=='std':
            roll = roll.std()
        elif metric=='mean':
            roll = roll.mean()
        elif metric=='median':
            roll = roll.median()
        elif metric in 'count':
            roll = roll.count()[len(ahead):len(array)-len(atail)]   # Required to correct the tails
        elif metric=='min':
            roll = roll.min()
        elif metric=='max':
            roll = roll.max()
        elif metric=='I':
            # This is a nonzero_count function as we change 0-->NaN in array variable, check up. 
            roll = roll.count()[len(ahead):len(array)-len(atail)]   # Required to correct the tails
        elif metric=='dens':
            _is  = windows(original_array, window_size=window_size, bins=bins, metric='I', circular=circular, overlapping=overlapping)
            roll = _is/window_size
        elif metric=='RI':
            _is  = windows(original_array, window_size=window_size, bins=bins, metric='I', circular=circular, overlapping=overlapping)
            _rs  = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = _rs*replace_zeroes(_is)
        elif metric=='CPM':
            roll = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = (roll/sum(original_array))*1e6
        elif metric=='RPKM':
            roll = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = roll/((window_size/1000.0)*(sum(original_array)/1e6))
        if metric not in ['dens', 'RI', 'RPKM', 'CPM']:
            roll.dropna(inplace=True)    # Remove tails

    else:
        if type(bins)!=int:
            bins = int(float(len(array))/window_size)

        if bins and not window_size:
            window_size = int(len(array)/bins)

        if metric=='R':
            roll = binned_statistic(np.arange(len(array)), array, bins=bins, statistic='sum')[0]
        elif metric=='I':
            roll = binned_statistic(np.arange(len(nz_array)), nz_array, bins=bins, statistic='sum')[0]
        elif metric=='dens':
            roll = binned_statistic(np.arange(len(nz_array)), nz_array, bins=bins, statistic='sum')[0]
            roll/=window_size
        elif metric=='RI':
            _is  = windows(original_array, window_size=window_size, bins=bins, metric='I', circular=circular, overlapping=overlapping)
            _rs  = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = _rs*_is
        elif metric=='CPM':
            roll = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = (roll/sum(original_array))*1e6
        elif metric=='RPKM':
            roll = windows(original_array, window_size=window_size, bins=bins, metric='R', circular=circular, overlapping=overlapping)
            roll = roll/((window_size/1000.0)*(sum(original_array)/1e6))
        else:
            roll = binned_statistic(np.arange(len(array)), array, bins=bins, statistic=metric)[0]

    return np.array(roll)

def sample_reads(coords, areads, frame=1, ignore_strand=True):
    """ INPUT = [dataset, genename, frame] """
    start, end, strand = coords
    gene_reads = areads[start-1:end]
    if strand=='-' or strand==0:
        gene_reads = gene_reads[::-1]
    if frame!=0:
        gene_reads = gene_reads[frame-1::3]
    if ignore_strand and (strand=='-' or strand==0):
        gene_reads = gene_reads[::-1]
    return gene_reads

def get_f_metric(reads, annotation, thr=0):
    ''' Calculates total ribosome counts by frame and in total '''
    start, end, strand = annotation
    if strand=='+':
        frame1 = [i for i in range(start-1,end, 3)]
        frame2 = [i for i in range(start  ,end, 3)]
        frame3 = [i for i in range(start+1,end, 3)]
        strd = 'pos'
    else:
        frame3 = [i for i in range(start-1,end, 3)]
        frame2 = [i for i in range(start  ,end, 3)]
        frame1 = [i for i in range(start+1,end, 3)]
        strd = 'neg'
    
    rs = {}
    c = 0
    for fr, frame_id in zip([frame1, frame2, frame3], [1,2,3]):
        subreads = reads[fr]
        rs[c]   = [frame_id, 'R', sum(subreads)]
        rs[c+1] = [frame_id, 'I', len(subreads[subreads>thr])]
        rs[c+2] = [frame_id, 'coverage', len(subreads[subreads>thr])/((annotation[1]-annotation[0])/3)]
        c += 3
    
    rs = pd.DataFrame.from_dict(rs, orient='index', columns=['frame', 'metric', 'value'])
    return rs


def psplot(coords,
           extra=50, 
           profile=None, 
           tracks=[1,1,1,1,1,1,1], log=True, 
           figsize=[5,1], _annmode='ps', _repr='scatter', start_codons=['ATG', 'GTG', 'TTG'], **kwargs):
    """
    coords = [st, en, strand], 
    tracks = [omics, insertions, tm, annotation]

    kwargs:
        nt_seq = 'ATGCT...TATTAA'
        aa_seq = 'MK...RCT'
        piles = {'sampleA':{'+':np.array}}
        selection = ['Cm_B_15_1', 'Cm_B_5_2', etc.]
        query = ['Cm'] --> includes all Cm samples // [['Cm'], ['10', 15']] --> sample has to be 'Cm' and any of '10' and '15'
        annotation = orfs dataframe
        ide = ide to plot or to base the annotation reference overlaps

        # Booleans
        alt_start = 1 # Plot alternative start codons 
        tm = transmembrane estimation
    """
    # Define coords
    st, en, strand = coords
    xaxis = np.array(range(st-extra, en+extra))
    if 'colord' in kwargs:
        colord = kwargs['colord']
    else:
        colord = {1:'darkorange', 2:'dodgerblue', 3:'purple'}
    
    # Checkers and setting up the plotting environment
    axtracks = {}
    c = 0
    if 'piles' in kwargs and tracks[0]==1:
        if _annmode=='ps':
            axtracks['piles'] = c
            c+=1
        else:
            axtracks['piles+'] = c
            c+=1
            axtracks['piles-'] = c
            c+=1
    if 'peptides' in kwargs and tracks[1]==1 and kwargs['peptides']:
        axtracks['peptides'] = c
        c+=1
    if tracks[2]==1 and profile:
        if type(profile)!=dict:
            profile = {'sample':profile} # Assumed an array is passed
        if 'query' in kwargs:
            if type(kwargs['query'])==str:
                profile = return_sample(profile, select_all=[kwargs['query']])
            elif type(kwargs['query'][0])==str:
                profile = return_sample(profile, select_all=kwargs['query'])
            elif type(kwargs['query'][0])==list:
                profile = return_sample(profile, select_all=kwargs['query'][0], select_any=kwargs['query'][1])
            else:
                sys.exit('provide accepted query')
        for k in sorted(list(profile.keys())):
            if _repr!='comb' or k[0]!='T':
                axtracks[k] = c
                c+=1
    if 'aa_seq' in kwargs and tracks[3]==1:
        axtracks['tm'] = c
        c+=1
    if 'wdw_size' in kwargs and tracks[4]==1:
        axtracks['cpc'] = c
        c+=1
    if 'domains' in kwargs and tracks[5]==1:
        for k in sorted(list(kwargs['domains'].keys())):
            axtracks['domains '+k] = c
            c+=1
    if 'annotation' in kwargs and tracks[6]==1:
        axtracks['annotation'] = c
        c+=1

    if len(axtracks)==1:
        f, ax = plt.subplots(c, 1, figsize=(figsize[0], figsize[1]*2), sharex=True)   # Define the figure
        axtracks = {k:ax for k, v in axtracks.items()}
    else:
        f, ax = plt.subplots(c, 1, figsize=(figsize[0], figsize[1]*c), sharex=True)   # Define the figure
        axtracks = {k:ax[v] for k, v in axtracks.items()}
    
    if strand=='+':
        ref = 1*st
    else:
        ref = 1*en

    # Plot tracks
    maximum = 0
    for k in axtracks.keys():
        if k.startswith('piles'):
            if len(k)==5:
                pilestrand = strand
            else:
                pilestrand = k[-1]
            palette={k:_color for k, _color in zip(kwargs['piles'].keys(), sns.color_palette('mako', len(kwargs['piles'])))}
            for experiment_ide, pile in kwargs['piles'].items():
                expr_ys = np.log2(pile[pilestrand][xaxis-1])
                axtracks[k].plot(xaxis, expr_ys, label=experiment_ide, alpha=0.5, color=palette[experiment_ide])
                if max(expr_ys)>maximum:
                    maximum = max(expr_ys)
            axtracks[k].set_ylabel('log2(CPM)\nstrand {}'.format(pilestrand))
    for k in axtracks.keys():
        if k.startswith('piles'):
            axtracks[k].set_ylim(0, maximum+0.5)
    maximum=0
           
    if 'peptides' in axtracks and len(kwargs['peptides'])>0:
        c=1
        for pept in sorted(kwargs['peptides']):
            if strand=='-':
                a, b = [ref-(i*3)-3 for i in pept]
            else:
                a, b = [ref+(i*3)-3 for i in pept]
            axtracks['peptides'].plot([a, b],[c/4,c/4], c="seagreen", linewidth=4, alpha=0.5)   # alt color = purple
            c+=1
            if c==9:
                c=1
        axtracks['peptides'].set_ylim(0, 2.25)
        axtracks['peptides'].set_yticks([])
        axtracks['peptides'].set_ylabel('UTPs')
        
    if profile and tracks[2]==1:
        np.seterr(divide = 'ignore')
        maximum = 0
        for sample, reads in profile.items():
            if log:
                subreads = np.log2(reads[st-extra-1:en+extra])
                nxaxis = np.arange(st-extra, en+extra+1)
                nxaxis = nxaxis[subreads>=0]
                subreads = subreads[subreads>=0]
            else:
                print('nolog')
                subreads = np.array(reads[st-extra-1:en+extra])
                nxaxis = np.arange(st-extra, en+extra+1)
                nxaxis = nxaxis[subreads>0]
                subreads = subreads[subreads>0]
            
            # update max
            if len(subreads)>0 and max(subreads)>maximum:
                maximum = max(subreads)
           
            # colord = {1:np.array([sns.color_palette("rocket",3)[1]]), 2:np.array([sns.color_palette("mako",3)[0]]), 3:np.array([sns.color_palette("mako",3)[0]])}
            if _annmode=='ps':
                if strand =='+':
                    colors = [colord[check_frame(i, st, 'pos')] for i in nxaxis]
                else:
                    colors = [colord[check_frame(i, en, 'pos')] for i in nxaxis]
            else:
                colors = colord[sample]
                
            # Set some parameters
            if kwargs.get('s', None):
                s = kwargs['s']
            else:
                s = 20
            if kwargs.get('width', None):
                width = kwargs['width']
            else:
                width = 8
                if _repr=='comb':
                    width = 12
            # Plot 
            if _repr=='comb':
                alpha=0.75
                if sample[0]=='T':
                    sample = sample.replace('T', 'P')
                    barh = 1
                    barbottom = 1
                    axtracks[sample].axhline(y=2, c='dimgrey', linewidth=0.5)
                else:
                    barh = 1
                    barbottom = 0
            else:
                alpha = 1
                barh = 1
                barbottom = 0

            if _repr=='bar':
                axtracks[sample].bar(nxaxis, subreads, color=colors, label='reads', edgecolor="none", width=width)
            else:
                axtracks[sample].scatter(nxaxis, subreads, color=colors, label='reads', marker='v', s=s, alpha=alpha)
                if _repr!='scatter':
                    #sns.rugplot(x=nxaxis, y=subreads, ax=axtracks[sample])        
                    axtracks[sample].bar(nxaxis, np.where(subreads > 0, barh, subreads), color=colors, label='reads', edgecolor="none", width=width, bottom=barbottom)
            
            # Set labels and ticks
            #sample_patch = mpatches.Patch(color='black', label='{} [log2(reads)]'.format(sample))
            #axtracks[sample].legend(handles=[sample_patch], loc='center left', bbox_to_anchor=(1, 0.5))
            if log:
                axtracks[sample].set_ylabel('{}\nlog2(reads)'.format(''.join(sample.split('_'))), rotation=0, labelpad=25)
            else:
                axtracks[sample].set_ylabel('{}\nreads count'.format(''.join(sample.split('_'))))    
        #for sample in profile.keys():
        #    axtracks[sample].set_ylim(-0.2, maximum+0.2)

    if 'tm' in axtracks:
        # subset annotation
        tmd = {'i':1, 'M':3, 'o':5, 'O':5}
        nttm = tm2nttm(kwargs['aa_seq'])
        if strand=='-':
            nttm = nttm[::-1]
        nttm = [tmd.get(i, 0) for i in nttm]

        axtracks['tm'].axhline(y=4, linestyle='--', c='dimgrey', label='outer leaflet')
        axtracks['tm'].axhline(y=2, linestyle='--', c='darkgrey', label='inner leaflet')
        axtracks['tm'].plot(xaxis, ([np.nan]*(extra+1))+nttm+([np.nan]*(extra+1)), color='seagreen', linewidth=4, label='TMHMM', alpha=0.7)
        axtracks['tm'].set_ylim(0, 6)
        axtracks['tm'].set_ylabel('TMHMM')
        axtracks['tm'].set_yticks((1,3,5))
        axtracks['tm'].set_yticklabels(('i', 'M', 'o'))    
        
    if 'cpc' in axtracks and profile:
        # extract in frame reads and compute windows
    
        # Run cpc
        wdw_profile = {k:windows(minmaxstandarization(sample_reads([st, en+kwargs['wdw_size']*3, strand], v, frame=1, ignore_strand=True)), window_size=kwargs['wdw_size'], metric='mean', circular=False, overlapping=True) for k, v in profile.items()}
        signal = np.array(list(wdw_profile.values())).T
        algo = rpt.Pelt(model="rbf").fit(signal)
        result = algo.predict(pen=5)
        result = [i*3+st+kwargs['wdw_size'] for i in result][:-1]  # last result is the closure
        
        if kwargs.get('aa_seq'):
            nttm = tm2nttm(kwargs['aa_seq']).lower()
            if strand=='-':
                nttm = nttm[::-1]
            topo_trans = [i.start(0)+st for i in re.finditer(r"im|mo|om|mi",nttm)]
        
        subxaxis = np.arange(st, en+1+kwargs['wdw_size']*3)[::3][:signal.shape[0]]
        axtracks['cpc'].plot(subxaxis, signal, color='lightgrey')
        axtracks['cpc'].plot(subxaxis, np.mean(signal, axis=1), color=colord[1])
        for i in result:
            axtracks['cpc'].axvline(x=i, linestyle='--', color='grey')
        if kwargs.get('aa_seq'):
            for i in topo_trans:
                axtracks['cpc'].axvline(x=i, linestyle='--', color='mediumseagreen')
        axtracks['cpc'].set_ylim(0, np.max(signal+0.01))
        axtracks['cpc'].set_ylabel('log2(reads)')
        
    if 'domains' in kwargs and tracks[5]==1:        
        for k, doms in kwargs['domains'].items():
            span = set(range(st-extra, en+extra+1))
            subann =  doms[((doms['start']<st) & (doms['end']>en))|
                           ((doms['start'].isin(span)) | ( doms['end'].isin(span)))].copy()
            x = []
            y = []
            for index, row in subann.iterrows():
                x+=[row['start'], row['end']]
                y+=[row['mean_R'], row['mean_R']]
            y = np.log2(np.array(y)+1)
            axtracks['domains '+k].plot(x, y, c='darkcyan')
            axtracks['domains '+k].set_ylabel(k+'\n log2(reads)')
            axtracks['domains '+k].set_ylim(0, np.max(y)+0.3)

            for index, row in subann.iterrows():
                if row['class']=='E':
                    axtracks['domains '+k].axvspan(row['start'], row['end'], facecolor='grey', alpha=0.2)
    
    if 'annotation' in axtracks:
        if _annmode=='ps':
            span = set(range(st-extra, en+extra+1))
            subann =  kwargs['annotation'][(kwargs['annotation']['strand']==strand) & 
                                          # (kwargs['annotation'].index!=kwargs.get('ide', 'region')) &
                                           (((kwargs['annotation']['start']<st) & (kwargs['annotation']['end']>en) & (kwargs['annotation']['alt_ann'].str.startswith('M')))|
                                           ((kwargs['annotation']['start'].isin(span)) | ( kwargs['annotation']['end'].isin(span))))].copy()
            if strand=='+':
                subann['frame'] = [check_frame(a, st, 'pos') for a in subann['start']]
            else:
                subann['frame'] = [check_frame(a, en, 'pos') for a in subann['end']]
            subann['color'] = subann['frame'].map(colord)
            ys = []
            for altname, fr in zip(list(subann['alt_ann']), list(subann['frame'])):
                if str(altname)=='0':
                    ys.append(fr+3)
                elif str(altname)=='-1':
                    ys.append(fr)
                else:
                    ys.append(fr+6)
            subann['y'] = ys
                       
            # Plot
            for i, row in subann.iterrows():
                axtracks['annotation'].plot([row['start'], row['end']], [row['y'], row['y']], c=row['color'], linewidth=4, alpha=0.7)   # alt color = purple
                axtracks['annotation'].set_ylim(0, 10)
                axtracks['annotation'].set_ylabel('Overlaps')
                axtracks['annotation'].set_yticks((2,5,8))
                axtracks['annotation'].set_yticklabels(('Micro (<10 aa)', 'smORFs', 'Annotated CDS')) 
            axtracks['annotation'].axhline(y=3.5, c='lightgrey')
            axtracks['annotation'].axhline(y=6.5, c='lightgrey')

           
    # Adjust general
    for k, subax in axtracks.items():
        subax.set_xlim(min(xaxis), max(xaxis))

        # Plot annotation
        if strand == '+':
            subax.axvline(x=st, c='dimgrey', label='Start')
            subax.axvline(x=en, c='grey', label='Stop')
        else:
            subax.axvline(x=en, c='dimgrey', label='Start')
            subax.axvline(x=st, c='grey', label='Stop')
    
    # General legend
    if _annmode=='ps' and tracks[2]:
        f1_patch = mpatches.Patch(color=colord[1], label='f1')
        f2_patch = mpatches.Patch(color=colord[2], label='f2')
        f3_patch = mpatches.Patch(color=colord[3], label='f3')
        axtracks[list(profile.keys())[0]].legend(handles=[f1_patch, f2_patch, f3_patch], loc='center left', bbox_to_anchor=(1, 0.5))
            
    # Set labels and ticks
    plt.xlim(min(xaxis), max(xaxis))
    plt.xlabel('Genome base position \n {} - {}..{}{}'.format(kwargs.get('ide', 'region'), st, en, ' (complement)' if strand=='-' else ''))
    f.subplots_adjust(hspace = .1)
    f.tight_layout()
    
    f.subplots_adjust(hspace=.0)
    
    if kwargs.get('save', None):
        f.savefig(kwargs['save'])
    return f


######

def metagene_by_frame(zreads, several_annotations, frame=1, metric='R', step=100, gl=816394, filt=1):
    """ extend_mode can be 'perc' or 'base' """
    metaprofile = []
    for annotation in several_annotations:
        start, end, strand = sorted(annotation[:2])+[annotation[-1]]
        gene_reads = zreads[start-1:end]
        if strand=='-' or strand==0:
            gene_reads = gene_reads[::-1]
        if frame!=0:
            gene_reads = gene_reads[frame-1::3]
        gene_reads = np.array([0 if i<filt else i for i in gene_reads])
        wdws = windows(gene_reads, bins=step, metric=metric, circular=False, overlapping=False)
        metaprofile.append(wdws)
    return metaprofile

######

def window_filter(a, n=3, func='max'):
    new_a = list(a)
    for i in range(len(new_a)):
        sub_a = new_a[i:i+3]
        dd = {0:0, 1:0, 2:0}
        ii = sub_a.index(max(sub_a))
        if func=='max':
            dd[ii] = sub_a[ii]  # keep the maximum
        elif func=='sum':
            dd[ii] = sum(sub_a)
        
        if i<len(new_a):
            new_a[i] = dd[0]
        if i+1<len(new_a):
            new_a[i+1] = dd[1]
        if i+2<len(new_a):
            new_a[i+2] = dd[2]
    return new_a

def max_win(d):
    new_reads = window_filter(d.zreads)
    new_profile = {i+1:new_reads[i] for i in range(0, 816394) if new_reads[i]>0}
    v = d.copy('swallow')
    v._update_stats_(new_profile)
    return v

def sum_win(d):
    new_reads = window_filter(d.zreads, func='sum')
    new_profile = {i+1:new_reads[i] for i in range(0, 816394) if new_reads[i]>0}
    v = d.copy('swallow')
    v._update_stats_(new_profile)
    return v

def annotation_overlaps(dic, strand=False, threshold=2):
    positions = []
    for k, v in dic.items():
        if strand and v[-1]==strand:
            positions+=list(range(v[0], v[1]+1))        
        else:
            positions+=list(range(v[0], v[1]+1))                    
    positions = Counter(positions)
    if threshold>0:
        counted = [k for k, v in positions.items() if v>=threshold]
    else:
        counted = {}
        for k, v in positions.items():
            if v in counted:
                counted[v].append(k)
            else:
                counted[v]=[k]
    return counted