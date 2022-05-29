#!/usr/bin/env python


import functools
import pandas as pd





import pydna
import Bio
import Bio.SeqFeature
import Bio.SeqRecord
import Bio.SeqIO
import copy
#%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import pandas as pd
import datetime
import sys
import math
import pathlib
import textwrap as _textwrap
from pydna._pretty import pretty_str as _pretty_str
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord
import csv


''' This part of the design module is used for USER cloning

-----
HELPER FUNCTIONS: 
-----------------

'''


def CAS9_cutting(gRNA_record, background_record):
    """ Simulates cutting by CAS9 given a gRNA
    
    Parameters
    ----------
    gRNA_sequence: string
    
    background_record: pydna.dseqrecord (pydna.assembly)
    
    up_surrounding: int (number of bp)

    dw_surrounding: int (number of bp)

    Returns
    -------
    
    Two dseqrecords

    Examples
    --------
    """
    gRNA_sequence = gRNA_record.seq.watson.upper()
    
    background_sequence = background_record.seq.upper()
    
    gRNA_strand = 1
    if background_sequence.find(gRNA_sequence) == -1:
    #    print('not on +1')
        gRNA_strand = -1
        gRNA_sequence = reverse_complement(gRNA_sequence)
        if background_sequence.find(gRNA_sequence) ==-1:
            print("not on -1, CAN'T FIND THE CUT SITE IN YOUR SEQUENCE")
    #print(gRNA_strand)
    #print(gRNA_sequence)
    
    if gRNA_strand == 1:
        cut_pos_rel_start = 17
    else:
        cut_pos_rel_start = 3
    
    gRNA_start = background_sequence.find(gRNA_sequence)
    
    cut_site = gRNA_start + cut_pos_rel_start
 
    up = background_record[0 : cut_site]
    dw = background_record[cut_site : len(background_sequence)]
    
    up.name = 'UP' + '_' + gRNA_record.name + '_' + background_record.name
    dw.name = 'DW' + '_' + gRNA_record.name + '_' + background_record.name
    
    up_feature = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(0, len(up)), type="misc_feature", strand=+1)
    up_feature.qualifiers["label"] = up.name
    up.features.append(up_feature)
    
    dw_feature = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(0, len(dw)), type="misc_feature", strand=+1)
    dw_feature.qualifiers["label"] = dw.name
    dw.features.append(dw_feature)
    
    up = pydna.dseqrecord.Dseqrecord(up,)
    dw = pydna.dseqrecord.Dseqrecord(dw)
    
    # UPS more than one cut site?
    if dw.seq.find(gRNA_sequence) != -1:
        print("OBS", gRNA_sequence, "cuts more than one location!")
        
    return(up,dw)

def CRIPSIR_knockout(gRNA_record, insertion_site, repair_DNA):
    """Cuts the insertion site with CAS9_cutting and assebmle knockout"""
    #Create fragments after CAS9 cut
    IS_UP, IS_DW = CAS9_cutting(gRNA_record, insertion_site)
    
    #create list of parts and assemble to knockout sequence
    assmeble_parts = [IS_UP, repair_DNA, IS_DW]
    assembled_knockout = pydna.assembly.Assembly(assmeble_parts, limit=200).assemble_linear()[0]
    
    return assembled_knockout


def extract_gRNAs(template, name):
    gRNAs = []
    for feature in template.features:
        if name in feature.qualifiers.get('name',""):
            
            gRNA = template[feature.location.start:feature.location.end]
            gRNA.name = feature.qualifiers.get('name',"")
            gRNAs.append(gRNA)
    
    return(gRNAs)




def remove_features_with_negative_loc(record):
    for i, feature in enumerate(record.features):
        if feature.location.start < 0 or feature.location.start < 0:
            del record.features[i]
            print(feature.qualifiers.get("label","") + " deleted")


def casembler2(bg_strain, site_names=None, gRNAs=None, parts=None, assembly_limits=None, assembly_names=None, verbose=False, to_benchling=False):
    

    """
    Simulate in vivo assembly and integration
    
    Parameters
    ----------
    bg_strain      
    site_names      list of names                      e.g. [X-3, XI-3]    
    gRNAs           list of 20 bp seqrecords           e.g. [ATF1_gRNA, CroCPR_gRNA]
    parts           list of list of parts              e.g. [[ATF1_repair_template],[CPR_repair_template]]
    assembly_limits list of numbers of bp              e.g. [200,400]
    assembly_names  list of names of DNA post assembly e.g. ["X_3_tADH1_P2_pPGK1", "XI_3_UP_DW"]
    verbose         write DNA                          e.g. False
    to_benchling    upload DNA                         e.g. False
    
    Returns
    -------
    
    One dseqrecord

    Examples
    --------
    """
    

    #initialisation
    list_of_list = [site_names, gRNAs, parts, assembly_limits, assembly_names]
    not_empty =[]
    
    #Make sure part is dseqrecord
    if list_of_list[2] != None:
        for part in list_of_list[2]:
            part = pydna.dseqrecord.Dseqrecord.from_SeqRecord(part)
        

    #Iterate through list to differentate None/list values
    for list_no in range(0,len(list_of_list)):
        if list_of_list[list_no] != None:
            not_empty.append(list_of_list[list_no])
        else:
            list_of_list[list_no] = []

    #Check list lengths match and assign variables. 
    if len({len(lists) for lists in not_empty}) != 1:
        print("The provided lists's lengths do not match.")
        for value in not_empty:
            print(str(len(value)))
 

    if site_names is None:
        site_names = []
    else:
        site_names = site_names
    
    if gRNAs is None:
        gRNAs = []
    else:
        gRNAs = gRNAs
          
    if parts is None:
        parts = []
    else:
        parts = parts
            
    if assembly_limits is None:
        assembly_limits = []
    else:
        assembly_limits = assembly_limits
        
    if assembly_names is None:
        assembly_names = []
    else:
        assembly_names = assembly_names
        

    
    
    assemblies = []
    for int_no in range(0,len(gRNAs)):  
        
        gRNA = gRNAs[int_no]

        site_name = site_names[int_no]
        
        site_extract = extract_sites([site_name], [bg_strain], [site_name])[0]
        
        UP, DW = CAS9_cutting(gRNA, site_extract)
        
        fragments = [UP, pydna.dseqrecord.Dseqrecord.from_SeqRecord(parts[int_no]), DW]
        
        assembly = pydna.assembly.Assembly(fragments).assemble_linear()[0]
        
        # sometimes pydna.assembly.Assembly distorts the start, end location of features to become negative which produces and error when printing. This function is created as a workaround. 
        # CPR assembly gives DW_XI_3 annotation called "DW_XI_3" a negative start location.
        # A quick workaround is to remove featuere
        remove_features_with_negative_loc(assembly)
        
        assembly.name = assembly_names[int_no]
        
        assembly_feat = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(0, len(assembly), strand=1), type="misc_feature")
        assembly_feat.qualifiers['name']= site_names[int_no]
        assembly_feat.qualifiers['label']= site_names[int_no]
        assembly.features.append(assembly_feat)
        
        if verbose:
            DNAs = [UP] + [assembly] + [DW]
            for DNA in DNAs:
                DNA.write("./" + DNA.name + ".gb") # "../data/processed/"
        
        if to_benchling:
            to_benchling(assembly, "to_benchling")
        
        assemblies.append(assembly)
        
    return functools.reduce(lambda x, y: x+y, assemblies)



def UPandDW(strain, isite_name):
    """
    take strain and grna
    output UP and DW seqrecord
    """ 
    
    # load lookup table
    gRNAtable = pd.read_csv("../data/raw/gRNAtable.csv", index_col="name")
    
    chromosome_no      = gRNAtable.loc[isite_name,'chromosome']
    
    # load chromosome
    PathToChromosomeSeq = "../data/raw/" + strain + "/" + str(chromosome_no).zfill(2) + ".fa"
    ChromosomeSeq      = Bio.SeqIO.read(PathToChromosomeSeq, "fasta").seq
    
    
    # define homology region location with respect to gRNA sequence
    # f_hom is the length of the UP homology
    # e_hom is the length of the DW homology
    # f_dist is distance from end of UP homology to the first base in isite_sequence
    # e_dist is distance from end of the isite_sequence to DW homology
    f_dist = gRNAtable.loc[isite_name,'f_dist']
    e_dist = gRNAtable.loc[isite_name,'e_dist']
    f_hom   = gRNAtable.loc[isite_name,'f_hom']
    e_hom   = gRNAtable.loc[isite_name,'e_hom']
    
    
    isite_sequence = gRNAtable.loc[isite_name,'sequence']
    isite_sequence = Bio.Seq.Seq(isite_sequence)
    
    # Determine gRNA sequence strand
    gRNA_strand = 1
    if ChromosomeSeq.find(isite_sequence)==-1:
        print('not on +1')
        gRNA_strand = -1
        isite_sequence = isite_sequence.reverse_complement()
        if ChromosomeSeq.find(isite_sequence)==-1:
            print('not on -1')
            print("CAN'T FIND THE CUT SITE IN YOUR SEQUENCE")
    
    # Locate UP and DW
    StartIndex = ChromosomeSeq.find(isite_sequence)
    EndIndex   = StartIndex 
    
    
    
    UPseq = ChromosomeSeq[StartIndex + f_dist - f_hom : StartIndex  + f_dist]
    DWseq = ChromosomeSeq[EndIndex + e_dist          : EndIndex + e_dist+ e_hom]

    UPrec = Bio.SeqRecord.SeqRecord(UPseq, name=isite_name + "UP")
    DWrec = Bio.SeqRecord.SeqRecord(DWseq, name=isite_name + "DW")
    
    # Annotate 
    UP_feature = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(0, len(UPseq)), type="misc_feature", strand=+1)
    UP_feature.qualifiers["label"] = UPrec.name
    UPrec.features.append(UP_feature)
    
    DW_feature = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(0, len(DWseq)), type="misc_feature", strand=+1)
    DW_feature.qualifiers["label"] = DWrec.name
    DWrec.features.append(DW_feature)
    
    return(([UPrec], [DWrec]))


def extract_template_amplification_sites(templates, names, terminator):
    template_amplification_sites = []
    for name, template in zip(names, templates):
#        print(name)
        for feature in template.features:
            if feature.qualifiers['name'] == name:
                CDS_strand = feature.strand
                if CDS_strand == 1:
                    start = feature.location.start
                else:
                    end = feature.location.end
            if feature.qualifiers['name'].endswith(terminator):
                terminator_strand = feature.strand
                if terminator_strand == 1:
                    end = feature.location.end
                else:
                    start = feature.location.start
        
        template_amplification_site = template[start:end]
        
        template_amplification_site.annotations['template_location'] = template.annotations['location']

#        if CDS_strand == terminator_strand:
#            print("CDS and term on same strand")
#        else: 
#            print("OBS! CDS and term not on same strand")
        
        template_amplification_sites.append(template_amplification_site)
    return template_amplification_sites

def extract_template_amplification_sites1(templates, names, terminator):
    template_amplification_sites = []
    for name, template in zip(names, templates):
        for feature in template.features:
            if feature.qualifiers['name'] == name:
                CDS_strand = feature.strand
                if CDS_strand == 1:
                    start = feature.location.start
                else:
                    end = feature.location.end
            if feature.qualifiers['name'].endswith(terminator):
                terminator_strand = feature.strand
                if terminator_strand == 1:
                    end = feature.location.end
                else:
                    start = feature.location.start
        
        template_amplification_site = template[start:end]
        
        if 'batches' in template_amplification_site.annotations.keys():
            template_amplification_site.annotations["batches"].append(template.annotations['batches'][0]) #{'location': template.annotations['batches'][0]['location']}
        else: 
            template_amplification_site.annotations["batches"] = []
            template_amplification_site.annotations["batches"].append(template.annotations['batches'][0]) #{'location': template.annotations['batches'][0]['location']}

#        if CDS_strand == terminator_strand:
#            print("CDS and term on same strand")
#        else: 
#            print("OBS! CDS and term not on same strand")
        
        template_amplification_sites.append(template_amplification_site)
    return template_amplification_sites


def extract_sites1(annotations, templates, names):
    sites = []
    for anno, template, name in zip(annotations, templates, names):
        for feature in template.features:
            if feature.qualifiers['name'] == anno:
                site = template[feature.location.start:feature.location.end]
        site.name = name
        if 'batches' in site.annotations.keys():
            site.annotations["batches"].append(template.annotations['batches'][0])
        else: 
            site.annotations["batches"] = []
            site.annotations["batches"].append(template.annotations['batches'][0])
        sites.append(site) 
    return(sites)



def seq_to_annotation(seqrec_from, seqrec_onto, aType):
    
    ''' Anotate sequences from '''
    seq_from = seqrec_from.seq.watson.upper()
    seq_onto = seqrec_onto.seq.watson.upper()
    
    strand = 1
    match_index = seq_onto.find(seq_from)
    
    if match_index !=-1:
        start = match_index
    #    print("+strand")
        end =  start + len(seq_from)
    #    print("start: ", start)
    #    print("end: ", end)
    else:
    #    print("-strand")
        rev_match_index = reverse_complement(seq_onto).find(seq_from)
        
        if rev_match_index != -1:
            strand = -1
            reclength = len(str(seqrec_onto.seq))
            end = reclength - rev_match_index
            
            start = end - len(seq_from)
    #        print("start: ", start)
    #        print("end: ", end)
        else:
            print("no match! seq:" + seq_from + "\nnot annealing to:" + seq_onto)
    
    feature = Bio.SeqFeature.SeqFeature(Bio.SeqFeature.FeatureLocation(start, end), type=aType, strand=strand)
    feature.qualifiers["label"] = seqrec_from.id
    seqrec_onto.features.append(feature)
    
    def USER_enzyme(amplicon):
        """
        Take pydna.amplicon.Amplicon
        Return Dseqrecord trimmed according to pcr primers
        """
        fw_U_idx = amplicon.forward_primer.seq.find("U")
        #fw_U_idx

        rv_U_idx = amplicon.reverse_primer.seq.find("U")
        #rv_U_idx

        digested_watson = amplicon.seq.watson[fw_U_idx:]
        #digested_watson

        digested_crick = amplicon.seq.crick[rv_U_idx+1:]
        #digested_crick

        digested_pcr = pydna.dseqrecord.Dseqrecord(pydna.dseq.Dseq(watson=digested_watson, crick = digested_crick),
                                 features = amplicon.features, annotations = amplicon.annotations)
        #digested_pcr.seq
        return digested_pcr