#!/usr/bin/env python


import functools
import pandas as pd
import pydna
import Bio
import Bio.SeqFeature
import Bio.SeqRecord
import Bio.SeqIO
import copy

''' This part of the design module is used for USER cloning

-----
HELPER FUNCTIONS: 
-----------------
USER_enzyme
CAS9_cutting
CRIPSIR_knockout
extract_gRNAs
remove_features_with_negative_loc
casembler2
extract_template_amplification_sites1
UPandDW
extract_sites1
seq_to_annotation

'''


def CAS9_cutting(gRNA_record, background_record):
    from Bio.Seq import Seq

    
    """ Simulates cutting by CAS9 given a gRNA
    
    Parameters
    ----------
    gRNA_sequence: pydna.dseqrecord
    
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

        gRNA_strand = -1
        gRNA_sequence = Seq(gRNA_sequence)

        gRNA_sequence = (gRNA_sequence).reverse_complement()
        if background_sequence.find(gRNA_sequence) ==-1:
            print("not on -1, CAN'T FIND THE CUT SITE IN YOUR SEQUENCE")

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
    from pydna.assembly import Assembly
    """Cuts the insertion site with CAS9_cutting and assebmle knockout"""
    #Create fragments after CAS9 cut
    IS_UP, IS_DW = CAS9_cutting(gRNA_record, insertion_site)
    
    #create list of parts and assemble to knockout sequence
    assmeble_parts = IS_UP, repair_DNA, IS_DW
    assembled_knockout = Assembly(assmeble_parts).assemble_linear()[0]
    
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
    ''' Removes a features if it is negative'''
    for i in range(len(record.features)):
        if record.features[i].location.start < 0 or record.features[i].location.end < 0:
            del record.features[i]

""" ### the old function for removing negative features
def remove_features_with_negative_loc1(record):
    for i, feature in enumerate(record.features):
        if feature.location.start < 0 or feature.location.start < 0:
            del record.features[i]
            print(feature.qualifiers.get("label","") + " deleted")
 """



def extract_template_amplification_sites(templates, names, terminator):
    '''
    This function takes in 3 parameters: 
    
    '''
    template_amplification_sites = []
    for name, template in zip(names, templates):
        for feature in template.features:
            if feature.qualifiers['name'] == name:
                
                CDS_strand = feature.strand
                if CDS_strand == 1:
                    start = feature.location.start
                    print(start)
                else:
                    end = feature.location.end
                    print(end)

            if feature.qualifiers['name'].endswith(terminator):
                terminator_strand = feature.strand
                if terminator_strand == 1:
                    end = feature.location.end
                else:
                    start = feature.location.start
                    print(start)

        
        template_amplification_site = template[start:end]
        
        if 'batches' in template_amplification_site.annotations.keys():
            template_amplification_site.annotations["batches"].append(template.annotations['batches'][0]) #{'location': template.annotations['batches'][0]['location']}
        else: 
            template_amplification_site.annotations["batches"] = []
            template_amplification_site.annotations["batches"].append(template.annotations) #{'location': template.annotations['batches'][0]['location']}

        template_amplification_sites.append(template_amplification_site)
    return template_amplification_sites

def extract_sites(annotations, templates, names):
    ''' This function extracts the sequences from annotated sequences based on they names
    
    '''
    sites = []
    for anno, template, name in zip(annotations, templates, names):
        for feature in template.features:
            site = template[feature.location.start:feature.location.end]
            if feature.qualifiers['name'] == anno:
                site.name = name

            # If there is an batch anotation we can save it 
            if 'batches' in site.annotations.keys():
                site.annotations["batches"].append(template.annotations['batches'][0])
            else: 
                site.annotations["batches"] = []
                site.annotations["batches"].append(template.annotations)

            sites.append(site) 
    return(sites)



def seq_to_annotation(seqrec_from, seqrec_onto, aType):
    from Bio.Seq import Seq

    ''' Anotate an amplicon object from another amplicon object'''
    seq_from = seqrec_from.seq.watson.upper()
    seq_onto = seqrec_onto.seq.watson.upper()
    
    strand = 1
    match_index = seq_onto.find(seq_from)
    print(match_index)
    
    
    # if there is match
    if match_index !=-1:
        start = match_index
        end =  start + len(seq_from)
    else:
        # If no match we look at the reverse complement
        seq_onto = Seq(seq_onto)
        rev_match_index = seq_onto.reverse_complement().find(seq_from)
        
        # if we get a match here
        if rev_match_index != -1:
            strand = -1
            reclength = len(str(seqrec_onto.seq))
            end = reclength - rev_match_index
            
            start = end - len(seq_from)
            print("start: ", start)
            print("end: ", end)
        else:
            print("no match! seq:" + str(seqrec_from.name) + "\nnot annealing to:" + str(seqrec_onto.name))
    
    # add the feature to the amplicon
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

    digested_watson = amplicon.seq.watson[fw_U_idx+1:]
    #digested_watson

    digested_crick = amplicon.seq.crick[rv_U_idx+1:]
    #digested_crick

    digested_pcr = pydna.dseqrecord.Dseqrecord(pydna.dseq.Dseq(watson=digested_watson, crick = digested_crick),
                                features = amplicon.features, annotations = amplicon.annotations)
    #digested_pcr.seq
    return digested_pcr