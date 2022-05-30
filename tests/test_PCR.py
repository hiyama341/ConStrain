#!/usr/bin/env python

# Test PCR module

from pydna.dseqrecord import Dseqrecord
from pydna.readers import read
from pydna.amplify import pcr
from pydna.primer import Primer
from Bio.SeqRecord import SeqRecord


# Importing the module we are  testing
from constrain.lab.PCR import *



def test_pcr_volumes(): 
    vol_p_reac = 10,
    
    #no_of_reactions = 6,
    #standard_reagents = ["DNA","Buffer, Cutsmart","H20","Enz, USER"],
    # standard_volumes = [1,1,7,1]
    #pcr_reaction_df = pcr_volumes(vol_p_reac,no_of_reactions ,standard_reagents, standard_volumes)
    
    pass
    
def test_det_proc_speed(): 
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    
    # the tests
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    det_proc_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 60
    
    amplicon.annotations['polymerase'] = "Q5 Hot Start"
    det_proc_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

    amplicon.annotations['polymerase'] = "Phusion"
    det_proc_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

def test_det_elon_time():
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    det_proc_speed(amplicon)
    #test
    det_elon_time(amplicon)
    assert amplicon.annotations['elongation_time'] == 4
    
    
    # initialize 2
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    det_proc_speed(amplicon)
    det_elon_time(amplicon)

    #tests
    assert amplicon.annotations['elongation_time'] == 128
    
    
    # initialize 2
    middle = 'a'*3000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    det_proc_speed(amplicon)
    det_elon_time(amplicon)

    #test
    assert amplicon.annotations['elongation_time'] == 188

    

def test_PCR_program(): 
    # no test atm
    pass


def test_grouper():
    elong_times = [60,60, 46, 60, 45, 30, 200, 100]
    elong_times.sort()
    elong_time_max_diff = 10
    groups = dict(enumerate(grouper(elong_times,elong_time_max_diff), 1))

    assert groups == {1: [30], 2: [45, 46], 3: [60, 60, 60], 4: [100], 5: [200]}





def test_det_no_of_thermal_cyclers():
    # no test atm
    pass


def test_pcr_locations():
    # not test - hard to test a df
    pass


def test_nanophotometer_concentrations(): 
    pass
