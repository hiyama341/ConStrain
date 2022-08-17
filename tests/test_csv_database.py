#!/usr/bin/env python

# Test csv_database
import pandas as pd
import numpy as np
import pytest
import math

# Import the csv_database modules
from constrain.lims.csv_database import *

def test_get_unique_id(): 
    unique_id = get_unique_id(path = '../ConStrain/tests/files_for_testing/csv_database_tests')
    assert unique_id == 10041

def test_get_database(): 
    ds_dna_box = get_database('ds_dna_box', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')
    assert ds_dna_box.iloc[0]['ID'] == 10011.0


def test_get_plate(): 
    plasmid_plates = get_database('plasmid_plates', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')
    box1 = get_plate(0,plasmid_plates )
    assert box1.iloc[0]['ID'] == 10000
    

def test_get_box(): 
    ds_dna_box = get_database('ds_dna_box', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')
    box1 = get_box(0,ds_dna_box )
    assert box1.iloc[0]['name'] == 'UP_XI-2'


def test_add_unique_ids():
    test_templates = []
    for seq_record in SeqIO.parse('../ConStrain/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    
    add_unique_ids(test_templates, path= '../ConStrain/tests/files_for_testing/csv_database_tests')

    assert test_templates[0].id == 10041.0
    assert test_templates[1].id == 10042.0
    assert test_templates[2].id == 10043.0
    assert test_templates[3].id == 10044.0


def test_add_annotations(): 
    test_templates = []
    for seq_record in SeqIO.parse('../ConStrain/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    test_templates

    add_annotations(test_templates, concentration = 100)
    assert test_templates[0].annotations['concentration'] == 100


def test_get_dna_from_plate_name():
    my_dna = get_dna_from_plate_name('pRS416.gb','plasmid_plates',database_path =  '../ConStrain/tests/files_for_testing/csv_database_tests/')
    assert len(my_dna.seq) == 4898

def test_get_dna_from_box_name():
    my_dna = get_dna_from_box_name('AanCPR_tCYC1','ds_dna_box', database_path= '../ConStrain/tests/files_for_testing/csv_database_tests/')
    assert len(my_dna.seq) == 2298




def test_change_row(): 
    #databse
    plasmid_plates = get_database('plasmid_plates', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')

    # new insert
    test_templates = []
    for seq_record in SeqIO.parse('../ConStrain/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    biopython_object = test_templates[0]
    biopython_object.id = 999999
    biopython_object.annotations['concentration'] = ''
    biopython_object.annotations['comments'] = ''
    biopython_object.annotations['reference'] = ''
    biopython_object.annotations['volume'] = ''
    change_row(0, plasmid_plates,  biopython_object)


    assert plasmid_plates.iloc[0]['ID'] == 999999



def test_delete_row_df(): 

    #databse
    plasmid_plates = get_database('plasmid_plates', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')
    delete_row_df(0, plasmid_plates)
    
    assert str(plasmid_plates.iloc[0]['ID']) == 'nan'




def test_add_sequences_to_dataframe(): 
    test_templates = []

    for seq_record in SeqIO.parse('../ConStrain/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'): 
        test_templates.append(seq_record) 
    add_annotations(test_templates)
    add_unique_ids(test_templates, path= '../ConStrain/tests/files_for_testing/csv_database_tests')

    plasmid_plates = get_database('plasmid_plates', path= '../ConStrain/tests/files_for_testing/csv_database_tests/')

    add_sequences_to_dataframe(test_templates,plasmid_plates, index= 20 )

    assert str(plasmid_plates.iloc[20]['name']) == 'pCCW12'
    assert str(plasmid_plates.iloc[21]['name']) == 'pTPI1'
    assert str(plasmid_plates.iloc[22]['name']) == 'pCYC1'
    assert str(plasmid_plates.iloc[23]['name']) == 'pENO2'





def test_update_database(): 
    # making a test to test if pandas writes a csv correctly
    pass