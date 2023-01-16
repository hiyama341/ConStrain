#!/usr/bin/env python

# Test fetch_sequences module
from constrain.design.fetch_sequences import retrieve_sequences_from_ncbi, read_fasta_files, retrieve_sequences_from_PDB, fetch_Promoter


from Bio import SeqIO

def test_retrieve_sequences_from_ncbi():
    # random acc_number
    acc_numbers = ['Q05001']
    # call the function: 
    assert retrieve_sequences_from_ncbi(acc_numbers, '../tests/files_for_testing/test_fetch.fasta') == None

    # check if it is the correct sequences
    Q05001_seq = SeqIO.read('../ConStrain/tests/files_for_testing/test_fetch.fasta', 'fasta')
    assert str(Q05001_seq.seq[:10]) == 'MDSSSEKLSP'
    assert Q05001_seq.id == 'sp|Q05001.1|NCPR_CATRO'
 
def test_read_fasta_files(): 
    sequences = read_fasta_files('../ConStrain/tests/files_for_testing/test_fetch.fasta')
    assert sequences[0].seq[:10] == 'MDSSSEKLSP'
    assert sequences[0].id == 'sp|Q05001.1|NCPR_CATRO'

def retrieve_sequences_from_PDB(): 
    acc_numbers = ['Q1PQK4']
    sequences = retrieve_sequences_from_PDB(acc_numbers)
    assert str(sequences[0][0].seq[:10]) == 'MQSTTSVKLS'
    assert str(sequences[0][0].id) == 'sp|A0A2U1LIM9|NCPR1_ARTAN'

def test_fetch_Promoter(): 
    cyc1 = fetch_Promoter('CYC1')
    assert cyc1[:20] == 'GAGGCACCAGCGTCAGCATT'
