#!/usr/bin/env python

# Test USER cloning module

# Importing the module we are  testing
from pydna.amplify import pcr
from constrain.design.cloning import *


def test_USER_enzyme(): 
    
    # inititalize
    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'

    gRNA1_pcr_prod = pcr(primerF,primerR, template)

    Digested_USER = USER_enzyme(gRNA1_pcr_prod)
    assert len(Digested_USER.seq.watson) == 417
    assert Digested_USER.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'


def test_remove_features_with_negative_loc():
    from Bio.SeqRecord import SeqRecord

    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'
    rec_vector = pcr(primerF,primerR, template)
    rec_vector.add_feature(-40,60, label =['gRNA'] )
    assert len(rec_vector.features) == 3

    # remove negative featureS
    remove_features_with_negative_loc(rec_vector)
    assert len(rec_vector.features) == 2


def test_CAS9_cutting():
    from pydna.dseqrecord import Dseqrecord

    template = Dseqrecord('TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT')
    gRNA = Dseqrecord('TCTAGATTTTGTAGTGCCCT')
    up, dw = CAS9_cutting(gRNA, template)

    assert len(up)== 125
    assert up.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGC'

    assert len(dw) == 292
    assert dw.seq.watson == 'CCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'



def test_CRIPSIR_knockout():
    # initialize
    from pydna.dseqrecord import Dseqrecord
    insertion_site = Dseqrecord('TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT')
    gRNA = Dseqrecord('TCTAGATTTTGTAGTGCCCT')
    repair_template = Dseqrecord('CGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGT')

    # call the function
    Knock_out = CRIPSIR_knockout(gRNA,insertion_site, repair_template)
    assert Knock_out.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'

def test_extract_gRNAs():
    from Bio.SeqRecord import SeqRecord
    from pydna.amplify import pcr

    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'

    rec_vector = pcr(primerF,primerR, template)

    # Adding a feature
    rec_vector.add_feature(40,60, name ='gRNA' , label = ['gRNA'])
    gRNA = extract_gRNAs(rec_vector, 'gRNA')

    assert gRNA[0].seq.watson == 'ACTCATATTTATACAGAAAC'



def test_extract_template_amplification_sites():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # Initializing
    features1 = [SeqFeature(FeatureLocation(1, 100, strand=1), type='CDS')]
    features2 = [SeqFeature(FeatureLocation(101, 300, strand=1), type='terminator')]
    dict1 = {'name': 'ATF'}
    dict2 = {'name': 'ATF_tTEF'}
    features1[0].qualifiers = dict1
    features2[0].qualifiers = dict2

    TEMPLATES = [SeqRecord(seq = Seq('ATGATGA'*1000), id='seq_DPrPIMvy', name='ATF', description='<unknown description>', dbxrefs=[], features = features1+ features2)]

    # the test
    names_of_sites_we_want_to_amplify = ['ATF']
    name_of_terminator_site_to_be_incorporated= 'tTEF'

    extractions_sites = extract_template_amplification_sites(TEMPLATES, names_of_sites_we_want_to_amplify, name_of_terminator_site_to_be_incorporated)
    
    assert len(extractions_sites[0].seq) == 299
    assert len(extractions_sites[0].features) == 2
    assert str(extractions_sites[0].seq) == 'TGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATG'
    assert extractions_sites[0].name == 'ATF'
    


def test_extract_sites(): 
    """ This function will extract sequences from feature locations """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    features1 = [SeqFeature(FeatureLocation(2, 100, strand=1), type='promoter1')]
    features2 = [SeqFeature(FeatureLocation(200, 300, strand=1), type='promoter2')]
    features1[0].qualifiers['name'] =  'pCYC1'
    features2[0].qualifiers['name'] =  'pCYC1'

    prom_seq = Seq('CAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAAATATAGATAGATGATGATGATA'*50)

    prom_template = [SeqRecord(seq = prom_seq, id='seq_DPrPIMvy', name='pCYC1', description='<unknown description>', dbxrefs=[], features = features1+ features2)]
    p_annotations = 'pCYC1'
    p_names = 'pCYC1'

    promoter_sites =  extract_sites(p_annotations, prom_template, p_names)
    assert promoter_sites[0].seq == 'GCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAAATATAGATAGATGATGATGATACAGCATTTTCAAAGGTGTGTTCTT'
    assert promoter_sites[0].name == 'pCYC1'
    assert promoter_sites[1].seq =='ATGAAAATATAGATAGATGATGATGATACAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAAATATAGATAGATGATGAT'



def test_seq_to_annotation(): 
    template1 = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF1 = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR1 = 'ACCTGCACUTAACTAATTACATGACTCGA'
    rec_vector1 = pcr(primerF1,primerR1, template1)

    template2 = 'ATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTG'
    primerF2 = 'ATATTTATACAGAAACTTGATG'
    primerR2 = 'CAGAACATTGTAGGGTGTGAAAAAAT'
    rec_vector2 = pcr(primerF2,primerR2, template2)
    
    # The actual test
    seq_to_annotation(rec_vector2, rec_vector1, 'AMPLICON2')
    assert rec_vector1.features[2].type == 'AMPLICON2'
    assert rec_vector1.features[2].location.start.position == 44
    assert rec_vector1.features[2].location.end.position == 186

def test_UPandDW(): 
    pass 

def test_casembler():
    pass