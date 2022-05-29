#!/usr/bin/env python

# Test Transformation module

import sys
import pytest

# Importing the module we are  testing
from constrain.lab.transformation import *


def test_ng_to_nmol():
    
    ''' Testing that it calculates correctly'''
    
    #increasing ng
    assert ng_to_nmol(100, 1000) == 0.00015384615384615385
    assert ng_to_nmol(1000, 1000) == 0.0015384615384615385
    
    #increaseing basepair
    assert ng_to_nmol(100, 5000) == 3.076923076923077e-05
    assert ng_to_nmol(100, 10000) == 1.5384615384615384e-05
    
    # Non-valid input
    assert ng_to_nmol(-500, 8000) == 'non-valid_input'
    assert ng_to_nmol(500, -9000) == 'non-valid_input'
    
    # division by zero
    assert ng_to_nmol(0, 8000) == 'non-valid_input'
    assert ng_to_nmol(500, 0) == 'non-valid_input'



def test_ODtime(): 
    
    assert ODtime(1, 1) == 1.414

    assert ODtime(0, 1) == 0.0

    assert ODtime(1, 0) == 1.0
    assert ODtime(1, 10) == 32.0
    
    #non-valid input
    assert ODtime(-1, 1) == 'non-valid_input'
    assert ODtime(1, -1) == 'non-valid_input'
    
    
def test_time_to_inoculate():
    pass
    
    
def test_transformation_mix(): 

    pass
 
