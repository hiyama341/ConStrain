#!/usr/bin/env python

# Test RobotAssembly
import sys
import pytest
import pandas as pd



from constrain.design.robot_assembly import RobotAssembly

# Read in test data
test_data = pd.read_excel('test_files/Random_PCR_list.xlsx')

# Add test primers and templates
fwd_primers = ['F0069','F0078', 'F0083','F0088','F0093','F0098','F00103','F00108']
rev_primers = ['R00174','R00175','R00176','R00177']
templates = ['pRPL15B']

# Initiate the Robotassembly class
test_Robotassemly = RobotAssembly(test_data,fwd_primers ,rev_primers, templates)


def test_Robotassemly_length_of_actions():
    # there should be 32(rows)*3*colums*5actions instructions here
    len_of_picklist = (test_Robotassemly.picklist.to_flowbot_instructions()).split()

    assert len(len_of_picklist) == 32*3*5


def test_Robotassemly_correct_input():
    assert type(test_Robotassemly.pandas_PCR) == pd.core.frame.DataFrame
    assert type(test_Robotassemly.forward_primers) == list
    assert type(test_Robotassemly.reverse_primers )  == list
    assert type(test_Robotassemly.templates ) == list


def test_Robotassemly_something():
    pass
