#!/usr/bin/env python
# MIT License
# Copyright (c) 2022, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

""" A script that provides utility functions"""


def mean(lst:list):
    """Get mean from a list.
    Parameters
    ----------

    lst : list of floats.

    Return
    ------
    mean : float
    """
    mean = sum(lst) / len(lst)
    return mean

def counting_occurences(data_with_occurences:dict):
    """
    Count the occurences of each key in the input dict and returns the percentage of each key in the total values
    
    Parameters
    ----------
    data_with_occurences : dict
        The dictionary containing the data for counting occurences.
    Returns
    -------
    tuple
        A tuple containing two lists, the first one is the occurence percentage of each key, the second one is the list of keys.
    """
    columns = []
    data = []
    for key,value in data_with_occurences.items(): 

        values = data_with_occurences.values()
        total = sum(values)

        data.append((value/total)*100)
        columns.append(key)
        
    return data, columns