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

""" Easy to use benchling functions to fetch sequences and objects"""


#!/usr/bin/env python
import os
from benchlingapi import Session
import pandas as pd
import datetime
import Bio
from Bio.SeqFeature import SeqFeature
import pydna

# Fetching API key from env file  - the function can find the env file in the root folder
import os
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

# Initializing session
api_key = os.environ.get("API_KEY")

# Add your benchling home environment here
home_url = os.environ.get("HOME_url")

# print(home_url, api_key)

session = Session(api_key=api_key, home=home_url)


""" This part of the LIMS module is used for importing sequences and exporting sequences to benchling

-----------------
HELPER FUNCTIONS: 
-----------------
sequence_to_benchling
from_benchling
update_loc_vol_conc
unnest_dict
nest_dict
SE_to_CLoc
split_based_on_keys
rename_dict_keys
Compound_to_SE()
"""


def sequence_to_benchling(folder_name, oligo_name, oligo_bases, schema):

    """This function uploads sequences to benchlong

    Parameters
    ----------

    Returns
    -------

    """

    folder = session.Folder.find_by_name(folder_name)
    dna = ""
    dna = session.DNASequence(
        name=oligo_name, bases=oligo_bases, folder_id=folder.id, is_circular=False
    )

    # save the dna to Benchling account
    dna.save()

    schemas = (
        "Primer",
        "DNA Fragment",
        "Plasmid",
        "Gene",
        "gRNA",
        "Marker",
        "Promoter",
        "Terminator",
        "Tag",
        "Origin of Replication",
    )
    if schema in schemas:
        dna.set_schema(schema)

    dna.register()


def from_benchling(bname, schema=""):

    """Extract information of object on benchling.


    Parameters
    ----------

    Returns
    -------
    """
    # retrieve benchling sequence dict
    bench_dict = session.DNASequence.find_by_name(bname).dump()

    # Convert benchling sequence dict into Biopython/Genbank object
    ## rename keys from Benchling to Biopython/Genbank
    trans_dict = {
        "bases": "seq",
        "id": "id",
        "annotations": "features",
        "name": "name",
        "fields": "annotations",
    }
    trans_bench_dict = rename_dict_keys(bench_dict, trans_dict)

    ## Select most import information
    translated_bench_dict_sel, translated_bench_dict_other = split_based_on_keys(
        trans_bench_dict, trans_dict.values()
    )

    ## Add Benchling custom fields and isCircular to fields dict in Biopython annotations.
    translated_bench_dict_sel["annotations"].update(
        translated_bench_dict_other["customFields"]
    )
    translated_bench_dict_sel["annotations"].update(
        {"topology": translated_bench_dict_other["isCircular"]}
    )

    ##
    date = datetime.datetime.today().strftime("%d-%b-%Y").upper()
    # comment is renamed to commentary as Bio.IO.write provides an error when this key is in the annotation dict.
    comment = translated_bench_dict_sel["annotations"].pop("comment", None)

    translated_bench_dict_sel["annotations"].update(
        {
            "data_file_division": "PLN",
            "date": date,
            "molecule_type": "DNA",
            "location": "unknown",
            "commentary": comment,
        }
    )

    ## Convert values
    ### Seq
    translated_bench_dict_sel["seq"] = Bio.Seq.Seq(translated_bench_dict_sel["seq"])

    ### Features
    seq_length = len(translated_bench_dict_sel["seq"])

    #### Translate start end to compound locations object
    translated_bench_dict_sel["features"] = [
        SE_to_CLoc(feature_dict, seq_length)
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    #### Add "label" feature.qualifier dict for CLC workbench visualization
    for feature in translated_bench_dict_sel["features"]:
        feature.update({"label": feature.get("name", "")})

    #### Move non Biopython features to qualifier subdict
    translated_bench_dict_sel["features"] = [
        nest_dict(
            feature_dict,
            first_order_keys=["location", "type", "strand"],
            key_for_nested_dict="qualifiers",
        )
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    #### Create SeqFeatures
    translated_bench_dict_sel["features"] = [
        Bio.SeqFeature.SeqFeature(**feature_dict)
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    seqRecord = Bio.SeqRecord.SeqRecord(**translated_bench_dict_sel)
    seqRecord.name = bname
    seqRecord = update_loc_vol_conc(seqRecord)

    if schema == "Primer":
        seqRecord = pydna.primer.Primer(seqRecord)

    # return object
    return seqRecord


def update_loc_vol_conc(seqRecord):
    """Update with location volume and concentration
    information downloaded from benchling if possible.


    Parameters
    ----------

    Returns
    -------
    take Bio.SeqRecord.SeqRecord

    """

    DBpath = "../data/external/benchlingBoxes_sorted.csv"

    DB = pd.read_csv(DBpath)

    seqRecord.annotations["batches"] = []
    for i, row in DB.iterrows():
        if row["batchEntId"] == seqRecord.id:
            location = row["parentBoxPlateName"] + "_" + row["parentBoxPlatePos"]
            seqRecord.annotations["batches"].append(
                {
                    "box": row["parentBoxPlateName"],
                    "position": row["parentBoxPlatePos"],
                    "volume": int(row["volume"]),
                    "concentration": int(row["Concentration (ng/ul)"]),
                    "location": location,
                }
            )
        else:
            pass
    return seqRecord


def unnest_dict(dictionary, key_unnest_dict):
    """

    Parameters
    ----------

    Returns
    -------
    """

    nested_dict = dictionary.pop(key_unnest_dict)

    unnested_dictionary = {**dictionary, **nested_dict}

    if "label" in unnested_dictionary.keys():
        unnested_dictionary["name"] = unnested_dictionary.pop("label")

    return unnested_dictionary


def nest_dict(dictionary, key_for_nested_dict, first_order_keys=None):
    """

    Parameters
    ----------

    Returns
    -------
    """
    if first_order_keys is None:
        first_order_keys = []
    else:
        first_order_keys = first_order_keys

    first_order_dict, second_order_dict = split_based_on_keys(
        dictionary, first_order_keys
    )

    # Create qualifer dict
    dict_to_be_nested = {}
    dict_to_be_nested[key_for_nested_dict] = second_order_dict

    # Merge dictionaries by nesting qualifer dict
    nested_dictionary = {**first_order_dict, **dict_to_be_nested}

    return nested_dictionary


def SE_to_CLoc(dictionary, length):
    """Start and End Key Value pair to Compound Location Key Value pair.



    Parameters
    ----------

    Returns
    -------
    """
    start = dictionary.pop("start")
    end = dictionary.pop("end")
    start_pos = Bio.SeqFeature.ExactPosition(start)
    end_pos = Bio.SeqFeature.ExactPosition(end)
    if start_pos < end_pos:
        location = Bio.SeqFeature.FeatureLocation(start_pos, end_pos)
    else:
        f1 = Bio.SeqFeature.FeatureLocation(start_pos, length - 1)
        f2 = Bio.SeqFeature.FeatureLocation(0, end_pos)
        location = Bio.SeqFeature.CompoundLocation([f1, f2])

    dictionary["location"] = location
    return dictionary


def split_based_on_keys(dictionary, key_list):
    """Split dictionary based on key list"""
    # Split keys
    first_keys = key_list
    other_keys = [k for k in dictionary.keys() if k not in first_keys]

    # Create dicts
    first_dict = {k: dictionary[k] for k in first_keys}
    other_dict = {k: dictionary[k] for k in other_keys}

    return (first_dict, other_dict)


def rename_dict_keys(dictionary, trans_dictionary):
    """rename the keys of a dictionary using another dictionary


    Parameters
    ----------

    Returns
    -------
    """

    keys = dictionary.keys()
    values = dictionary.values()

    new_keys = [trans_dictionary.get(k, k) for k in keys]
    return dict(zip(new_keys, values))


def Compoundloc_to_SE(dictionary):
    """location :

    Parameters
    ----------

    Returns
    -------

     CompoundLocation to start : v ,
     end : v
     strand : v key value pairs"""

    Compoundloc = dictionary.pop("location")

    dictionary["start"] = Compoundloc.start.real
    dictionary["end"] = Compoundloc.end.real
    dictionary["strand"] = Compoundloc.strand

    return dictionary
