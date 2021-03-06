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

""" This part of the lab module is used for simulating and calculating PCR reactions"""
"""
HELPER FUNCTIONS:
-----------------
pcr_volumes (formerly known as volumes)
det_proc_speed
det_elon_time
PCR_program
det_no_of_thermal_cyclers
pcr_locations1
transf_locations1

"""

#!/usr/bin/env python
# standard libraries
import pathlib

import textwrap as _textwrap
import math
import csv
import pathlib
import json

# Extra
import pandas as pd
from pydna._pretty import pretty_str as _pretty_str
import requests


def primer_tm_neb(primer):
    """Calculates a single primers melting temp"""

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer]]

    input = {"seqpairs": seqpairs, "conc": 0.5, "prodcode": "q5-0"}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["tm1"]
    else:
        print("request failed")
        print(r["error"][0])


def primer_ta_neb(primer1, primer2):
    """Calculates primer pair melting temp"""

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer1, primer2]]

    input = {"seqpairs": seqpairs, "conc": 0.5, "prodcode": "q5-0"}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["ta"]

    else:
        print("request failed")
        print(r["error"][0])


def grouper(iterable, max_diff):
    """Groups objects into distinct groups based on differences"""
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= max_diff:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group


def pcr_volumes(
    vol_p_reac=0, no_of_reactions=1, standard_reagents=[], standard_volumes=[]
):
    """Can make a reaction scheme for PCR master mixes.


     Examples
     --------

     If this is used as input:
     ------------------------
    pcr_volumes(vol_p_reac = 10,
                no_of_reactions = 6,
                standard_reagents = ["DNA","Buffer, Cutsmart","H20","Enz, USER"],
                standard_volumes = [1,1,7,1])

     The following reaction scheme will be made:
     -------------------------
                     vol_p_reac	vol_p_x_reac
     DNA	                1.0	        6.0
     Buffer, Cutsmart	1.0	        6.0
     H20	                7.0	        42.0
     Enz, USER	        1.0	        6.0
     Total	            10.0      	60.0
     -------------------------

    """
    standard_total_volume = sum(standard_volumes)
    volumes_p_x = [val / standard_total_volume * vol_p_reac for val in standard_volumes]
    volumes_p_x_p_y_reactions = [val * no_of_reactions for val in volumes_p_x]

    volumes_p_x_plus_total = volumes_p_x + [sum(volumes_p_x)]
    volumes_p_x_p_y_reactions_plus_total = volumes_p_x_p_y_reactions + [
        sum(volumes_p_x_p_y_reactions)
    ]
    reagents_plus_total = standard_reagents + ["Total"]

    volumes_df = pd.DataFrame(
        data={
            "vol_p_reac": volumes_p_x_plus_total,
            "vol_p_"
            + str(no_of_reactions)
            + "_reac": volumes_p_x_p_y_reactions_plus_total,
        },
        index=reagents_plus_total,
    )
    return volumes_df


def det_proc_speed(amplicon):
    """
    This function determines process speed based on the which polymerase is used
    """

    if "proc_speed" in amplicon.forward_primer.annotations:
        print("proc_speed already set")
        return amplicon

    # proc_speed units are seconds/kb
    if amplicon.annotations["polymerase"] == "OneTaq Hot Start":
        proc_speed = 60
    elif amplicon.annotations["polymerase"] == "Q5 Hot Start":
        proc_speed = 30
    elif amplicon.annotations["polymerase"] == "Phusion":
        proc_speed = 30

    amplicon.annotations["proc_speed"] = proc_speed

    return amplicon


def det_elon_time(amplicon):

    """This function determines elongation time for an amplicon
    and add the elongation time to the amplicon annotations"""

    if "elongation_time" in amplicon.forward_primer.annotations:
        print("elongation_time already set")
        return amplicon

    # elongation_time units are seconds
    elongation_time = amplicon.annotations["proc_speed"] * len(amplicon) / 1000
    amplicon.annotations["elongation_time"] = math.ceil(elongation_time)

    return amplicon


def takeThird(elem):
    """Takes third element of a list"""
    return elem[2]


def det_no_of_thermal_cyclers(amplicons, polymerase, elong_time_max_diff=15):

    """Determines the number of thermalcyclers that is needed
    based on elongation time differences"""

    amp_names = [amplicon.name for amplicon in amplicons]
    elong_times = [amplicon.annotations["elongation_time"] for amplicon in amplicons]
    tas = [amplicon.annotations["ta " + polymerase] for amplicon in amplicons]
    order = list(range(0, len(amplicons)))

    list_of_tuples = list(zip(amp_names, tas, elong_times, order))

    list_of_tuples.sort(key=takeThird)

    groups = dict(enumerate(grouper(elong_times, elong_time_max_diff), 1))

    list_of_lists = [list(elem) for elem in list_of_tuples]

    for gNo, gTimes in groups.items():
        # print(gNo, gTimes)
        for idx, lst in enumerate(list_of_lists):
            if lst[2] in gTimes:
                list_of_lists[idx][2] = max(gTimes)

    thermal_cyclers = pd.DataFrame(
        list_of_lists, columns=["amplicons", "tas", "elong_times", "order"]
    )
    thermal_cyclers = thermal_cyclers.sort_values(["order"])
    thermal_cyclers = (
        thermal_cyclers.groupby(["tas", "elong_times"])["amplicons"]
        .apply(", ".join)
        .reset_index()
    )

    return thermal_cyclers


def pcr_locations(amplicons: list):
    """Obtain primer information for amplicons.

    Parameters
    ----------
    amplicon : list
        List of amplicon objects `pydna.amplicon`() # check this

    Returns
    -------
    pd.DataFrame
        Pandas dataframe with locations of your amplicons
    """
    # initialization
    product_loc = []
    product_names = []
    template_loc = []
    fw_loc = []
    rv_loc = []

    for i in range(0, len(amplicons)):
        product_names.append(amplicons[i].name)

        # Test if batches is present
        if (
            "batches" in amplicons[i].template.annotations.keys()
            and len(amplicons[i].template.annotations["batches"]) != 0
        ):
            product_loc.append(
                amplicons[i].template.annotations["batches"][0]["location"]
            )
            template_loc.append(
                amplicons[i].template.annotations["batches"][0]["location"]
            )
        elif (
            "batches" in amplicons[i].annotations.keys()
            and len(amplicons[i].annotations["batches"]) != 0
        ):
            product_loc.append(amplicons[i].annotations["batches"][0]["location"])
            template_loc.append(amplicons[i].annotations["batches"][0]["location"])
        else:
            product_loc.append("Empty")
            template_loc.append("Empty")
            print(
                "No batches were found for "
                + str(amplicons[i].name)
                + ". Please check the object."
            )

        # Save primer locations
        if (
            "batches" in amplicons[i].forward_primer.annotations.keys()
            and len(amplicons[i].forward_primer.annotations["batches"]) != 0
        ):
            fw_loc.append(
                amplicons[i].forward_primer.annotations["batches"][0]["location"]
            )
        else:
            fw_loc.append("Empty")
            print(str(amplicons[i].name) + ": Foward primer location was not found")
        if (
            "batches" in amplicons[i].reverse_primer.annotations.keys()
            and len(amplicons[i].reverse_primer.annotations["batches"]) != 0
        ):
            rv_loc.append(
                amplicons[i].reverse_primer.annotations["batches"][0]["location"]
            )
        else:
            rv_loc.append("Empty")
            print(str(amplicons[i].name) + ": Reverse primer location was not found")

    # Save information as dataframe
    df_pcr = pd.DataFrame(
        list(zip(product_loc, product_names, template_loc, fw_loc, rv_loc)),
        columns=["location", "name", "template", "fw", "rv"],
    )

    return df_pcr


def nanophotometer_concentrations(path: pathlib.PosixPath):

    """This function reads a CSV file with nanophotometer concentraions
    and returns the concentrations in a list"""
    concentrations = []
    with open(path, encoding="Latin1") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        next(reader)[4]
        for row in reader:
            conc = float(row[4].replace(",", "."))
            concentrations.append(conc)

    return concentrations


def simple_PCR_program(amplicon):

    """Simple PCR program designed to give a quick visual representations"""
    amplicon = det_elon_time(amplicon)
    amplicon = det_proc_speed(amplicon)

    # ta
    amplicon.annotations["ta Q5 Hot Start"] = primer_ta_neb(
        str(amplicon.forward_primer.seq), str(amplicon.reverse_primer.seq)
    )

    # tm forward and reverse
    amplicon.forward_primer.annotations["tm Q5 Hot Start"] = primer_tm_neb(
        str(amplicon.forward_primer.seq)
    )
    amplicon.reverse_primer.annotations["tm Q5 Hot Start"] = primer_tm_neb(
        str(amplicon.reverse_primer.seq)
    )

    r"""Returns a string containing a text representation of a suggested
    PCR program using Taq or similar polymerase.
    ::
     |98??C|98??C               |    |tmf:59.5
     |____|_____          72??C|72??C|tmr:59.7
     |30s |10s  \ 59.1??C _____|____|30s/kb
     |    |      \______/ 0:32|5min|GC 51%
     |    |       30s         |    |1051bp
    """

    formated = _textwrap.dedent(
        r"""
                            |98??C|98??C               |    |tmf:{tmf:.1f}
                            |____|_____          72??C|72??C|tmr:{tmr:.1f}
                            |30 s|10s  \ {ta:.1f}??C _____|____|{rate}s/kb
                            |    |      \______/{0:2}:{1:2}|2min|GC {GC_prod}%
                            |    |       20s         |    |{size}bp
                            """[
            1:-1
        ].format(
            rate=amplicon.annotations["proc_speed"],
            size=len(amplicon.seq),
            ta=amplicon.annotations["ta Q5 Hot Start"],
            tmf=amplicon.forward_primer.annotations["tm Q5 Hot Start"],
            tmr=amplicon.reverse_primer.annotations["tm Q5 Hot Start"],
            GC_prod=int(amplicon.gc()),
            *map(int, divmod(amplicon.annotations["elongation_time"], 60)),
        )
    )

    return _pretty_str(formated)
