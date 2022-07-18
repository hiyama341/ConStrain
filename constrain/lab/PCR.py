#!/usr/bin/env python
# standard libraries
import os
import pathlib
import itertools
import string
import numpy as np
import pandas as pd
import pydna
import math
import pathlib
import textwrap as _textwrap
from pydna._pretty import pretty_str as _pretty_str
import csv


""" This part of the lab module is used for PCRs

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
    Example:

    3-reactions:
    -----------------------
                   vol_p_reac	vol_p_3_reac
    Template       	 1.0	        3.0
    Primer 1	     2.5	        7.5
    Primer 2	     2.5	        7.5
    H20	             19.0	        57.0
    MM	             25.0	        75.0
    Total	         50.0         	150.0
    -----------------------

    If this is used as input:
    ------------------------
    (vol_p_reac = 10,
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
    volumes_p_x_p_y_reactions

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
    This function determines process speed dependent on the which polymerase is used
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


def PCR_program(
    amplicon, pgroup, touch_down=False, ta_tm_refresh=False, primer_con=500
):
    """Investigates ta and tm for amplicons and creates PCR program.

    param: amplicon             eg. Pydna amplicon object
    param: pgroup               eg.  product group - telling which polymerase is used

    """

    # Set product group & polymerase(/kit)
    # Product groups: Q5, Q5 Hot Start, OneTaq, OneTaq Hot Start, Phusion, LongAmp Taq, LongAmp Hot Start Taq,
    if pgroup == "OneTaq Hot Start":
        pol = "OneTaq Hot Start 2X Master Mix with Standard Buffer"
        (
            proc_speed,
            In_den_temp,
            melt_temp,
            melt_time,
            elong_temp,
            final_temp,
            final_time,
        ) = (60, 94, 94, 30, 68, 68, 10)
    elif pgroup == "Q5 Hot Start":
        pol = "Q5 Hot Start High-Fidelity 2X Master Mix"
        proc_speed, In_den_temp, melt_temp, melt_time, elong_temp, final_temp = (
            30,
            98,
            98,
            10,
            72,
            72,
        )
        final_time = " 2"
    elif pgroup == "Phusion":
        pol = "Phusion High-Fidelity DNA Polymerase (HF Buffer)"
        (
            proc_speed,
            In_den_temp,
            melt_temp,
            melt_time,
            elong_temp,
            final_temp,
            final_time,
        ) = (30, 98, 98, 10, 72, 72, 10)
    else:
        print("The provided polymerase is not accounted for in this code.")
        print("""Write either "OneTaq Hot Start", "Q5 Hot Start" or "Phusion".""")

    # Extension time.
    extension_time_taq = int(proc_speed * len(amplicon) / 1000)

    # check if keys already exist
    tm1_key = "tm " + pgroup
    tm2_key = "tm " + pgroup
    ta_key = "ta " + pgroup
    if (
        tm1_key not in amplicon.forward_primer.annotations
        and tm2_key not in amplicon.reverse_primer.annotations
        and ta_key not in amplicon.annotations
    ) or ta_tm_refresh:

        # Extract primer sequence
        p1_seq = Dseqrecord(amplicon.forward_primer.footprint).seq.watson
        p2_seq = Dseqrecord(amplicon.reverse_primer.footprint).seq.watson

        if "polymerase" in amplicon.annotations:
            pass
        else:
            amplicon.annotations["polymerase"] = pgroup

        # Open TM calculator
        chrome_webdriver_path = (
            os.path.normpath(os.getcwd() + os.sep + os.pardir) + "/chromedriver"
        )
        driver = webdriver.Chrome(chrome_webdriver_path)
        driver.get("http://tmcalculator.neb.com/#!/main")
        assert "Calculator" in driver.title

        # Wait for webpage to load
        driver.implicitly_wait(5)

        ## Select product group and polymerase/kit
        prod_group = Select(
            driver.find_element_by_xpath('//*[@id="input"]/div[1]/div/select[1]')
        )
        prod_group.select_by_visible_text(pgroup)
        polymerase = Select(
            driver.find_element_by_xpath('//*[@id="input"]/div[1]/div/select[2]')
        )
        polymerase.select_by_visible_text(pol)

        # Change primer concentration
        prim_con = driver.find_element_by_id("ct")
        prim_con.clear()  # Clear bar incase of input
        prim_con.send_keys(str(primer_con))

        # Input primers and extract tm
        def input_primer_and_extract_tm(prim_no, prim_seq):
            # Input primer
            bar_id = "p" + str(prim_no)
            prim_bar = driver.find_element_by_id(bar_id)
            prim_bar.clear()  # Clear bar incase of input
            prim_bar.send_keys(prim_seq)

            # Extract Tm
            result_id = "tm" + str(prim_no)
            tm_elem = driver.find_element_by_id(result_id)
            tm_text = tm_elem.text
            pattern = (
                "^Primer " + str(prim_no) + "\\n\d+\snt\\n\d+%\sGC\\nTm:\s(-?\d+).C$"
            )
            m = re.match(pattern, tm_text)
            tm = int(m.group(1))
            return tm

        # Save melting temperature
        tm1 = input_primer_and_extract_tm(1, p1_seq)
        tm2 = input_primer_and_extract_tm(2, p2_seq)

        # Save annealing temperature
        ta_elem = driver.find_element_by_id("ta")
        ta_text = ta_elem.text
        ta_pattern = "^Anneal at\\n(-?\d+)\s.C$"
        m3 = re.match(ta_pattern, ta_text)
        ta = int(m3.group(1))

        # Close tmcalculator
        driver.close()

        amplicon.forward_primer.annotations[tm1_key] = tm1
        amplicon.reverse_primer.annotations[tm2_key] = tm2
        amplicon.annotations[ta_key] = ta
    else:
        ta = amplicon.annotations[ta_key]
        tm1 = amplicon.forward_primer.annotations[tm1_key]
        tm2 = amplicon.reverse_primer.annotations[tm2_key]

    # ta messages
    if ta < 45:
        print("Annealing temperature is lower than the recommended minimum of 45 °C.")
    elif ta > 72:
        print(
            "Annealing temperature for experiments with this enzyme should typically not exceed 72°C."
        )
    if abs(tm1 - tm2) > 5:
        print("Tm difference is greater than the recommended limit of 5 °C.")

    # PCR program
    if touch_down:
        f = _textwrap.dedent(
            r"""
        | x1 |       -1°Cx6      |        x24        | 1x  |
        |{initial_den}°C|{melting_temp}°C               |{melting_temp}°C               |     |tmf:{tmf:.1f}
        |____|_____          {el_tem}°C|_____          {el_tem}°C|{f_temp}°C |tmr:{tmr:.1f}
        |30s |{melting_time}s  \ {ta_start:.1f}°C _____|{melting_time}s  \ {ta_end:.1f}°C _____|_____|{rate}s/kb
        |    |      \______/{0:2}:{1:2}|      \______/{0:2}:{1:2}|{f_time}min|GC {GC_prod}%
        |    |       30s         |       30s         |     |{size}bp
                            """[
                1:-1
            ].format(
                initial_den=In_den_temp,
                melting_temp=melt_temp,
                melting_time=melt_time,
                el_tem=elong_temp,
                f_temp=final_temp,
                f_time=final_time,
                rate=proc_speed,
                size=len(amplicon.seq),
                ta_start=ta + 3,
                ta_end=ta - 3,
                tmf=tm1,
                tmr=tm2,
                GC_prod=int(amplicon.gc()),
                *map(int, divmod(extension_time_taq, 60)),
            )
        )
    else:
        f = _textwrap.dedent(
            r"""
                                |{initial_den}°C|{melting_temp}°C               |     |tmf:{tmf:.1f}
                                |____|_____          {el_tem}°C|{f_temp}°C |tmr:{tmr:.1f}
                                |30s |{melting_time}s  \ {ta:.1f}°C _____|____ |{rate}s/kb
                                |    |      \______/{0:2}:{1:2}|{f_time}min|GC {GC_prod}%
                                |    |       30s         |     |{size}bp
                                """[
                1:-1
            ].format(
                initial_den=In_den_temp,
                melting_temp=melt_temp,
                melting_time=melt_time,
                el_tem=elong_temp,
                f_temp=final_temp,
                f_time=final_time,
                rate=proc_speed,
                size=len(amplicon.seq),
                ta=ta,
                tmf=tm1,
                tmr=tm2,
                GC_prod=int(amplicon.gc()),
                *map(int, divmod(extension_time_taq, 60)),
            )
        )

    return _pretty_str(f)


def takeThird(elem):

    return elem[2]


def det_no_of_thermal_cyclers(amplicons, polymerase, elong_time_max_diff=15):

    """Determines the number of thermalcyclers that is needed based on elongation time differences"""

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


def pcr_locations(amplicons):
    """Save primer information for amplicons."""
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
    df = pd.DataFrame(
        list(zip(product_loc, product_names, template_loc, fw_loc, rv_loc)),
        columns=["location", "name", "template", "fw", "rv"],
    )

    return df


def nanophotometer_concentrations(path: pathlib.PosixPath):

    """This function reads a CSV file with nanophotometer concentraions and returns the concentrations in a list"""
    concentrations = []
    with open(path, encoding="Latin1") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        next(reader)[4]
        for row in reader:
            conc = float(row[4].replace(",", "."))
            concentrations.append(conc)

    return concentrations
