#MIT License
#Copyright (c) 2022, Technical University of Denmark (DTU)
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

'''For making Combinatorial libraries of DNA fragments'''

#!/usr/bin/env python
# standard libraries
import os
import pathlib
import itertools
import string
import numpy as np
import pandas as pd
import pydna

# For primer TM function
import requests
import json

# Biopython standards
import Bio
import Bio.SeqRecord
import Bio.Seq
import Bio.SeqFeature

# Pydna for the molecular bio
import pydna.editor
import pydna.primer
import pydna.dseqrecord
import pydna.amplify
import pydna.assembly
import pydna.gel
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.design import assembly_fragments
from pydna.assembly import Assembly
import itertools
from pydna.tm import tm_default as _tm_default

# For anealing temp
from Bio.SeqUtils import MeltingTemp as mt
from pydna.tm import tm_default as _tm_default
from Bio.Seq import Seq


def primer_tm_neb(primer):
    '''Calculates a single primers melting temp'''


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
    '''Calculates primer pair melting temp'''

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


""" This part of the design module is used for making combinatorial libraries

Parameters
----------
param: list_of_seqs:list of list      e.g. DNA sequences for the parts
param: list_of_names:list             e.g. Names for the parts
param: pad:str                        e.g. A string that can be incorporated into the build
param: position_of_pad:int            e.g. At wich position the pad should be inserted (zero-indexed)
-------

Returns:
list of all primers needed for the build
list of all PCRs needed to make the build
graphical representations of the library

-----
HELPER FUNCTIONS: 
-----------------
primer_tm_neb
primer_ta_neb
CombinatorialListMaker
systematic_names_function
EmptyListMaker
SimpleAmpliconMaker
GetPrimers
AssemblyMaker
UniquePrimers
UniqueAmplicons
MakingAssemblyObjects
MakingAssembledContigs
"""


def Combinatorial_List_Maker(listOflist_that_is_being_made_into_all_combinations):
    """
    Makes all possible combinations from a list of list

    Parameters
    ----------

    Returns:


    """
    combinations = list(
        itertools.product(*listOflist_that_is_being_made_into_all_combinations)
    )

    return combinations


def systematic_names_function(List_of_list_parts: list):
    """Returns a list of list with systematic names i.e [1,1,1], [1,2,1]... etc

    Parameters
    ----------
    """
    # The number of parts of each fragment
    no_parts = [int(len(l)) for l in List_of_list_parts]

    ### For naming the strains systematically ### basicly making a list from the number of parts with indexes
    list_of_systematic = []
    midlertidiglist = []
    for parts in no_parts:
        for j in range(0, parts):
            midlertidiglist.append(j + 1)
        list_of_systematic.append(midlertidiglist)
        midlertidiglist = []

    # Then we use itertools to make the right combinations
    combinatorial_list_of_indexes = list(itertools.product(*list_of_systematic))

    return combinatorial_list_of_indexes


def EmptyListMaker(list_of_sequences: list):
    """returns empty list in the length of seqs

    Parameters
    ----------


    """
    EmptyList = [[] for i in range(len(list_of_sequences))]

    return EmptyList


def SimpleAmpliconMaker(list_of_seqs: list, list_of_names: list):
    """Calculates amplicons, updates their names


    Parameters
    ----------

    """
    # Start by making an empty list
    list_of_amplicons = [[] for i in range(len(list_of_seqs))]
    list_of_amplicon_primers = [[] for i in range(len(list_of_seqs))]
    list_of_amplicon_primer_temps = [[] for i in range(len(list_of_seqs))]

    ### HERE WE CALCULATE Amplicons, primers, and their temperatures

    # Then we calculate the primers with the NEB calculator
    for i in range(0, len(list_of_seqs)):
        for j in range(0, len(list_of_seqs[i])):
            # Append Amplicons
            amplicons = primer_design(
                list_of_seqs[i][j], tm_func=_tm_default, target_tm=56.0, limit=13
            )  ############## Can add NEB Calculator here: primer_TM ################# _tm_default i.e tm_func = _tm_default,

            # Updating names
            amplicons.name = list_of_names[i][j]
            list_of_amplicons[i].append(amplicons)

            # Save the primers
            primers = (amplicons.forward_primer.seq, amplicons.reverse_primer.seq)
            list_of_amplicon_primers[i].append(primers)

            # Save melting temps
            ############## Can add NEB Calculator here: primer_TM #############################
            melting_temps = (
                _tm_default(amplicons.forward_primer.seq),
                _tm_default(amplicons.reverse_primer.seq),
            )
            list_of_amplicon_primer_temps[i].append(melting_temps)

    return list_of_amplicons, list_of_amplicon_primers, list_of_amplicon_primer_temps


def GetPrimers(
    List_of_assemblies: list,
    combinatorial_list_of_names: list,
    combinatorial_list_of_primer_tm: list,
):
    """Returns a list of ALL primers from the combinatorial library

    Parameters
    ----------


    """

    primers_temporary = []
    primers = []

    counter = 0
    for i in range(0, len(List_of_assemblies)):
        for j in range(0, len(List_of_assemblies[i])):
            counter += 1
            # Names
            List_of_assemblies[i][j].name = combinatorial_list_of_names[i][j]
            # Primers
            # description ------ DESCRIBES what other part it overlaps-------------
            if j == 0:  # START OF THE ASSEMBLY
                List_of_assemblies[i][
                    j
                ].forward_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                )
                List_of_assemblies[i][j].reverse_primer.description = (
                    "Anneals to "
                    + str(List_of_assemblies[i][j].name)
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j + 1].name)
                )
            if j > 0 and j < len(List_of_assemblies[i]) - 1:  #      # THE rest:
                List_of_assemblies[i][
                    j
                ].forward_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j - 1].name)
                )
                List_of_assemblies[i][
                    j
                ].reverse_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j + 1].name)
                )
            if j == len(List_of_assemblies[i]) - 1:  # THE END OF THE ASSEMBLY
                List_of_assemblies[i][j].forward_primer.description = (
                    "Anneals to "
                    + str(List_of_assemblies[i][j].name)
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j - 1].name)
                )
                List_of_assemblies[i][
                    j
                ].reverse_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                )

            # template it aneals to
            List_of_assemblies[i][j].forward_primer.name = str(
                List_of_assemblies[i][j].name
            )
            List_of_assemblies[i][j].reverse_primer.name = str(
                List_of_assemblies[i][j].name
            )

            # Primer tm
            List_of_assemblies[i][j].forward_primer.features = round(
                float(combinatorial_list_of_primer_tm[i][j][0]), 2
            )
            List_of_assemblies[i][j].reverse_primer.features = round(
                float(combinatorial_list_of_primer_tm[i][j][1]), 2
            )

            fwd_rev_primers = [
                List_of_assemblies[i][j].forward_primer,
                List_of_assemblies[i][j].reverse_primer,
            ]
            primers_temporary.append(fwd_rev_primers)

        primers.append(primers_temporary)
        primers_temporary = []

    return primers


def AssemblyMaker(combinatorial_list_of_amplicons: list):
    """Assembles Amplicons with pad and makes new overlapping primers

    Parameters
    ----------

    """

    List_of_assemblies = []
    for i in range(0, len(combinatorial_list_of_amplicons)):
        List_of_assemblies.append(
            assembly_fragments(
                combinatorial_list_of_amplicons[i], overlap=35, maxlink=40
            )
        )

    return List_of_assemblies


def UniquePrimers(primers: list, list_of_assemblies):

    """ Finds unique primers from a list of assemblies
    Parameters
    ----------

    """

    unikke_F_primers = []
    unikke_R_primers = []
    length_of_unique_primers = 0
    counter = 0
    primer_list = []

    for i in range(0, len(primers)):
        for j in range(0, len(primers[i])):
            counter += len(primers[i][j])
            if primers[i][j][0] not in unikke_F_primers:
                unikke_F_primers.append(primers[i][j][0])
            if primers[i][j][1] not in unikke_R_primers:
                unikke_R_primers.append(primers[i][j][1])

    counter = 0
    unique_forward_primers = []
    unique_reverse_primers = []

    ### CHANGING THE NAMES OF THE PRIMERS
    for i in range(len(unikke_F_primers)):
        counter += 1
        unikke_F_primers[i].id = "F{number:03}".format(number=counter)
        length_of_unique_primers += len(unikke_F_primers[i].seq)
        U_f_primers = [
            unikke_F_primers[i].id,
            unikke_F_primers[i].name,
            unikke_F_primers[i].seq,
            unikke_F_primers[i].features,
            len(unikke_F_primers[i].seq),
            len(unikke_F_primers[i].seq) * 1.8,
        ]
        unique_forward_primers.append(U_f_primers)

    for i in range(len(unikke_R_primers)):
        counter += 1
        unikke_R_primers[i].id = "R{number:03}".format(number=counter)
        length_of_unique_primers += len(unikke_R_primers[i].seq)
        U_r_primers = [
            unikke_R_primers[i].id,
            unikke_R_primers[i].name,
            unikke_R_primers[i].seq,
            unikke_R_primers[i].features,
            len(unikke_R_primers[i].seq),
            len(unikke_R_primers[i].seq) * 1.8,
        ]
        unique_reverse_primers.append(U_r_primers)

    primer_list = (
        unique_forward_primers + unique_reverse_primers
    )  # cOULD CONCATONATE THEM INTO: unique_forward_primers + unique_reverse_primers

    ### Updating primer names and removing duplicates?
    for i in range(0, len(list_of_assemblies)):
        for j in range(0, len(list_of_assemblies[i])):
            for l in range(0, len(unikke_F_primers)):
                if (
                    list_of_assemblies[i][j].forward_primer.seq
                    == unikke_F_primers[l].seq
                ):
                    list_of_assemblies[i][j].forward_primer = unikke_F_primers[l]
            for m in range(0, len(unique_reverse_primers)):
                if (
                    list_of_assemblies[i][j].reverse_primer.seq
                    == unikke_R_primers[m].seq
                ):
                    list_of_assemblies[i][j].reverse_primer = unikke_R_primers[m]

    return primer_list


def UniqueAmplicons(list_of_assemblies: list):

    """ Finds Unique amplicons from a list of assemblies
    Parameters
    ----------
    """
    ### Unique amplicons
    unique_amplicons = []
    for i in range(0, len(list_of_assemblies)):
        for j in range(0, len(list_of_assemblies[i])):
            if list_of_assemblies[i][j] not in unique_amplicons:
                unique_amplicons.append(list_of_assemblies[i][j])

    return unique_amplicons


def MakingAssemblyObjects(list_of_assemblies: list):
    """
    Assembling amplicons into assembling class that shows fragments, limit,
    nodes and which algorithm that was used for assembling. 


    Parameters
    ----------


    """
    list_of_assembly_objects = []
    for i in range(0, len(list_of_assemblies)):
        list_of_assembly_objects.append(Assembly((list_of_assemblies[i]), limit=35))

    return list_of_assembly_objects


def MakingAssembledContigs(list_of_assembly_objects: list):
    """ Assembles a list of assembly object into contigs
    
    Input: 
    :param list of assembly objects
    ----------
    
    Returns: 
    list of contigs
    """
    contigs_assembled = []
    for j in range(0, len(list_of_assembly_objects)):
        contigs_assembled.append(list_of_assembly_objects[j].assemble_linear())

    return list_of_assembly_objects


class DesignAssembly:
    """Class representing a combinatorial design library of chosen constructs.

    Input:
    :param list_of_seqs: A list of list of a construct of choice.
    :param list_of_names: A list of list of the names wanted for the construct of choice.
    :param pad: Volume to be transferred, expressed in liters.
    :param position_of_pad: A dict containing any useful information about the transfer.
        This information can be used later e.g. as parameters for the printing out primers,PCRs needed
        for the construction of the library

    Design_output:
    The rest of the parameters are generated from the initial such as:
    :param list_of_amplicons
    :param systematic_names
    :param combinatorial_list_of_amplicons
    :param combinatorial_list_of_names
    :param combinatorial_list_of_primer_tm
    :param list_of_assemblies
    :param primers
    :param unique_primers
    :param unique_amplicons

    """

    def __init__(
        self, list_of_seqs: list, list_of_names: list, pad: str, position_of_pad: int
    ):

        ###  1.INITIALIZING ##
        self.list_of_seqs = list_of_seqs
        self.list_of_names = list_of_names
        self.pad = pad
        self.position_of_pad = position_of_pad

        ### 2. Amplicons, primers, and their temperatures
        (
            self.list_of_amplicons,
            self.list_of_amplicon_primers,
            self.list_of_amplicon_primer_temps,
        ) = SimpleAmpliconMaker(self.list_of_seqs, self.list_of_names)

        # Systematic names
        self.systematic_names = systematic_names_function(self.list_of_seqs)

        ### 3. COMBINATORIAL LISTS
        self.combinatorial_list_of_amplicons = Combinatorial_List_Maker(
            self.list_of_amplicons
        )
        self.combinatorial_list_of_names = Combinatorial_List_Maker(self.list_of_names)
        self.combinatorial_list_of_primer_tm = Combinatorial_List_Maker(
            self.list_of_amplicon_primer_temps
        )

        # Making the combinations into a list so we can insert PADS later (They are tuples at this stage, and insert doesnt work for tuples)
        for i in range(0, len(self.combinatorial_list_of_amplicons)):
            self.combinatorial_list_of_amplicons[i] = list(
                self.combinatorial_list_of_amplicons[i]
            )

        #### 4. Adding PAD ###
        for i in range(0, len(self.combinatorial_list_of_amplicons)):
            self.combinatorial_list_of_amplicons[i].insert(
                self.position_of_pad, self.pad
            )

        ### 5. Assembling and making overlapping primers
        self.list_of_assemblies = AssemblyMaker(self.combinatorial_list_of_amplicons)

        ### 6. GETTING all primers, annotating, adding features
        self.primers = GetPrimers(
            self.list_of_assemblies,
            self.combinatorial_list_of_names,
            self.combinatorial_list_of_primer_tm,
        )

        ### 7. Getting Unique primers and re-annotating list_assemblies to get right names
        self.unique_primers = UniquePrimers(self.primers, self.list_of_assemblies)

        ### 8. Unique amplicons
        self.unique_amplicons = UniqueAmplicons(self.list_of_assemblies)

    ##########################################################
    ###################  CLASS METHODS  ######################

    def ShowContigs(self):
        """Returns a string of the contigs generated by the assembly"""
        print("Template, Primer, tm")
        for i in range(0, len(self.list_of_assemblies)):
            print("\nContig" + str(self.systematic_names[i]))
            for j in range(0, len(self.list_of_assemblies[i])):
                print(
                    "Template: ", self.list_of_assemblies[i][j].name[0:15]
                )  # , '\t', self.primers[i][j][0].name,'\t',self.primers[i][j][0].features)
        return

    def ShowVariantsLibDF(self):
        """Returns a dataframe of all the variants"""
        combinatorial_lib_variants_df = pd.DataFrame(self.combinatorial_list_of_names)
        systematic_names = self.systematic_names
        combinatorial_lib_variants_df["Systematic_name"] = systematic_names
        combinatorial_lib_variants_df["Variant"] = np.arange(
            len(combinatorial_lib_variants_df)
        )

        return combinatorial_lib_variants_df

    def print_primer_list(self):
        """Return the list of transfers in human-readable format."""
        for primers in self.unique_primers:
            print(primers)

    def primer_list(self):
        """Return the list of transfers in human-readable format."""
        primer_list = []
        for primers in self.unique_primers:
            primer_list.append(primers)

        return primer_list

    def primer_list_to_dataframe(self):
        """Return a pandas dataframe with list of primers."""
        df = pd.DataFrame(self.unique_primers)
        df.columns = [
            "ID",
            "Anneals to",
            "Sequence",
            "Annealing temperature",
            "Length",
            "Price(DKK)",
        ]

        return df

    def print_PCR_list(self):
        """Prints PCR_list"""
        print("PCR#, Template,forward_primer, reverse primer, F_tm, R_tm")
        for i in range(0, len(self.unique_amplicons)):
            print(
                "PCR{number}".format(number=i + 1),
                ",",
                self.unique_amplicons[i].name,
                ",",
                self.unique_amplicons[i].forward_primer.id,
                ",",
                self.unique_amplicons[i].reverse_primer.id,
                ",",
                self.unique_amplicons[i].forward_primer.features,
                ",",
                self.unique_amplicons[i].reverse_primer.features,
            )

    def PCR_list(self):
        """Returns a PCR_list"""
        pcr_list = []
        for i in range(0, len(self.unique_amplicons)):
            PCR = [
                "PCR{number}".format(number=i + 1),
                self.unique_amplicons[i].name,
                self.unique_amplicons[i].forward_primer.id,
                self.unique_amplicons[i].reverse_primer.id,
                self.unique_amplicons[i].forward_primer.features,
                self.unique_amplicons[i].reverse_primer.features,
            ]
            pcr_list.append(PCR)

        return pcr_list

    def PCR_list_to_dataframe(self):
        """Prints PCR_list into a pandas dataframe"""
        dataframe_list = []
        for i in range(0, len(self.unique_amplicons)):
            lst = [
                "PCR{number}".format(number=i + 1),
                self.unique_amplicons[i].name,
                self.unique_amplicons[i].forward_primer.id,
                self.unique_amplicons[i].reverse_primer.id,
                self.unique_amplicons[i].forward_primer.features,
                self.unique_amplicons[i].reverse_primer.features,
            ]
            dataframe_list.append(lst)

        df = pd.DataFrame(dataframe_list)
        df.columns = [
            "PCR#",
            "Template",
            "forward_primer",
            "reverse_primer",
            "F_tm",
            "R_tm",
        ]

        return df

    def GraphicalRepresentationOfAssemblies(self):
        """
        Takes in the assembly object and returns graphical report of the fragments assembled
        """
        graphical_representation = [
            self.assembly_object[x].assemble_linear()[0].figure()
            for x in range(0, len(self.assembly_object))
        ]

        return graphical_representation
