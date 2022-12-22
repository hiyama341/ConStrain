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

""" This part of the design module is used fetching sequences"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import pearsonr
import scipy as sp
from scipy import stats
import matplotlib as mpl


# for phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
import pylab
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib as mpl


def plot_stacked_barplot(
    pd_dataframe_cross_tab_prop,
    colorDict: dict,
    save_pdf=True,
    path="",
    xlabel="",
    ylabel="",
) -> None:
    """Plotting stacked barplots from a pandas dataframe cross tab df

    Parameters
    ----------
    pd_dataframe_cross_tab_prop : pd.DataFrames
        pd.crosstab
    colorDict : Dict
        key as coloumn name, value hex color code
    save_pdf : bool
    path : string

    Returns
    -------
    stacked barplot
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Plot
    ax = pd_dataframe_cross_tab_prop.plot(
        kind="bar", stacked=True, figsize=(20, 6), color=colorDict, width=1.0
    )

    plt.legend(loc="upper right", ncol=2)
    plt.xlabel(xlabel, size=25, fontname="Helvetica", fontweight="bold")
    plt.ylabel(ylabel, size=25, fontname="Helvetica", fontweight="bold")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    plt.xscale("linear")

    # remove y axis labes
    plt.yticks([])
    plt.xticks(rotation=0, fontweight="bold", size=20)  # changing x scale by own

    ## size matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(25, 15)

    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=120, bbox_inches="tight")

    plt.show()


def plot_ml_learning_curve(
    x_partitioned_data: list,
    y_training: list,
    y_cv: list,
    training_sd,
    cv_sd: list,
    save_pdf=True,
    path="",
) -> None:
    """Plotting a learning curve from partitioned dataframes.

    Parameters
    ----------
    x_partitioned_data : list
        x coordinates
    y_training : list
        y coordinates for trained models
    y_cv : list
        y coordinates for cross-validated models
    training_sd : list
        standard-deviation for the models
    cv_sd : list
        standard-deviation for cross-validated models

    save_pdf : bool
    path : str

    Returns
    -------
    Learning curve
    """
    # making sd into numpy array
    training_sd = np.array(training_sd)
    cv_sd = np.array(cv_sd)

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # CV_mae
    plt.plot(x_partitioned_data, y_cv, color="#986e42")
    plt.fill_between(x_partitioned_data, y_cv - cv_sd, y_cv + cv_sd, color="#ffe7b5")

    # Model plotted
    plt.plot(x_partitioned_data, y_training, color="Blue")
    plt.fill_between(
        x_partitioned_data,
        y_training - training_sd,
        y_training + training_sd,
        color="#ADD8E6",
    )

    # Remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Add labels and titel
    ax.set_xlabel(
        "Length of the partitioned data",
        size=20,
        fontname="Helvetica",
        fontweight="bold",
    )
    ax.set_ylabel("MAE", size=20, fontname="Helvetica")
    ax.set_title("Learning curve on partitioned data", size=30, fontname="Helvetica")
    ax.legend(
        [
            "Cross-validation mean MAE",
            "Cross-validation standard-deviation",
            "Model performance MAE",
            "Model standard-deviation",
        ],
        loc="upper right",
        shadow=True,
        fontsize="small",
    )

    # Set color
    ax.set_facecolor("white")

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(20, 10)

    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=120, bbox_inches="tight")

    # show
    plt.show()


def bar_plot(
    x: list,
    y: list,
    error_bar: list = None,
    horisontal_line=True,
    save_pdf=True,
    color="white",
    path="",
    title=None,
    x_label=None,
    y_label=None,
) -> None:
    """Plotting a bar_plot .

    Parameters
    ----------
    x : list
        x coordinates
    y : list
        y coordinates
    error_bar : list  (Optional)
        lits of errorbars matching the x coordinates
    horisontal_line : bool
    save_pdf : bool
    path : str
    color : str
        can be matplotlib color names or hex color codes
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    bar_plot"""

    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # Plot
    plt.bar(x, y, edgecolor="black", color=color)  # white

    # Errorbar
    if error_bar != None:
        plt.errorbar(
            x,
            y,
            yerr=error_bar,
            fmt="o",
            color="black",
        )  # ms = 2)

    # add horisontal line
    if horisontal_line:
        plt.axhline(y=100, color="black", linestyle="-")

    # Title and labels
    if title is not None:
        ax.set_title(title, size=30, fontname="Helvetica")
    if x_label is not None:
        ax.set_xlabel(x_label, size=20, fontname="Helvetica")
    if y_label is not None:
        ax.set_ylabel(y_label, size=20, fontname="Helvetica")

    # Set color
    ax.set_facecolor("white")

    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(25, 15)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def horisontal_bar_plot(
    x: list,
    y: list,
    vertical_line: bool = True,
    save_pdf: bool = True,
    path: str = "",
    color="white",
    title=None,
    x_label=None,
    y_label=None,
    size_height: int = 10,
    size_length: int = 8,
) -> None:
    """Plotting a horisontal_bar_plot .

    Parameters
    ----------
    x : list
        x coordinates
    y : list
        y coordinates
    vertical_line : bool
    save_pdf : bool
    path : str
    color : str
    can be matplotlib color names or hex color codes
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    horisontal_bar_plot"""
    #### How can I export a matplotlib figure as a vector graphic with editable text fields?
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Create Figure and Axes instances
    fig, ax = plt.subplots(1)

    # Plot
    plt.barh(x, y, edgecolor="black", color=color)  # white # '#deebf7'

    # Change x labels rotation
    ax.tick_params(rotation=90)

    # remove gridlines
    ax.grid(False)

    # Title and labels
    if title is not None:
        ax.set_title(title, size=30, fontname="Helvetica")
    if x_label is not None:
        ax.set_xlabel(x_label, size=20, fontname="Helvetica")
    if y_label is not None:
        ax.set_ylabel(y_label, size=20, fontname="Helvetica")

    # Background
    ax.set_facecolor("white")

    # add horisontal line
    if vertical_line:
        plt.axvline(x=100, color="black", linestyle="-")

    # rotate sticks
    ax.tick_params(rotation=0)

    ## adding the labels on the bar
    for c in ax.containers:
        ax.bar_label(c, padding=10)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(size_height, size_length)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def correlation_plot(dataframe, x: str, y: str, save_pdf=True, path="") -> None:
    """Plotting a correlation_plot.

    Parameters
    ----------
    dataframe : pd.DataFrame
    x : str
        x coordinates
    y : str
        y coordinates
    save_pdf : bool
    path : str

    Returns
    -------
    correlation_plot"""

    # set seaborn plotting aesthetics as default
    sns.set_context("paper", font_scale=2.0, rc={"lines.linewidth": 1.5})
    dataframe["color"] = "black"
    g = sns.lmplot(
        data=dataframe,
        x=x,
        y=y,
        hue="color",
        palette=["#000000"],
        fit_reg=True,
        height=10,
        line_kws={"color": "black"},
        ci=False,
        legend=False,
    )
    r, p = stats.pearsonr(dataframe[x], dataframe[y])

    ax = plt.gca()

    ax.set_facecolor("white")
    plt.suptitle(
        f"R-squared = {r:.3f} \n  P-value = {p:.3E}",
        y=0.8,
        x=0.3,
        size=15,
        fontname="Helvetica",
        fontweight="bold",
    )
    # SIze matters
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(8, 8)

    if save_pdf and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def bar_plot_w_hue(
    dataframe,
    x: str,
    y: str,
    save_pdf=True,
    path="",
    hue: str = "category",
    title="",
    x_label="",
    y_label="",
    horisontal_line: bool = True,
) -> None:
    """Plotting a correlation_plot.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Dataframe with categories in seperate column from x and y
    x : str
        x coordinates
    y : str
        y coordinates
    save_pdf : bool
    path : str
        path to folder with name of the file
    title : str
    x_label : str
    y_label : str

    Returns
    -------
    bar_plot_w_hue"""

    ax = sns.barplot(x=x, y=y, hue=hue, data=dataframe, palette="pastel")

    ax = plt.gca()
    ax.set_xlabel(title, size=20, fontname="Helvetica")
    ax.set_ylabel(x_label, size=20, fontname="Helvetica")
    ax.set_title(y_label, size=30, fontname="Helvetica")

    # white background
    ax.set_facecolor("white")
    plt.xscale("linear")

    if horisontal_line:
        # normalized line
        ax.axhline(100)

    if save_pdf == True and path != "":
        plt.savefig(path + ".pdf", format="pdf", dpi=300)


def color_dict():
    ### Color codes ###

    yellow = [
        "#a59a00",
        "#aa9e00",
        "#b0a300",
        "#b5a800",
        "#bbad00",
        "#c0b200",
        "#c5b700",
        "#cbbc00",
        "#d0c100",
        "#d6c600",
        "#dbcb03",
        "#e1d011",
        "#e6d51b",
        "#ecda24",
        "#f2df2b",
        "#f7e432",
        "#fde938",
        "#ffef40",
        "#fff647",
    ]
    cpr_color = yellow[::2][::-1]

    orange = [
        "#cd6511",
        "#d26916",
        "#d76d1b",
        "#dc7120",
        "#e17624",
        "#e67a28",
        "#eb7e2d",
        "#f08231",
        "#f58635",
        "#fa8a39",
        "#fe8f3e",
        "#ff9544",
        "#ff9c4b",
        "#ffa351",
        "#ffaa58",
        "#ffb15e",
        "#ffb764",
        "#ffbd6a",
        "#ffc470",
        "#ffca76",
        "#ffd07c",
        "#ffd682",
    ]
    pcpr_color = orange[::6][::-1]

    blue = [
        "#2d89bc",
        "#348dc0",
        "#3b92c5",
        "#4197ca",
        "#479bcf",
        "#4da0d4",
        "#52a5d9",
        "#58aade",
        "#5daee3",
        "#63b3e8",
        "#68b8ee",
        "#6ebdf3",
        "#73c2f8",
        "#78c7fd",
        "#7eccff",
        "#84d1ff",
        "#8ad7ff",
        "#8fdcff",
        "#95e2ff",
        "#9be7ff",
        "#a1edff",
        "#a6f3ff",
        "#acf8ff",
        "#b2feff",
    ]
    g8h_color = blue[::3][::-1]

    green = [
        "#24b161",
        "#2bb565",
        "#32ba69",
        "#38bf6d",
        "#3ec371",
        "#44c876",
        "#4acd7a",
        "#4fd27e",
        "#54d683",
        "#5adb87",
        "#5fe08b",
        "#64e590",
        "#69ea94",
        "#6eee99",
        "#73f39d",
        "#78f8a2",
        "#7efda7",
        "#91ffb9",
    ]
    pg8h_color = green[::5][::-1]

    # List index
    g8h_list_of_index = [
        "CacuG8H",
        "OpumG8H",
        "CroG8H",
        "VminG8H",
        "SmusG8H",
        "RsepG8H",
        "OeuG8H",
        "CcalG8H",
    ]
    cpr_list_of_index = [
        "CroCPR",
        "AanCPR",
        "AraCPR",
        "CloCPR",
        "RseCPR",
        "AhuCPR",
        "AniCPR",
        "CacCPR",
        "OeuCPR",
        "CpoCPR",
    ]
    pg8h_list_of_index = ["CYC1", "ENO2", "PCK1", "RPL15B"]
    pcpr_list_of_index = ["CCW12", "TPI1", "MLS1", "URE2"]
    all_parts = (
        g8h_list_of_index + pg8h_list_of_index + pcpr_list_of_index + cpr_list_of_index
    )

    G8H_color_codes = dict(zip(g8h_list_of_index, g8h_color))
    pG8H_color_codes = dict(zip(pg8h_list_of_index, pg8h_color))
    pCPR_color_codes = dict(zip(pcpr_list_of_index, pcpr_color))
    CPR_color_codes = dict(zip(cpr_list_of_index, cpr_color))
    # making a final dict with colors
    colorDict = dict(
        G8H_color_codes, **pG8H_color_codes, **pCPR_color_codes, **CPR_color_codes
    )

    return colorDict


def plot_phylo_tree(alignment_file, save_pdf=True, path="", height=10, wideness=8):
    """Plotting a phylogenetic tree.

    Parameters
    ----------
    alignment_file : Bio.Align.MultipleSeqAlignment
        A multiple sequence alignemt made into a Bio.Align.MultipleSeqAlignment obejct
        Dataframe with categories in seperate column from x and y
    save_pdf : bool
    path : str
        path to folder with name of the file

    Returns
    -------
    bar_plot_w_hue"""

    # calculate distances in protein identety:
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment_file)

    # Construct tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Plot tree
    # set the size of the figure
    plt.rc("font", family="Helvetica")
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    fig = plt.figure(figsize=(height, wideness), dpi=300)

    plt.rcParams.update({"font.size": 10})
    axes = fig.add_subplot(1, 1, 1)

    fig1 = plt.gcf()
    Phylo.draw(tree, axes=axes, do_show=False, branch_labels=None)
    pylab.axis("off")
    # homologs
    plt.rcParams.update({"font.size": 10})

    # Saving the plot
    if save_pdf and path != "":
        ## save pdf
        plt.savefig(path + ".pdf", format="pdf", dpi=120, bbox_inches="tight")

    plt.show()
