#!/usr/bin/env python3

""" Make Plots

This script is called to generate the plots of the video analysis. It provides
three options:

    1) --log  (plots the results on a logarithmic scale)
    2) --fit  (plots the curve fitting vs the real data)
    3) --hist (plots the histogram of the number of detected spots per frame)

These options can be given in any order and operate on any number of input
directories, assuming that they contain the right file types (csv and json).

Examples:

(1) $ python3 make_plots.py --log --fit --hist ../data/GM/SV1/

    (will run all three plot functions on the data files that are included
    in the "../data/GM/SV1/" directory)

(2) $ python3 make_plots.py --log ../data/GM/SV1/ ../data/GM/SV2/

    (will run the 'log' plot function on the data files that are included
    in the "../data/GM/SV1/" and "../data/GM/SV2/" directories)

(3) $ python3 make_plots.py --fit --log ../data/GM/SV1/ ../data/GM/SV2/

    (will run the 'fit' and 'log' plot functions on the data files that are
    included in the "../data/GM/SV1/" and "../data/GM/SV2/" directories)

(4) etc.

"""

# Load the required packages.
import os
import sys
from pathlib import Path
import plotting_routines as plr

# INFO:
__author__ = "Michalis Vrettas, PhD"
__email__ = "michail.vrettas@gmail.com"

# Define the main function.
def main(data_dir=None, kind=None):
    """
        Description:
            Main function that handles the input (list) directory along with
            the kind of plots we want to generate and calls iteratively the
            right function(s).

        Args:
             data_dir (list): list of directories that include the input files
             to generate the requested plots.

             kind (str): available options include ('log', 'fit', 'hist').
    """

    # Check if an input directory is given.
    if data_dir is None:
        sys.exit(" Error: No input parameters.")
    # _end_if_

    # Scan all given directories.
    for dk in data_dir:

        # Convert string to path object.
        dk = Path(dk)

        # Check the validity of the input parameter.
        if not dk.is_dir():
            sys.exit(" Error: Input '{0}' is not a valid directory.".format(dk))
        # _end_if_

        if kind == "hist":
            # Only the csv files contain the spot counts.
            fnames = [f for f in dk.iterdir() if f.name.endswith(".csv")]
        else:
            # Only the json files contain the probabilities.
            fnames = [f for f in dk.iterdir() if f.name.endswith(".json")]
        # _end_if_

        # If there are not files exit.
        if not fnames:
            sys.exit(" Error: No input files were found.")
        # _end_if_

        # Process all json files in the directory.
        for fname in fnames:

            # Call the function according to the requested type.
            try:
                if kind == "hist":
                    plr.hist_plot(fname)
                elif kind == "fit":
                    plr.fit_plot(fname)
                elif kind == "log":
                    plr.log_plot(fname)
                else:
                    pass
            except RuntimeError as e:
                print(e)
                continue
            # _end_try_
        # _end_for_

    # _end_for_
# _end_def_


# Run the script.
if __name__ == "__main__":

    # Check if we have given input parameters.
    if len(sys.argv) > 1:
        # Get the list of all arguments.
        args_list = sys.argv[1:]

        # Check if we want help.
        if ("-h" in args_list) or ("--help" in args_list):
            print("\n usage: make_plots.py [--log] [--fit] [--hist] [ directories ]\n")
            print(" Optional arguments : ")
            print("                    : [--log]  creates the log plots ")
            print("                    : [--fit]  creates the curve fitting plots ")
            print("                    : [--hist] creates the histogram plots")
            print("\n make_plots.py [--help, -h] shows this message and exits.\n")
            print("\n Description : ")
            print(" It accepts a list of directories, with the required analysis \n"
                  " 'csv' and 'json' files, and creates the plots. The options   \n"
                  " can be combined in any order to produce the plots. Note, that\n"
                  " if the plots exist they will be overwritten without warning.")
            sys.exit(0)
        # _end_if_

        # Initialize options to false.
        make_log = False
        make_fit = False
        make_hist = False

        # Check which options are given.
        if "--log" in args_list:
            make_log = True
            args_list.remove("--log")
        # _end_if_

        if "--fit" in args_list:
            make_fit = True
            args_list.remove("--fit")
        # _end_if_

        if "--hist" in args_list:
            make_hist = True
            args_list.remove("--hist")
        # _end_if_

        if not args_list:
            sys.exit("Error: Directory list is empty.")
        # _end_if_

        # Now run all the options.
        if make_log:
            print(" >> Making log plots: ")
            main(args_list, kind="log")
            print(" ")
        # _end_if_

        if make_fit:
            print(" >> Making fit plots: ")
            main(args_list, kind="fit")
            print(" ")
        # _end_if_

        if make_hist:
            print(" >> Making hist plots: ")
            main(args_list, kind="hist")
            print(" ")
        # _end_if_

        # Final message.
        print(" Done! ")
    else:
        sys.exit("Error: Not enough input parameters.")
    # _end_if_
# _end_program_
