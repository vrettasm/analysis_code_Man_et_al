#!/usr/bin/env python3

""" Make Probabilities

This script is called to calculate the probabilities (or the autocorrelation)
of the detected spots appearing in times 't' and 't+tau'. The algorithm uses
the csv generated files from the ImageJ analysis that include, among others,
the (x,y) positions and the frame (f) for all spots on all frames, and counts
the total number of spots that appear both at frames f(t) and f(t+tau).

The implementation is done using parallel execution, thus operates on multiple
files simultaneously. The results from all parallel executions are grouped in
a dictionary, where each processed filename is used as a key for the data.

The program, except from the directories that uses as input, accepts also a few
other parameters for the autocorrelation function. These are:

    1) --num_frames (or the 'tau')  (maximum number of frames to look ahead of
    the current time 't')
    2) --radius (in pixel units to look for similarity between the coordinates
    x and y for the spots)

These options can be given in any order and operate on any number of input
directories, assuming that they contain the right file types.

Examples:

(1) $ python3 make_probabilities.py --dir ../data/Mutants/A30P/SV1/

    (will run the computation on the data files that are included in
    the "../data/Mutants/A30P/SV1/" directory only. Default values for
    num_frames and radius are used.)

(2) $ python3 make_probabilities.py --dir ../data/Mutants/A30P/SV1/
                                          ../data/Mutants/A30P/SV2/
                                          ../data/Mutants/A30P/SV3/

    (will run the computation on the data files that are included in
    the "../data/Mutants/A30P/SV1/", "../data/Mutants/A30P/SV2/" and
    "../data/Mutants/A30P/SV3/" directories. Default values for num_frames
    and radius are used.)

(3) $ python3 make_probabilities.py --dir ../data/Mutants/A30P/SV1/
                                    --num_frames 500 --radius 1.5

    (will run all the computation on the data files that are included in
    the "../data/Mutants/A30P/SV1/" directory only. Values for num_frames
    and radius are given explicitly.)

(4) etc.

"""

# Load the required packages.
import os
import sys
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path
import concurrent.futures

# INFO:
__author__ = "Michalis Vrettas, PhD"
__email__ = "michail.vrettas@gmail.com"

# Define a function to calculate the probability.
def calc_probability(filename, T=100, radius=1.0):
    """
        Description:
            Calculates the autocorrelation function (or the probability) for a
            given filename and input parameters.

        Args:
            filename (Path): path object for the csv file that contains the
            detected spot data (x,y,f).

            T (int): maximum number of frames to look ahead in the autocor-
            relation function. Default value is 100 (frames).

            radius (float): in pixel units. Default value is 1.0 (pixel).

        Returns:
            A tuple (filename, prob_tau), where "prob_tau" is the dictionary
            with the estimated autocorrelation values for each "tau".
    """
    # Display which file is being processed.
    print(" Processing file: {0}".format(filename.stem))

    # Load the data from the .csv file.
    all_spots = pd.read_csv(filename)

    # Copy only the relevant columns.
    df_spots = all_spots[[ "FRAME", "POSITION_X", "POSITION_Y"]].copy()

    # Sort the row by frame.
    df_spots.sort_values(by=["FRAME", 'POSITION_X', "POSITION_Y"], inplace=True)

    # Get the maximum number of frames.
    num_frames = df_spots["FRAME"].max()+1

    # Holds the spots (x, y) for each frame.
    frame = {}

    # Extract the spots (x, y) in numpy array.
    for i in range(0, num_frames):
        frame[i] = df_spots[df_spots["FRAME"] == i][[ "POSITION_X", "POSITION_Y"]].to_numpy()
    # _end_for_

    # Declare a dictionary to hold the probability results.
    prob_tau = {}

    # Start the timer.
    t0 = time.process_time()

    # Repeat for all tau values.
    for tau in range(0, T+1):

        # Counts the spots that satisfy the condition.
        count_spots = 0

        # Counts all the spots.
        total_spots = 0

        # Scan all the frames dataframe.
        for t in range(0, num_frames-tau):

            # Set the "future" frame.
            f = t + tau

            # Number of spots in the frame.
            num_spots_t = frame[t].shape[0]

            # If there are no spots detected, continue to the next frame.
            if num_spots_t == 0:
                continue
            # _end_if_

            # Add the current spots in the total counter.
            total_spots += num_spots_t

            # Special case: "tau == 0".
            if tau == 0:
                count_spots += num_spots_t
            else:

                # Count how many spots from frame "t" also exist in frame "t + tau".
                for x_i, y_i in frame[t]:

                    # If any of the points in the spots_f list exists inside
                    # the circle (x_i, y_i, r), then we increase the counter.
                    if np.any(np.sqrt((x_i - frame[f][:,0])**2 +
                                      (y_i - frame[f][:,1])**2) <= radius):
                        count_spots += 1
                    # _end_if_
                # _end_for_
            # _end_if_
        # _end_for_

        # Add the total counts (for a given tau).
        prob_tau[tau] = count_spots/total_spots
    # _end_for_

    # Stop timer.
    t_elapsed = time.process_time() - t0

    # Display a final message.
    print(" File: {0} completed in {1:.2f} seconds.".format(filename.stem, t_elapsed))

    # Return the probability: P(tau)
    return (filename.stem, prob_tau)
# _end_def_

# Define the main calculation function.
def main(data_dir=None, T=100, radius=0.5):
    """
        Description:
            Main function that handles the input (list) directory along with
            the rest of the parameters (T and radius).

        Args:
             data_dir (list): list of directories that include the input files
             to generate the requested plots.

             T (int): maximum number of frames to look ahead in the autocor-
             relation function. Default value is 100 (frames).

             radius (float): in pixel units. Default value is 0.5 (pixel).
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

        # Only the csv files contain the spot per frame.
        fnames = [f for f in dk.iterdir() if f.name.endswith(".csv")]

        # If there are not files exit.
        if not fnames:
            sys.exit(" Error: No input files were found.")
        # _end_if_

        # Stores the probabilities from all analysis.
        prob_all = {}

        # Analysis in parallel.
        with concurrent.futures.ProcessPoolExecutor() as ppe:

            # Run the analysis in parallel.
            results = [ppe.submit(calc_probability, fname, T=T, radius=radius)
                        for fname in fnames]

            # Display info.
            print("\n ... waiting for the process pool to finish ...\n")

            # Get the results.
            for proc in concurrent.futures.as_completed(results):
                # Key: is the filename that we analyzed
                # Value: is the estimated probability.
                key, value = proc.result()

                # Store them in a separate dictionary.
                prob_all[key] = value
            # _end_for_
        # _end_with_

        # Construct a filename.
        prob_name = "prob_" + dk.parent.name + "_" + dk.name + ".json"

        # Save the results to json file format.
        with open((dk/prob_name), 'w') as fp:
            json.dump(prob_all, fp, indent=4)
        # _end_with_

    # _end_for_

# _end_def_


# Run the script.
if __name__ == "__main__":

    # Check if we have given input parameters.
    if len(sys.argv) > 1:
        # Local import.
        import argparse

        # Create a parser object
        parser = argparse.ArgumentParser()

        # Number of frames (tau)
        parser.add_argument("--num_frames", help="number of frames for the"
                            " autocorrelation function (probability)",
                            default=625)
        # Search radius.
        parser.add_argument("--radius", help="radius of the circle around"
                            " the spot to look for similarity",
                            default=0.5)

        # List of directories to analyze.
        parser.add_argument("--dir", nargs = '*', help = " list of directories"
        " with the csv files that contain the (x,y,f) from the detected spots.")

        # Parse the arguments.
        args = parser.parse_args()

        # Call the main function.
        main(args.dir, args.num_frames, args.radius)

        # Display final info.
        print(" Done.")
    else:
        sys.exit("Error: Not enough input parameters.")
    # _end_if_

# _end_program_
