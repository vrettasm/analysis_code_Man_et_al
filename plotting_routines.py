""" Plotting Routines

This script includes the basic plotting functions of the video analysis results.
Since the analysis results are stored in comma separated values (csv) and Java
Script Object Notation (json) files, the packages that handle these files (i.e.
pandas and json) are assumed to have been present in the python installation.

The packages includes three public functions:

    1) fit_plot
    2) log_plot
    3) hist_plot

and one private (d_exp_func). Since the curve fitting is fast we do not store
the data to separate files, but rather re-run the fit before plotting. This
package is included from the make_plots script.

Note:
    The curve fitting algorithm "scipy.optimize.curve_fit" might throw some
    warning messages on the output (screen), during execution. These can be
    ignored for the current input data files.
"""

# Load the required packages.
import json
import numpy as np
from pathlib import Path
from matplotlib import gridspec
import matplotlib.pyplot as plt

# Public interface.
__all__ = ["fit_plot", "hist_plot", "log_plot"]

# INFO:
__author__ = "Michalis Vrettas, PhD"
__email__ = "michail.vrettas@gmail.com"

# Define helper function.
def d_exp_func(t, a1, b1, a2, b2, c0):
    """
        Description:
            Double exponential function. It is used to fit the probability data
            (autocorrelation func.) that are generated from the video analysis.

        Args:
             t (float): independent variable representing time.

            a1 (float): scale factor of the first exponential function.
            b1 (float): decaying factor of the first exponential.

            a2 (float): scale factor of the second exponential function.
            b2 (float): decaying factor of the second exponential.

            c0 (float): offset (bias) parameter.

        Returns:
            f(t) = a1*e^(-b1*t) + a2*e^(-b2*t) + c0.

        Note:
            Since this is only a helper function, we do not have to include it
            in the public interface "__all__".
    """
    return a1*np.exp(-b1*t) + a2*np.exp(-b2*t) + c0
# _end_def_


def fit_plot(filename=None):
    """
        Description:
            Fit curve plotting function. It uses the autocorrelation data from
            the video analysis, to fit a double exponential function and plots
            the fit vs the real data, along with the error between the curves.

            The fitted parameters, are given in the legend information.

        Args:
             filename (Path): path object for the json file that contains the
             probability data.

        Returns:
            Saves the plot to the same directory as the input file.

        Raises:
            ValueError if the input filename:
                1) is not given
                2) does not exist
                3) is not json file
    """
    # Load the curve_fit locally.
    from scipy.optimize import curve_fit

    # Check for empty input.
    if filename is None:
        raise ValueError("fit_plot: No input.")
    # _end_if_

    # Check for wrong input.
    if not filename.exists():
        raise ValueError("fit_plot: File does not exist.")
    # _end_if_

    # Check for wrong file type.
    if filename.suffix != ".json":
        raise ValueError("fit_plot: Filename is not a json file.")
    # _end_if_

    try:
        # Load the results in dictionary format.
        with open(filename, 'r') as json_file:
            data = json.load(json_file)
        # _end_with_

        # Plot all directory entries in one plot.
        for n, k in enumerate(data.keys()):

            # Max limit in the plot.
            max_limit = 0

            # First export the data from the dictionary.
            x_k = np.asarray(list(data[k].keys()), dtype=np.float32)
            y_k = np.asarray(list(data[k].values()), dtype=np.float32)

            try:
                # Do the non-linear curve fitting.
                p_opt, p_cov = curve_fit(d_exp_func, x_k, y_k, method="trf")
            except RuntimeError as e0:
                print(e0, end="\n")
                continue
            # _end_try_

            # Estimate the error (std) of the parameters.
            p_err = np.sqrt(np.diag(p_cov))

            # Calculate the function with the new (estimated) parameters.
            y_plot = d_exp_func(x_k, *p_opt)

            # Compute the root-mean-square-deviation.
            rmsd = np.sqrt(np.sum((y_k-y_plot)**2) / len(y_k))

            # The last key entry should have the maximum value.
            max_limit = np.maximum(max_limit, int(list(data[k].keys())[-1]))

            # Make a new file.
            plt.figure(n, figsize=(8, 6))
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])

            # First plot.
            ax0 = plt.subplot(gs[0])

            # Plot the data vs the fitting.
            ax0 = plt.plot(x_k, y_k, 'k-', label='data')
            ax0 = plt.plot(x_k, y_plot, 'r--',
                           label="fit (value $\pm$ error)\n"
                           "a1={0:.2f} $\pm$ {1:.2f}\n"
                           "b1={2:.2f} $\pm$ {3:.2f}\n"
                           "a2={4:.2f} $\pm$ {5:.2f}\n"
                           "b2={6:.2f} $\pm$ {7:.2f}\n"
                           "c0={8:.2f} $\pm$ {9:.2f}".format(p_opt[0], p_err[0],
                                                             p_opt[1], p_err[1],
                                                             p_opt[2], p_err[2],
                                                             p_opt[3], p_err[3],
                                                             p_opt[4], p_err[4]))
            # Setup title and labels.
            plt.title("{0} - RMSD = {1:.3f}".format(k, rmsd))
            plt.ylabel(r"P($\tau$) vs $\hatP(\tau)$")

            # Set the x-ticks with the correct dt.
            plt.xticks(np.linspace(0, max_limit, 6),
                       np.round(np.linspace(0, max_limit, 6)*0.04, 3))

            # The legend contains the fitting parameters (+/-) error.
            plt.legend()
            plt.grid(True)

            # Second plot.
            ax1 = plt.subplot(gs[1])
            ax1.plot(y_k-y_plot, '.k')
            plt.grid(True)
            plt.xticks(np.linspace(0, max_limit, 6),
                       np.round(np.linspace(0, max_limit, 6)*0.04, 3))
            plt.ylabel(r"$\Delta P(\tau)$")
            plt.xlabel(r"$\tau$ (sec)")

            # Save the file to the same directory.
            plot_name = filename.parent / Path("fit_" + k + ".png")
            plt.savefig(plot_name, dpi=300)

            # Display a message.
            print(" Create:", plot_name)

            # Close the figure.
            plt.close(n)
        # _end_for_

    except IOError as e1:
        print(e1)
    # _end_try_
# _end_def_


def hist_plot(filename=None):
    """
        Description:
            Histogram plotting function. It uses the detected spots data from
            the video analysis, to plot a histogram of the number of detected
            spots from all frames. In addition it fits a KDE on the same data
            and contrasts the PDF over the histogram.

        Args:
             filename (Path): path object for the csv file that contains the
             number of detected spots for each frame.

        Returns:
            Saves the plot to the same directory as the input file.

        Raises:
            ValueError if the input filename:
                1) is not given
                2) does not exist
                3) is not csv file
    """
    # Local imports.
    import pandas as pd
    from scipy.stats import iqr
    from sklearn.neighbors import KernelDensity

    # Check for empty input.
    if filename is None:
        raise ValueError("hist_plot: No input.")
    # _end_if_

    # Check for wrong input.
    if not filename.exists():
        raise ValueError("hist_plot: File does not exist.")
    # _end_if_

    # Check for wrong file type.
    if filename.suffix != ".csv":
        raise ValueError("hist_plot: Filename is not a csv file.")
    # _end_if_

    # Set the plot style.
    plt.style.use('ggplot')

    try:
        # Start a new figure.
        plt.figure(0)

        # Load the data from the .csv file.
        all_spots = pd.read_csv(filename)

        # Copy only the relevant columns.
        df_spots = all_spots[[ "FRAME"]].copy()

        # Sort the row by frame.
        df_spots.sort_values(by=["FRAME"], inplace=True)

        # Count the (detected) spots per frame.
        counts_by_frame = df_spots.groupby('FRAME').size()

        # Get the mean and std of the values.
        mean_std = " Mean:{0:.3f}\n Std:{1:.3f}".format(counts_by_frame.values.mean(),
                                                        counts_by_frame.values.std())
        # Freedman–Diaconis rule.
        bin_width = 2.0*iqr(counts_by_frame.values, axis=0,
        interpolation='midpoint')/(len(counts_by_frame)**(1.0/3.0))

        # Number of bins.
        n_bins = int((counts_by_frame.values.max()-counts_by_frame.values.min()/bin_width))

        # Create the histogram.
        plt.hist(counts_by_frame.values, bins=n_bins, density=True,
                 facecolor='tab:blue', alpha=0.5, label=mean_std)

        # Setup title and labels.
        plt.title(filename.stem)
        plt.xlabel("Number of detected spots")
        plt.ylabel("Density")

        plt.grid(True)
        plt.legend()

        # instantiate and fit the KDE model
        kde = KernelDensity(bandwidth=1.0, kernel='gaussian')
        kde.fit(counts_by_frame.values[:, None])

        # Plot x-range.
        x_d = np.linspace(0.0, counts_by_frame.values.max()+5, 1000)

        # score_samples returns the log of the probability density
        logprob = kde.score_samples(x_d[:, None])

        plt.plot(x_d, np.exp(logprob), color="black", alpha=0.95)

        plot_name = filename.parent / Path("hist_" + filename.stem + ".png")

        plt.savefig(plot_name, dpi=300)

        # Display a message.
        print(" Create:", plot_name)

        plt.close(0)
    except IOError as e:
        print(e)
    # _end_try_
# _end_def_


def log_plot(filename=None):
    """
        Description:
            Semi log-scale plotting function. It uses the autocorrelation data
            (probabilities) from the video analysis and plots the results on
            a logarithmic scale to contrast the control vs the rest cases.

        Args:
             filename (Path): path object for the json file that contains the
             number of detected spots for each frame.

        Returns:
            Saves the plot to the same directory as the input file.

        Raises:
            ValueError if the input filename:
                1) is not given
                2) does not exist
                3) is not json file
    """

    # Check for empty input.
    if filename is None:
        raise ValueError("log_plot: No input.")
    # _end_if_

    # Check for wrong input.
    if not filename.exists():
        raise ValueError("log_plot: File does not exist.")
    # _end_if_

    # Check for wrong file type.
    if filename.suffix != ".json":
        raise ValueError("log_plot: Filename is not a json file.")
    # _end_if_

    try:
        # Load the results in dictionary format.
        with open(filename, 'r') as json_file:
            data = json.load(json_file)
        # _end_with_

        # Make a new file.
        plt.figure(0)

        # Max limit in the plot.
        max_limit = 0

        # Plot all directory entries in one plot.
        for key in data:
            # Check if there is a "control" case.
            if "Control" in key:
                plt.semilogy(list(data[key].values()), 'k--', label=key)
            else:
                plt.semilogy(list(data[key].values()), label=key)
            # _end_if_

            # The last key entry should have the maximum value.
            max_limit = np.maximum(max_limit, int(list(data[key].keys())[-1]))
        # _end_for_

        # Final adjustments
        plt.grid(True)
        plt.xlabel(r"$\tau$ (sec)")
        plt.ylabel(r"$P(\tau)$")
        plt.xticks(np.linspace(0, max_limit, 6),
                   np.round(np.linspace(0, max_limit, 6)*0.04, 3))
        plt.legend()

        # Save the file to the same directory.
        plot_name = filename.parent / Path("log_" + filename.stem + ".png")
        plt.savefig(plot_name, dpi=300)

        # Display a message.
        print(" Create:", plot_name)

        # Close the figure.
        plt.close(0)

    except IOError as e:
        print(e)
    # _end_try_
# _end_def_
