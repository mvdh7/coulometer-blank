# coulometer-blank

This code accompanies a paper undergoing peer review.  The preprint can be found here:

> Humphreys, M. P., and Ossebaar, S.: Blank variability in coulometric measurements of dissolved inorganic carbon, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-3644, 2025. 

## Requirements

To run the code here, you need to install koolstof v1.0.0-b.3 and PyArrow (to open Parquet files).  You can do this with

    pip install -r requirements.txt

using the [requirements.txt](requirements.txt) file in the repo.

## Scripts

[read_datasets.py](read_datasets.py) imports and processes the data files (from the [data](data) directory).  You don't need to run this script directly, it's imported by the other scripts when required.

[print_stats.py](print_stats.py) prints out a bunch of statistics which are reported in the manuscript.

To recreate the figures in the manuscript:

  * Figure 1 - [plot_increments.py](plot_increments.py)
  * Figure 2 - [plot_sessions_tidy4.py](plot_sessions_tidy4.py)
  * Figure 3 - [plot_histograms.py](plot_histograms.py)
  * Figure 4 - [plot_uncertainty.py](plot_uncertainty.py)
  * Supp. Fig. S1 - [plot_nuts_dic.py](plot_nuts_dic.py)
  * Other supp. figs. - [plot_sessions_all.py](plot_sessions_all.py)

The figures are saved in .png format in the [figures](figures) directory.
