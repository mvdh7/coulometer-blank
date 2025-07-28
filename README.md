# coulometer-blank

## Requirements

To run the code here, you need to install koolstof v1.0.0-b.3 and PyArrow (to open Parquet files).  You can do this with

    pip install -r requirements.txt

## Scripts

read_datasets.py imports and processes the data files (from the data directory).  You don't need to run this script directly, it's imported by the other scripts when required.

print_stats.py prints out a bunch of statistics which are reported in the manuscript.

To recreate the figures in the manuscript:

  * Figure 1 - plot_increments.py
  * Figure 2 - plot_sessions_tidy4.py
  * Figure 3 - plot_histograms.py
  * Figure 4 - plot_uncertainty.py
  * Supp. Fig. S1 - plot_nuts_dic.py
  * Other supp. figs. - plot_sessions_all.py

The figures are saved in .png format in the figures directory.
