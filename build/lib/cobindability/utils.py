import logging
import pandas as pd
import numpy as np
from scipy.stats import zscore


def config_log(switch, logfile=None):
    """
    Configureing the logging module.

    Parameters
    ----------
    switch : bool
        Debugging switch.
    Returns
    -------
    None.

    """
    if switch is True:
        if logfile is None:
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
        else:
            logging.basicConfig(
                filename=logfile,
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
    else:
        if logfile is None:
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
        else:
            logging.basicConfig(
                filename=logfile,
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)


def cal_zscores(infile, outfile):
    """calculate z-score of the six collocation measurements
    TF_name C       J       SD      SS      PMI     NPMI
    RAD21   0.1446  0.0224  0.0438  0.9326  2.0074  0.3417
    SMC3    0.143   0.0214  0.042   0.9525  2.0285  0.3428
    SMC1A   0.1413  0.0211  0.0413  0.9462  2.0219  0.3407
    TRIM22  0.14    0.0214  0.0419  0.9127  1.9858  0.3355
    STAG1   0.1368  0.0191  0.0375  0.9787  2.0556  0.3407
    """
    logging.info("Calculate Z-scores from \"%s\"" % infile)
    df = pd.read_csv(infile, index_col=0, sep="\t", engine='python')
    print(df)
    col_names = ['C', 'J', 'SD', 'SS', 'PMI', 'NPMI']
    input_names = df.columns
    if(all(x in input_names for x in col_names)):
        numeric_cols = col_names
    else:
        numeric_cols = df.select_dtypes(include=[np.number]).columns
    df2 = df[numeric_cols].apply(zscore)
    z_scores = df2.sum(axis=1, numeric_only=True)/(len(numeric_cols)**0.5)
    df2['Zscore'] = z_scores
    logging.info("Save Z-scores to \"%s\"" % outfile)
    df2.to_csv(outfile, sep="\t")
    print(df2)
