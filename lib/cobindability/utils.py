import logging


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
