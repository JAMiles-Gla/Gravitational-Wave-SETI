import numpy as np
import pandas as pd

# TODO: this is a temporary version, needs to be adapted when full search is performed


def dataLoader(fgfile, bgfile):
    # load background data
    bgdata = pd.read_csv(bgfile)

    # get first 100 * 30 triggers for testing
    bgdata = bgdata.head(15000)

    # split dataframe into 100 equal parts
    bgslices = np.array_split(bgdata, 500)

    # find only data chunks that are sufficiently long
    # TODO: change this so that times are equal instead of taking equal numbers of triggers
    trueBgSlices = []
    for sl in bgslices:
      totalTime = (sl['time0'].max() - sl['time0'].min()) / (3600 * 24)
      if 15 < totalTime < 35:
        trueBgSlices.append(sl)

    # load foreground data
    fgdata = pd.read_csv(fgfile)

    # pick first 30 triggers to test on
    fgdata = fgdata.head(120)
    fgslices = np.array_split(fgdata, 4)

    return fgslices, trueBgSlices
