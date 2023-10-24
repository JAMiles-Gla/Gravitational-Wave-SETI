# local imports
from load_files import dataLoader
from likelihood_calculations import likelihood
from plotting_functions import plotStatHist
from sequence_functions import getPrimes

# TODO: make this file less messy, clean up function calls, correct data segments

# foreground and background data files
fgfile = "GW_data/wave_O3_K99_C01_LH_BurstLF_BKG_run1_M2_V_hvetoLH_foreground.csv"
bgfile = "GW_data/wave_O3_K99_C01_LH_BurstLF_BKG_run1_M2_V_hvetoLH_background.csv"

# load data
fgslices, bgslices = dataLoader(fgfile, bgfile)

# test parameters
maxSeq = 5
distanceWindow = 100 # degrees
timeWindow = 500 # seconds, around 8 minutes

# initialise storage lists
logLikelihoods = []
maxLogLikelihoods = []

# loop over data segments
for i in range(len(bgslices)):

    # print progress
    print('Running on data segment {0}...'.format(i+1))

    # call algorithm
    L, maxL = likelihood(bgslices[i], distanceWindow, timeWindow, getPrimes, maxSeq, plot=False, verbose=False)

    # save to lists
    logLikelihoods.append(L)
    maxLogLikelihoods.append(maxL)

# foreground test run
fgL, fgMaxL = likelihood(fgslices[0], distanceWindow, timeWindow, getPrimes, maxSeq, plot=False, verbose=False)

# plot histogram
plotStatHist(maxLogLikelihoods, fgMaxL)
print('Absolute maximum background:', max(maxLogLikelihoods))
print('Foreground:', fgMaxL)
