import numpy as np
import matplotlib.pyplot as plt

# for LaTeX font on plots
plt.rcParams.update({'font.size': 15, "font.family": "Times New Roman", "text.usetex": True})


def histListLengths(timeLocs):
  '''
  Plots a histogram of the lengths of each sequence list.

  Params:  timeLocs:  list of time sequences
  '''

  # find lengths of each sublist
  listLengths = []
  for l in timeLocs:
    listLengths.append(len(l))

  # show max list length
  print('Max list length: ', max(listLengths))

  # initialise figure
  plt.figure(figsize=(8, 6))

  # define histogram bins in log space
  logbins = np.geomspace(min(listLengths), max(listLengths), 20)

  # plot
  plt.hist(listLengths, bins=15, color='k')
  plt.xscale('log')
  plt.yscale('log')

  # plot configs
  plt.xlabel('List length', fontsize=20)
  plt.ylabel('Number of lists', fontsize=20)
  plt.show()


def plotTimeLocs(trigger1Time, trigger2Time, timeLocs, sequenceLocs, timeWindows):
  '''
  Plots a graph of the time locations in each sequence list.

  Params:  trigger1:      time of first trigger
           trigger2:      time of second trigger
           timeLocs:      list of time sequences
           sequenceLocs:  sequence numbers for annotation
           timeWindows:   list of uncertainty windows around each time location
  '''

  # initialise figure
  plt.figure(figsize=(8, 6))

  # loop over outer list dimension
  for i in range(1, len(timeLocs)+1):

    # set y axis to list index
    index = i * np.ones(len(timeLocs[i-1]))

    # scatter time locations
    plt.scatter(trigger1Time, i, s=5, c='b')
    plt.scatter(trigger2Time, i, s=5, c='b')
    #plt.errorbar(timeLocs[i-1], index, xerr=timeWindows[i-1], ls='none', capsize=3, mew=0.8, lw=0.8, c='k')
    plt.scatter(timeLocs[i-1], index, s=5, c='r')

    # create full list of sorted times
    fullTimes = [trigger1Time, trigger2Time]
    fullTimes.extend(timeLocs[i-1])
    fullTimes.sort()

    # test prints
    #print(fullTimes)
    #print(sequenceLocs[i-1])

    # annotate each trigger by sequence numbers
    for j in range(len(fullTimes)):
      plt.annotate(str(sequenceLocs[i-1][j]), (fullTimes[j], i+0.1), fontsize=10)

  # plot configs
  plt.xlabel('Time location [GPS]', fontsize=20)
  plt.ylabel('Sequence', fontsize=20)
  plt.show()


def plotStatHist(maxLogLikelihoods, fgLogLikelihood):

  # initialise figure
  plt.figure(figsize=(8, 6))

  # plot
  plt.hist(maxLogLikelihoods, bins=35, color='k') # density=True for normalised version
  plt.hist(fgLogLikelihood, bins=1, color='r')

  # plot configs
  plt.xlabel('Log likelihood', fontsize=20)
  plt.ylabel('Number of events', fontsize=20)
  plt.show()