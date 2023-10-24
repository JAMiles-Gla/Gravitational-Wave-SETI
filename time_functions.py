import numpy as np

# local imports
from uncertainty_windows import windows


def timeLocations(trigger1Time, trigger2Time, minTime, maxTime, sequence, maxSeq, minDelta, timeWindow, plot=False, verbose=False):
  '''
  Finds locations of time points that need to be checked.

  Params:     trigger1Time:   time of first trigger in pair
              trigger2Time:   time of second trigger in pair
              minTime:        minimum time in data, used to find minimum time
                              to be checked in the sequence
              maxTime:        maximum time in data, used to find maximum time
                              to be checked in the sequence
              sequence:       sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:         maximum number of times in sequence, after which
                              the sequence restarts
              minDelta:       minimum value of the time unit to assure
                              reasonable lengths of each trigger sequence
              timeWindow:     uncertainty window around times

  Optional:   plot:             boolean, determines if lists necessary for plots
                                need to be filled, False by default

  Returns:    timeLocs:       two dimensional list containing time locations
                              to be checked for
              timeWindows:    two dimensional list containing uncertainty window
                              sizes for each possible time
              sequenceLocs:   two dimensional list containing the numbers in
                              sequence that each trigger corresponds to
  '''

  # time step between triggers
  step = trigger2Time - trigger1Time

  # initialise list of time sequences and windows
  timeLocs = []
  timeWindows = []
  if plot == True:
    sequenceLocs = []

  # loop over possible sequence steps
  for i in range(1, maxSeq+1):
    for j in range(i, maxSeq+1):

      # find subsequence joining the two triggers
      seq = sequence(i, j)

      # define delta as one sequence unit (e.g., the number 2 in a sequence
      # would be marked as 2 * delta)
      delta = step / sum(seq)

      # check if delta is within allowed bounds
      if delta > minDelta:

        # generate list of time values forwards in sequence
        forwardSequence, forwardLocs = locForward(delta, j, sequence, maxSeq,
                                        trigger2Time, maxTime, plot=plot)

        # generate list of time values backwards in sequence
        backwardSequence, backwardLocs = locBackward(delta, i, sequence, maxSeq,
                                          minTime, trigger1Time, plot=plot)

        # generate list of time values between two initial triggers
        midSequence, midLocs = locBetween(delta, i, j, sequence, maxSeq,
                                trigger1Time, trigger2Time, plot=plot)

        # concatenate lists
        fullSequence = backwardSequence + midSequence + forwardSequence
        if plot == True:
          fullLocs = backwardLocs + midLocs + forwardLocs

        # calculate uncertainty windows for the full sequence
        fullWindow = windows(trigger1Time, trigger2Time, fullSequence, timeWindow, delta)

        # add sequence and window to list of sequence lists to be tested
        if len(fullSequence) != 0:
          timeLocs.append(fullSequence)
          timeWindows.append(fullWindow)
          if plot == True:
            sequenceLocs.append(fullLocs)

      # sanity check for smaller deltas
      else:
        if verbose == True:
          print('Delta too small:', delta)

  # add sequence location list to returns if needed for plotting
  if plot == True:
    return timeLocs, timeWindows, sequenceLocs

  # return two dimensional list of possible positions to check for and their
  # associated uncertainties
  return timeLocs, timeWindows, list()


def locForward(delta, loc2, sequence, maxSeq, trigger2Time, maxTime, plot=False):
  '''
  Finds locations of time points forward in sequence that need to be checked.

  Params:     delta:            sequence unit, (e.g., the number 2 in a sequence
                                would be marked as 2 * delta)
              loc2:             sequence location of second trigger
              sequence:         sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:           maximum number of times in sequence, after which
                                the sequence restarts
              trigger2Time:     time of second trigger in pair
              maxTime:          maximum time in data, used to find maximum time
                                to be checked in the sequence

  Optional:   plot:             boolean, determines if lists necessary for plots
                                need to be filled, False by default

  Returns:    forwardSequence:  list containing time locations to be checked for
              forwardLocs:      list containing the numbers in sequence that
                                each forward trigger corresponds to
  '''

  # initialise sequence list
  forwardSequence = []

  # initialise sequence location list needed for plotting
  if plot == True:
    forwardLocs = []

  # initialise time location
  time = trigger2Time

  # initialise sequence location
  seqLoc = loc2 + 1

  # define exit condition
  while time < maxTime:

    # reset sequence if it reaches the max step
    if seqLoc > maxSeq:
      seqLoc = 1

    # add to time until maximum is reached
    time += delta * sequence(seqLoc, seqLoc)[0]

    # add time location to list, if less than max time
    if time < maxTime:
      forwardSequence.append(time)

      # add to sequence location list
      if plot == True:
        forwardLocs.append(sequence(seqLoc, seqLoc)[0])

    # add final sequence location to list, if already greater than max time
    if plot == True and time > maxTime:
      forwardLocs.append(sequence(seqLoc, seqLoc)[0])

    # update sequence location
    seqLoc += 1

  # include sequence locations in return if needed for plotting
  if plot == True:
    return forwardSequence, forwardLocs

  # return full list
  return forwardSequence, list()


def locBackward(delta, loc1, sequence, maxSeq, minTime, trigger1Time, plot=False):
  '''
  Finds locations of time points forward in sequence that need to be checked.

  Params:     delta:            sequence unit, (e.g., the number 2 in a sequence
                                would be marked as 2 * delta)
              loc1:             sequence location of first trigger
              sequence:         sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:           maximum number of times in sequence, after which
                                the sequence restarts
              minTime:          minimum time in data, used to find minimum time
                                to be checked in the sequence
              trigger1Time:     time of first trigger in pair

  Optional:   plot:             boolean, determines if lists necessary for plots
                                need to be filled, False by default

  Returns:    backwardSequence:  list containing time locations to be checked for
              backwardLocs:      list containing the numbers in sequence that
                                 each backward trigger corresponds to
  '''

  # initialise sequence list
  backwardSequence = []

  # initialise sequence location list needed for plotting
  if plot == True:
    backwardLocs = []

  # initialise time location
  time = trigger1Time

  # initialise sequence location
  seqLoc = loc1 - 1

  # define exit condition
  while time > minTime:

    # reset sequence if it reaches the max step
    if seqLoc < 1:
      seqLoc = maxSeq

    # add to time until maximum is reached
    time -= delta * sequence(seqLoc, seqLoc)[0]

    # add time location to list, if more than min time
    if time > minTime:
      backwardSequence.append(time)

      # add to sequence location list
      if plot == True:
        backwardLocs.append(sequence(seqLoc, seqLoc)[0])

    # update sequence location
    seqLoc -= 1

  # reverse lists to respect chronology
  backwardSequence.reverse()
  if plot == True:
    backwardLocs.reverse()

  # include sequence locations in return if needed for plotting
  if plot == True:
    return backwardSequence, backwardLocs

  # return full list
  return backwardSequence, list()


def locBetween(delta, loc1, loc2, sequence, maxSeq, trigger1Time, trigger2Time, plot=False):
  '''
  Finds locations of time points forward in sequence that need to be checked.

  Params:     delta:            sequence unit, (e.g., the number 2 in a sequence
                                would be marked as 2 * delta)
              loc1:             sequence location of first trigger
              loc2:             sequence location of second trigger
              sequence:         sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:           maximum number of times in sequence, after which
                                the sequence restarts
              trigger1Time:     time of first trigger in pair
              trigger2Time:     time of second trigger in pair

  Optional:   plot:             boolean, determines if lists necessary for plots
                                need to be filled, False by default

  Returns:    midSequence:      list containing time locations to be checked for
              midLocs:          list containing the numbers in sequence that
                                each intermediate trigger corresponds to
  '''

  # initialise sequence list
  midSequence = []

  # initialise sequence location list needed for plotting
  if plot == True:
    midLocs = []

  # initialise time location
  time = trigger1Time

  # initialise sequence location
  seqLoc = loc1

  # define exit condition
  while seqLoc < loc2:

    # add to time until maximum is reached
    time += delta * sequence(seqLoc, seqLoc)[0]

    # add time location to list
    midSequence.append(time)

    # add to sequence location list
    if plot == True:
      midLocs.append(sequence(seqLoc, seqLoc)[0])

    # update sequence location
    seqLoc += 1

  # add final sequence location
  if plot == True:
    midLocs.append(sequence(loc2, loc2)[0])

  # include sequence locations in return if needed for plotting
  if plot == True:
    return midSequence, midLocs

  # return full list
  return midSequence, list()