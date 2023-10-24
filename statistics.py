import math


def combinations(numberOfSequencePoints, numberOfTriggers):
  '''
  Finds the number of different possibilities to assign triggers to
  candidate locations for a given sequence.

  Params:   numberOfSequencePoints:  number of candidate time locations in a
                                     given sequence
            numberOfTriggers:        total number of triggers searched over

  Returns:  number of possibilities
  '''

  # define internal variables
  m = numberOfSequencePoints
  n = numberOfTriggers

  # initialise summing variable
  combin = 0

  # calculate number of possibilities
  for i in range(m-1):
    combin += math.comb(m, i) * math.perm(n, m-i)

  # return sum
  return combin