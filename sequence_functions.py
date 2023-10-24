import sympy as sp

def getFibonacci(a, b):
  '''
  Generates ath to bth Fibonacci numbers. For example, getFibonacci(3, 5)
  returns [2, 3, 5].

  Params:     a, b:       range of Fibonacci numbers to generate
  Returns:    FiboList:   list of first n Fibonacci numbers
  '''

  FiboList = []

  # initialise previous terms
  f2 = 0
  f1 = 1

  for i in range(b+1):

    if i <= 1:
      fi = i
    else:
      fi = f2 + f1

      # redefine terms
      f2 = f1
      f1 = fi

    # append to list within desired bounds
    if i >= a:
      FiboList.append(fi)

  return FiboList


def getPrimes(a, b):
  '''
  Generates ath to bth prime numbers. For example, getPrimes(3, 5)
  returns [5, 7, 11].

  Params:     a, b:        range of prime numbers to generate
  Returns:    primesList:  list of first n prime numbers
  '''

  primesList = []

  for i in range(a, b+1):
    primesList.append(sp.prime(i))

  return primesList