
# a dict that gives "not " if the input is False and "" otherwise
vocalize = {False: "not ", True: ""}

# numerical properties:
# r = number of reticulations
# lvl = level
# rsi = maximum indegree of a connected component consisting only of reticulations
# ur = number of unstable roots
# urb = max. number of unstable roots per block
# nd = nesting depth
# sp = max. over all reticulations r of the shortest path connecting the parents of r that avoids r
# srh = max. small reticulation height of any reticulation (shortest distance to a lowest vertex that is ancestor of both parents of r)
# brh = max. big reticulation height of any reticulation (farthest distance to a lowest vertex that is ancestor of both parents of r)
numerical_properties = ['r', 'lvl', 'rsi', 'ur', 'urb', 'nd', 'sp', 'srh', 'brh', 'ua']
network_types = ['tc', 'ntc', 'gs', 'ts', 'rv', 'cv', 'cp', 'ns', 'tb']
# implication of network types (NOTE: only the transitive reduction is given)
# step 1: positive implications: if a network is key, then it is also value
positive_implications = {'tc' : ['ntc', 'ns'],
                         'ntc': ['gs'],
                         'gs' : ['ts', 'rv'],
                         'rv' : ['cv'],
                         'ns' : ['cv'],
                         'gt' : ['rv']}
# step 1: negative implications: if a network is not key, then it is also not value
negative_implications = dict((x,[y for y in positive_implications if x in positive_implications[y]]) for x in network_types)

# given a list X of lists of the same size, return a list L of indices of X such that L[i] is pareto minimal in coordinate i
def pareto_mins(X):
  # return -1 if x is pareto smaller than y, 1 if y is pareto smaller than x, 0 otherwise
  def pareto_cmp(x,y):
    result = cmp(x[0], y[0])
    for i in xrange(1, len(x)):
      comp = cmp(x[i], y[i])
      # if comp == 0, everything stays as it is
      if comp != 0:
        # if comp is not 0 but the inverse of result, noone is pareto-anything
        if comp + result == 0:
          return 0
    return result
  
  result = []
  for L in X:
    do_add = True
    for i in xrange(len(result)-1,-1,-1):
      comp = pareto_cmp(result[i], L)
      
      if comp == 1:
        # if L is pareto smaller than result[i], remove result[i] (replace with last result)
        tmp = result.pop()
        if i < len(result):
          result[i] = tmp
      elif comp == -1:
        # if someone in result is pareto smaller than L, then forget about L and get the next list
        do_add = False
        i = 0
    # if we didn't find anyone that is pareto better than L, then we add L to the results
    if do_add:
      result.append(L)
      
  print("pareto mins of %s:\n %s" % (str(X), str(result)))
  return result
        

# error class for encountering vertices not matching indeg == 1 XOR outdeg == 1
class DegreeError(ValueError):
    def __init__(self, arg):
        self.strerror = arg
        self.args = {arg}

# error class for exiting multiple encapsulated loops
class Break(ValueError):
    pass

class NetworkProperty:
  # long description of the property
  short = None

  # short abbreviation of the property
  long = None

  # pointer to the host network the property is of
  network = None

  # the actual value of the property (or None if not yet computed)
  val = None

  # a logging function
  log = None

  # NOTE: the actual computation of the value is done in self.check() which is not implemented in this class.
  #       Each child class is supposed to implement their own check()
  def __init__(self, _long, _short, _network):
    self.long = _long
    self.short = _short
    self.network = _network
    self.log = _network.log
    self.init()

  def init(self):
      pass

  # set the value
  def set(self, _val = True):
    if self.val != _val:
      self.log('setting ' + self.long + ' to ' + str(_val))
      self.val = _val
      if _val:
        for impl in positive_implications.get(self.short, []):
          self.network.properties[impl].set()
      else:
        for impl in negative_implications.get(self.short, []):
          self.network.properties[impl].unset()

  # unset the value
  def unset(self):
    self.set(False)


  # sets the value only if it has not been previously set
  def set_if_none(self, _val = True):
    if self.val is None:
      self.val = _val


  # get the value of the property, computing it if needed (using self.check())
  def get(self):
    if self.val is None:
      self.check()
    return self.val


  # conversion to bool
  def __nonzero__(self):
    return self.get()

  # Python3 compat
  __bool__ = __nonzero__

  # report about a property
  def report(self):
    self.log("=== N is " + vocalize[self.get()] + self.long + " ===")
    return vocalize[self.get()] + self.short



class NumericalNetworkProperty(NetworkProperty):
  # report about a value
  def report(self):
    self.log("=== " + self.long + ": " + str(self.get()) + " ===")
    return str(self.get())   


