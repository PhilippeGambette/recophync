

import re
import random
import networkx
from collections import defaultdict, deque

from properties.common import numerical_properties, network_types
from properties.bool_prop import *
from properties.num_prop import *



class PhyloNetwork:
  # whether or not verbose output is required
  verbose = None

  # the actual network
  N = None

  # the root of the network
  root = None

  # set of reticulations
  reticulations = None

  # set of leaves
  leaves = None

  # immediate dominators of each vertex in N
  dominators = None

  # a dict mapping each vertex to ONE OF THE leaves it is stable on
  stability = None

  # network regularity (2 = binary, 0 = irregular)
  regular = None

  # the property objects that will compute the properties
  properties = None

  def __init__(self, _verbose = False):
    self.verbose = _verbose
    self.N = networkx.DiGraph()
    self.root = None
    self.reticulations = set()
    self.leaves = set()
    self.stability = dict()
    self.lowest_stable = dict()
    self.regular = None
     # map short descriptions to functions checking classes
    self.properties = {'tc' : TreeChildProp('tree child', 'tc', self),
                       'ntc': NearlyTreeChildProp('nearly tree child','ntc', self),
                       'gs' : GeneticallyStableProp('genetically stable','gs', self),
                       'ts' : TreeSiblingProp('tree sibling','ts', self),
                       'rv' : ReticulationVisibleProp('reticulation visible','rv', self),
                       'cv' : ComponentVisibleProp('component visible','cv', self),
                       'cp' : CompressedProp('compressed', 'cp', self),
                       'ns' : NearlyStableProp('nearly stable', 'ns', self),
                       'gt' : GalledTreeProp('galled tree', 'gt', self),
                       'tb' : TreeBasedProp('tree based','tb', self),
                       # numerical properties
                       'r'  : NumReticulations('#reticulations', 'r', self),
                       'lvl': Level('level', 'lvl', self),
                       'nd' : NestingDepth('nesting depth', 'nd', self),
                       'ur' : NumUnstableRoots('#unstable component roots', 'ur', self),
                       'urb': NumUnstableRootsPerBlock('#unstable components roots per block', 'urb', self),
                       'rsi': MaxReticulationSubgraphIndegree('max #incoming edges to any reticulation component', 'rsi', self),
                       'sp' : ShortestPath('max shortest path above reticulation', 'sp', self),
                       'srh': SmallReticulationHeight('min distance of any lowest vertex seeing both parents', 'srh', self),
                       'brh': BigReticulationHeight('max distance of any lowest vertex seeing both parents', 'brh', self),
                       'ua' : UncommonAncestors('max number of uncommon ancestors of parents of any reticulation', 'ua', self)
                      }

  # Input: log string
  # Effect: print log string if verbose flag is raised
  def log(self, s):
    if self.verbose:
      print(s)


  # Input: rooted network N, some vertex u
  # Output: the type of the node as string: 'root', 'leaf', 'tree', 'reticulation'
  def node_type(self, u):
    ideg = self.N.in_degree(u)
    if ideg == 1:
      return 'leaf' if self.N.out_degree(u) == 0 else 'tree'
    elif ideg > 1:
      return 'reticulation'
    else:
      return 'root'


  def is_reticulation(self, v):
    return self.N.in_degree(v) > 1


  def is_stable(self, v):
    if self.node_type(v) == 'leaf':
      return True
    else:
      return self.stability.get(v, False)


  # Input: some edge uv
  # Effect: contract uv in N unto v
  def contract_edge(self, uv, delete_node = True):
    u,v = uv
    self.N.remove_edge(u,v)
    self.N.add_edges_from((v,w) for w in self.N.successors_iter(u))
    self.N.add_edges_from((w,v) for w in self.N.predecessors_iter(u))
    if delete_node:
      self.N.remove_node(u)
    else:
      self.N.remove_edges_from((u,w) for w in self.N.successors_iter(u))
      self.N.remove_edges_from((w,u) for w in self.N.predecessors_iter(u))

  # Input: a vertex u
  # Effect: replace u with the join of predecessors(u) with successors(u)
  # NOTE: this will turn length-2 cycles into loops
  def shortcut_over(u):
    for v in self.N.predecessors_iter(u):
      for w in self.N.successors_iter(u):
        self.N.add_edge(v,w)
    self.N.remove_node(u)

  # Input: edge (u,v)
  # Output: subdivide uv and return the new vertex w, as well as a list of new edges
  def subdivide(self, uv):
    w = self.N.number_of_nodes()
    new_edges = [(uv[0],w),(w,uv[1])]
    self.N.remove_edge(uv[0], uv[1])
    self.N.add_edges_from(new_edges)
    return w, new_edges


  # Input: a rooted phylogenetic network N
  # Output: N where each indegree and outdegre 1 vertex is contracted with its parent
  def _contract(self):
    for v in self.N.nodes():
      #self.log('considering ' + str(v) + ' with indeg ' + str(self.N.in_degree(v)) + ' & outdeg ' + str(self.N.out_degree(v)))
      if self.N.in_degree(v) == self.N.out_degree(v) == 1:
        u = next(self.N.predecessors_iter(v))
        w = next(self.N.successors_iter(v))
        self.log(str(v) + " was indegree-1 outdegree-1, so I removed it.")
        self.N.add_edge(u, w)
        self.N.remove_node(v)


  # input: none
  # task: fill stability relation for all non-leaves
  def _update_stable(self):
    # compute stable vertices
    self.log("""== Computing the stable vertices ==""")

    # compute immediate dominators
    self.dominators = networkx.immediate_dominators(self.N, self.root)

    # for each vertex, assign a leaf that they are stable on
    self.stability = dict((u,u) for u in self.leaves)
    X = list(self.leaves.copy())
    while X:
      u = X.pop()
      v = self.dominators[u]
      if v not in self.stability:
        self.stability[v] = u
        X.append(v)

    if self.verbose:
      self.log('immediate dominators:')
      for u in self.stability:
        if self.is_stable(u):
          self.log(str(u) + ': ' + str(self.stability[u]))


      

  # Input: none
  # Effect: set root to the first vertex of in-degree 0 and set regularity and reticulations and stable vertices
  def _update_infos(self):
    self.root = None
    self.regular = None
    self.log('checking the network for strangeness...')
    for i in self.N:
      indeg = self.N.in_degree(i)
      outdeg = self.N.out_degree(i)
      if indeg == 1 and outdeg == 1:
        raise DegreeError(str(i) + " has indegree & outdegree 1")
      if indeg > 1 and outdeg > 1:
        raise DegreeError(str(i) + " has indegree & outdegree > 1")
      if indeg == 0:
        if self.root is None:
          self.root = i
        else:
          raise DegreeError("more than one root in the network: " + str(i) + " & " + str(self.root))
      else:
        if outdeg > 0:
          if self.regular is None:
            self.regular = indeg * outdeg
          elif self.regular != indeg * outdeg:
            self.regular = 0
      if indeg > 1:
        self.reticulations.add(i)
      if outdeg == 0:
        self.leaves.add(i)
    if self.root is None:
      raise DegreeError("no root in network")

    if self.N.out_degree(self.root) != self.regular:
      self.regular = 0

    self.log('N has ' + str(self.N.number_of_nodes()) + " vertices and " \
        + str(self.N.number_of_edges()) + " edges." \
        + " Its root is " + str(self.root) + "."
        + " N is " + ("not " if self.regular == 0 else str(self.regular) + "-") + "regular")
    self.log('reticulations: ' + str(self.reticulations))
    self.log('leaves: ' + str(self.leaves))


  # Input: text file containing a list of arcs
  # Effect: read network from the file
  def open(self, filename):
    with open(filename) as fd:
      for line in fd.readlines():  # match patterns of type "vertex1 vertex2"
        res = re.search("^([^ ]*) ([^ \n\r]*)$", line)
        if res:
          self.N.add_edge(res.group(1), res.group(2))
    self._contract()
    self._update_infos()
    self._update_stable()


  # Input: n, regularity
  # Effect: create regular spanning tree on n vertices if possible
  def create_regular_tree(self, n, regularity = 2):
    if (n - 1) % regularity != 0:
      raise DegreeError('no ' + str(regularity) + '-regular tree with ' + str(n) + ' nodes exist')

    self.log('creating ' + str(regularity) + '-regular tree on ' + str(n) + ' nodes')
    rand = random.Random()
    current_leaves = [0]
    for i in xrange((n - 1)//regularity):
      leaf = current_leaves.pop(int(rand.uniform(0, len(current_leaves))))
      for j in xrange(regularity):
        index = i * regularity + j + 1
        self.N.add_edge(leaf, index)
        current_leaves.append(index)
    return current_leaves


  # Input: n
  # Effect: create a random BINARY network with n nodes & reticulation number distributed exponentially with mean r_mean
  def create_binary(self, num_nodes, r_mean):
    # step 0: get eligible numbers of reticulations for n & regularity
    n = num_nodes
    if n % 2 == 0:
      raise DegreeError('no network with even number of nodes (' + str(n) + ') can be binary')

    rand = random.Random()
    retis = min(int(rand.expovariate(1.0 / r_mean)), (n-1)//2 - 1);

    # NOTE: we use m = sum of in-degrees = sum of out-degrees, assuming the root is counted in t:
    # (1) regularity * r + (n - r) = m + 1    NOTE: +1 needed for the root
    # (2) regularity * t + r = m
    
    # step 2: compute r, t, and l
    m = retis + n - 1
    # retis = (m - n + 1) / (regularity - 1)
    treenodes = (m - retis) // 2
    leaves = n - retis - treenodes
    self.log('creating random network with ' + str(treenodes) + ' tree nodes, ' + str(retis) + ' reticulations, and ' + str(leaves) + ' leaves ( = ' + str(num_nodes) + ' nodes)')

    # step 3: create binary spanning tree on treenodes - retis + leaves
    self.create_regular_tree(n - 2*retis)

    # step 4: repeatedly add an edge between 2 random edges, creating a new tree node and reticulation each time
    # NOTE: ensure time consistency by giving each node a time and each new node the mean of its old parent and child
    self.log('adding reticulations...')
    edges = self.N.edges()
    time = dict((x,x) for x in xrange(n - 2*retis))
    for i in xrange(retis):
      e_idx = rand.sample(xrange(len(edges)), 2)

      w = [0, 1]
      for i in (0,1):
        e = edges[e_idx[i]]
        w[i], new_edges = self.subdivide(e)
        time[w[i]] = (time[e[0]] + time[e[1]]) // 2
        # update edges
        edges[e_idx[i]] = new_edges.pop()
        edges.extend(new_edges)
      
      self.N.add_edge(w[0], w[1])
      edges.append((w[0], w[1]))

    self._update_infos()
    self._update_stable()


  # Input: n
  # Effect: create a random network with n nodes &
  #         reticulation number distributed exponentially with mean r_mean
  #         the network results from contracting a number of edges that is normally distributed around contraction_mean with sigma = n/10
  def create_random(self, num_nodes, r_mean, contractions_mean):
    self.log('creating random network with ' + str(num_nodes) + ' nodes')
    rand = random.Random()
    n = num_nodes
    contractions = int(rand.gaussvariate(non_binary_mean, n/10.0));
    if n + contractions % 2 == 0:
      contractions += 1

    # step 0: create random binary network with n + contractions nodes
    self.create_binary(n + contractions, r_mean)

    # step 1: contract 'contractions' edges in N
    # NOTE: we need to keep an eye on vertices with both in- and out-degree > 1. They will be uncontracted later
    invalid_nodes = set()
    edges = set(self.N.edges())
    i = 0
    while i < contractions + len(invaild_nodes):
      e = rand.randint(0, len(edges) - 1)
      if e[0] in invalid_nodes:
        invalid_nodes.discard(e[0])
        invalid_nodes.add(e[1])
      elif e[1] not in invalid_nodes:
        if self.is_reticulation(e[0]) != self.is_reticulation(e[1]):
          invalid_nodes.add(e[1])
      self.contract_edge(e)

    # step 2: uncontract invalid nodes
    for u in invalid_nodes:
      v = self.N.add_node()
      for w in self.N.successors_iter(u):
        self.N.remove_edge(u, w)
        self.N.add_edge(v, w)
      self.N.add_edge(u, v)

    self._update_infos()
    self._update_stable()
     


  # ===================== computing values =========================

  # Input: a type t (0,1,2,3)
  # Output: list of vertices of N of type t
  def get_nodes(self, t):
    return (u for u in self.N.nodes() if self.node_type(u) == t)

  # Input: a type t (0,1,2,3)
  # Output: list of vertices of N of type t
  def count_nodes(self, t):
    if t == 'root':
      return 1
    else:
      if self.regular:
        # we consider the root to be a tree node
        n = self.N.number_of_nodes()
        m = self.N.number_of_edges()
        retis = (m - n + 1) // (self.regular - 1)
        tree = (m - retis) // self.regular
        leaf = n - retis - tree
        return {'leaf': leaf,
                'tree': tree,
                'reticulation': retis}.get(t, 0)
      else:
        return len(self.get_nodes(t))

  # input: a list of nodes
  # output: the size of the largest sublist that intersects a block in N
  def max_per_block(self, v_set):
    biconn_comp = networkx.biconnected_components(self.N.to_undirected())
    lists_per_block = [v_set.intersection(B) for B in biconn_comp if len(B) > 2]
    self.log('per block: ' + str(lists_per_block))
    return max((len(x) for x in lists_per_block)) if lists_per_block else 0

  
  def write_dot(self, filename):
    networkx.write_dot(self.N, filename)





