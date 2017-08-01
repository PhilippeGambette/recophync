#!/usr/bin/python

import argparse
import os
import re
import networkx
import Queue
import sys
import random

#############################################################################
# RecoPhyNC - recognizing phylogenetic networks classes
# 
# This program tests whether a phylogenetic network, given as a file
# containing a list of arcs, is level-k (computeLevel), tree-child 
# (isTreeChild), nearly tree-child (isNearlyTreeChild), genetically
# stable (isGeneticallyStable), reticulation-visible (isReticulationVisible),
# tree-sibling (isTreeSibling), compressed (isCompressed), nearly stable
# (isNearlyStable).
# 
# To use it you need to install the NetworkX library:
# https://networkx.github.io/
# 
# You need to put all your network files in the data folder contained
# in the same folder as this file.
# 
#############################################################################
# Copyright (C) 2015-2016 Philippe Gambette, Anthony Labarre
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

#verbose = True
VERBOSE = False

print """
RecoPhyNC Copyright (C) 2015 Philippe Gambette
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions (GNU General Public License).
"""

# a dict that gives "not " if the input is False and "" otherwise
vocalize = {False: "not ", True: ""}

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
negative_implications = dict([(x,[y for y in positive_implications if x in positive_implications[y]]) for x in network_types])

# error class for encountering vertices not matching indeg == 1 XOR outdeg == 1
class DegreeError(ValueError):
  def __init__(self, arg):
    self.strerror = arg
    self.args = {arg}


class PhyloNetwork:
  # whether or not verbose output is required
  verbose = None

  # the actual network
  N = None

  # the root of the network
  root = None

  #list of reticulations
  reticulations = set()

  # list of leaves
  leaves = []

  # an array mapping each vertex to the set of leaves it is stable on
  stability = dict()

  # network regularity (2 = binary, 0 = irregular)
  regular = None

  # store identified classes and useful properties in a dictionary,
  # so we can use previously acquired knowledge to speed things up
  # (i.e., if we know that N is tree-child, then we know it's also
  # nearly tree-child so no need to ask a function for the result)
  properties = dict()


  def __init__(self, _verbose = False):
    self.verbose = _verbose
    self.N = networkx.DiGraph()
    self.properties = dict()
  
  # Input: log string
  # Effect: print log string if verbose flag is raised
  def log(self, s):
    if self.verbose:
      print s

  # Input: rooted network N, some vertex u
  # Output: 0 if u is a reticulation, 1 if u is a tree vertex, 2 if u is a leaf, 3 if u is the root
  def node_type(self, u):
    return {
        0: 'root',
        1: 'leaf' if self.N.out_degree(u) == 0 else 'tree'
        }.get(self.N.in_degree(u), 'reticulation')

  def is_reticulation(self, v):
    return self.N.out_degree(v) == 1

  def is_stable(self, v):
    if self.node_type(v) == 'leaf':
      return True
    else:
      return self.stability.get(v, False)

  # Input: a rooted phylogenetic network N
  # Output: N where each indegree and outdegre 1 vertex is contracted with its parent
  def _contract(self):
    """
    :return:
    """

    X = set(self.N.nodes()) #NOTE: we need a set here to avoid having a node multiple times in the queue
    while X:
      v = X.pop()
      #self.log('considering ' + str(v) + ' with indeg ' + str(self.N.in_degree(v)) + ' & outdeg ' + str(self.N.out_degree(v)))
      if self.N.in_degree(v) == self.N.out_degree(v) == 1:
        u = self.N.predecessors(v)[0]
        w = self.N.successors(v)[0]
        self.log(str(v) + " was indegree-1 outdegree-1, so I removed it.")
        self.N.add_edge(u, w)
        self.N.remove_node(v)
        X.add(u)
        X.add(w)

  # input: none
  # task: fill stability relation for all non-leaves
  def _update_stable(self):
    # compute stable vertices
    self.log("""== Computing the stable vertices ==""")
    #t0=datetime.datetime.now()

    # the current bottom-up front
    X = Queue.Queue()
    map(X.put, self.leaves)

    # the current number of vertices that can see x
    spread = dict([(x,1) for x in self.leaves])
    # for each vertex, the set of leaves it can reach
    reach = dict([(x,set()) for x in self.N.nodes()])
    reach.update(dict([(x,set([x])) for x in self.leaves]))
    # for each vertex, how many successors we have processed
    childs_seen = dict([(x,None) for x in self.leaves])

    while not X.empty():
      i = X.get(False)

      for j in reach[i]:
        spread[j] -= 1
      for p in self.N.predecessors(i):
        childs_seen[p] = childs_seen.get(p, 0) + 1
        for j in reach[i]:
          if j not in reach[p]:
            spread[j] = spread.get(j, 0) + 1
          reach[p].add(j)
        if childs_seen[p] == self.N.out_degree(p):
          X.put(p)
          self.stability[p] = set([x for x in reach[p] if spread[x] == 1])
      del reach[i]
      del childs_seen[i]
    #print "Time Stable vertices: "+str((datetime.datetime.now()-t0).microseconds)+"ms."
    
    if self.verbose:
      self.log('Stable vertices: ' + str([x for x in self.N.nodes() if self.is_stable(x) ]))

      

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
          raise DegreeError("more than one root in the network: " + str(i) + " & " + str(root))
      else:
        if outdeg > 0:
          if self.regular is None:
            self.regular = indeg * outdeg
          elif self.regular != indeg * outdeg:
            self.regular = 0
      if indeg > 1:
        self.reticulations.add(i)
      if outdeg == 0:
        self.leaves.append(i)
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

  # Input: edge (u,v)
  # Output: subdivide uv and return the new vertex w, as well as a list of new edges
  def subdivide(self, uv):
    w = self.N.number_of_nodes()
    new_edges = [(uv[0],w),(w,uv[1])]
    self.N.remove_edge(uv[0], uv[1])
    self.N.add_edges_from(new_edges)
    return w, new_edges

  # Input: n, regularity
  # Effect: create regular spanning tree on n vertices if possible
  def create_regular_tree(self, n, regularity = 2):
    if (n - 1) % regularity != 0:
      raise DegreeError('no ' + str(regularity) + '-regular tree with ' + str(n) + ' nodes exist')

    self.log('creating ' + str(regularity) + '-regular tree on ' + str(n) + ' nodes')
    rand = random.Random()
    current_leaves = [0]
    for i in xrange((n - 1)/regularity):
      leaf = current_leaves.pop(int(rand.uniform(0, len(current_leaves))))
      for j in xrange(regularity):
        index = i * regularity + j + 1
        self.N.add_edge(leaf, index)
        current_leaves.append(index)
    return current_leaves

  # Input: list of node types
  # Output: a random node of this type if include = True, otherwise a random node not of this type
  def random_node(self, allowed_types, include = True):
    rand = random.Random()
    u = rand.randint(0, self.N.number_of_nodes() - 1)
    while (self.node_type(u) in allowed_types) != include:
      u = rand.randint(0, self.N.number_of_nodes() - 1)
    return u

  # Input: m
  # Output: m edges picked uniformly at random
  def random_edges(self, num_edges):
    m = self.N.number_of_edges()
    if num_edges > m:
      raise ValueError('not enough edges in the network')
    
    rand = random.Random()
    return rand.sample(self.N.edges(), num_edges)
    

  # Input: n
  # Effect: create a random BINARY network with n nodes & reticulation number distributed exponentially with mean r_mean
  def create_binary(self, num_nodes, r_mean):
    # step 0: get eligible numbers of reticulations for n & regularity
    n = num_nodes
    if n % 2 == 0:
      raise DegreeError('no network with even number of nodes (' + str(n) + ') can be binary')

    rand = random.Random()
    retis = min(int(rand.expovariate(1.0 / r_mean)), (n-1)/2 - 1);

    # NOTE: we use m = sum of in-degrees = sum of out-degrees, assuming the root is counted in t:
    # (1) regularity * r + (n - r) = m + 1    NOTE: +1 needed for the root
    # (2) regularity * t + r = m
    
    # step 2: compute r, t, and l
    m = retis + n - 1
    # retis = (m - n + 1) / (regularity - 1)
    treenodes = (m - retis) / 2
    leaves = n - retis - treenodes
    self.log('creating random network with ' + str(treenodes) + ' tree nodes, ' + str(retis) + ' reticulations, and ' + str(leaves) + ' leaves ( = ' + str(num_nodes) + ' nodes)')

    # step 3: create binary spanning tree on treenodes - retis + leaves
    self.create_regular_tree(n - 2*retis)

    # step 4: repeatedly add an edge between 2 random edges, creating a new tree node and reticulation each time
    # NOTE: ensure time consistency by giving each node a time and each new node the mean of its old parent and child
    self.log('adding reticulations...')
    edges = self.N.edges()
    time = dict([(x,x) for x in xrange(n - 2*retis)])
    for i in xrange(retis):
      e_idx = rand.sample(xrange(len(edges)), 2)

      w = [0, 1]
      for i in (0,1):
        e = edges[e_idx[i]]
        w[i], new_edges = self.subdivide(e)
        time[w[i]] = (time[e[0]] + time[e[1]]) / 2
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

    # step 0: create random binary network with n + non_binary nodes
    self.create_binary(n + contractions, r_mean)

    # step 1: contract 'contractions' edges in N
    # keep an eye on vertices with both in- and out-degree > 1. They will be uncontracted later
    invalid_nodes = set()
    for i in xrange(contractions):
      u = rand.randint(0, self.N.number_of_nodes())
      
      

    self._update_infos()
    self._update_stable()
     


  # Input: some edge uv
  # Effect: contract uv in N
  def contractEdge(self, uv, deleteNode = True):
    u = uv[0]
    v = uv[1]
    self.N.remove_edge(u,v)
    self.N.add_edges_from([(v,w) for w in self.N.successors(u)])
    self.N.add_edges_from([(w,v) for w in self.N.predecessors(u)])
    if deleteNode:
      self.N.remove_node(u)
    else:
      self.N.remove_edges_from([(u,w) for w in self.N.successors(u)])
      self.N.remove_edges_from([(w,u) for w in self.N.predecessors(u)])

  # Input: a vertex u
  # Effect: replace u with the join of predecessors(u) with successors(u)
  # NOTE: this will turn length-2 cycles into loops
  def shortcut_over(u):
    for v in self.N.predecessors(u):
      for w in self.N.successors(u):
        self.N.add_edge(v,w)
    self.N.remove_node(u)


  # ===================== computing values =========================

  # Input: a type t (0,1,2,3)
  # Output: list of vertices of N of type t
  def get_nodes(self, t):
    return [u for u in self.N.nodes() if self.node_type(u) == t]

  # Input: a type t (0,1,2,3)
  # Output: list of vertices of N of type t
  def count_nodes(self, t):
    if self.regular:
      if t == 'root':
        return 1
      else:
        # we consider the root to be a tree node
        n = self.N.number_of_nodes()
        m = self.N.number_of_edges()
        retis = (m - n + 1) / (self.regular - 1)
        tree = (m - retis) / self.regular
        leaf = n - retis - tree
        return {'leaf': leaf,
                'tree': tree,
                'reticulation': retis}.get(t, 0)
    else:
      return len(self.get_nodes(t))

  # Output: #reticulations
  def count_reticulations(self):
    return self.count_nodes('reticulation')

  # output: the list of unstable component roots in N
  def list_unstable_components(self):
    for u in self.reticulations:
      if self.node_type(self.N.successors(u)[0]) == 'tree':
        self.log('component-root ' + str(u) + ' is stable on ' + str(self.stability.get(u, [])))
    unstable_roots = set([u for u in self.reticulations if not self.is_reticulation(self.N.successors(u)[0]) and not self.is_stable(u)])
    # if there are unstable roots, we cannot have component visibility
    self.set_property('cv', not unstable_roots)
    self.log("unstable roots: " + str(unstable_roots))
    return unstable_roots

  # output: the number of unstable component roots in N
  def count_unstable_components(self):
    return len(self.list_unstable_components())

  # input: a list of nodes
  # output: the size of the largest sublist that intersects a block in N
  def max_per_block(self, v_set):
    biconn_comp = networkx.biconnected_components(self.N.to_undirected())
    lists_per_block = [ v_set.intersection(B) for B in biconn_comp if len(B) > 2]
    self.log('per block: ' + str(lists_per_block))
    return max([len(x) for x in lists_per_block]) if lists_per_block else 0

  # Input: binary rooted network N
  # Output: level of N
  def compute_level(self):
    return self.max_per_block(self.reticulations)

  # output: the number of unstable component roots in N
  def count_ucr_per_block(self):
    return self.max_per_block(self.list_unstable_components())


  # output: the highest number of vertices per block that are unstable component roots in N
  def compute_ucr_per_block(self):
    return self.max_per_block(set(self.list_unstable_components()))


  # Input: binary rooted network N
  # Output: -1 if not nested, nested depth otherwise
  # Does not work yet
  def compute_nested_depth():
      """Returns the nested depth of a rooted binary network N, or -1 if N is not
       nested.

      :rtype : int
      """
      biconn_comp = networkx.biconnected_components(self.N.to_undirected())
      for B in biconn_comp:
          # compute the nested depth of each biconnected component
          for node in B:
              print node
              #...
      return -1

  # output: the highest number of incoming arcs into any connected component in N[R]
  def max_reti_subgraph_indeg(self):
    # contract all edges between reticulations
    parents_seen = dict()
    collective_indeg = dict()
    X = Queue.Queue()
    X.put(self.root)

    while not X.empty():
      v = X.get_nowait()
      for s in self.N.successors(v):
        if self.is_reticulation(s):
          collective_indeg[s] = collective_indeg.get(s,0) \
                              + (collective_indeg[v] if self.is_reticulation(v) else 1)
        parents_seen[s] = parents_seen.get(s, 0) + 1
        if parents_seen[s] == self.N.in_degree(s):
          X.put(s)
    return max(collective_indeg.values() or [0])


  # map short descriptions to functions computing values
  short_to_val = {'ret': ('#reticulations', count_reticulations),
                  'lvl': ('level', compute_level),
                  'rcl': ('max incoming arcs into any reticulation component', max_reti_subgraph_indeg),
                  'uc' : ('#unstable components', count_unstable_components),
                  'ucb': ('#unstable component roots / block', count_ucr_per_block)
                 }

  # report about a value
  def report_value(self, short):
    description = self.short_to_val[short][0]
    val = self.short_to_val[short][1](self)
    self.log("=== " + description + ": " + str(val) + " ===")
    return ";" + str(val)   


  # ===================== checking properties =========================

  def set_property(self, prop, val = True):
    self.properties[prop] = val
    if val:
      for impl in positive_implications.get(prop, []):
        if impl not in self.properties:
          self.log('deriving ' + impl + ' from ' + prop)
          self.set_property(impl, True)
    else:
      for impl in negative_implications.get(prop, []):
        if impl not in self.properties:
          self.log('deriving not-' + impl + ' from not-' + prop)
          self.set_property(impl, False)

  def unset_property(self, prop):
    self.set_property(prop, False)

  # Input: rooted network N
  # Output: True if N is tree-child, False otherwise
  def is_tree_child(self):
    """
    :return:
    """
    if 'tc' in self.properties:
      return self.properties['tc']

    successors = self.N.successors_iter
    for nonleaf in (v for v, degree in self.N.out_degree_iter() if degree):
        treechildNode = any(self.N.in_degree(c) == 1 for c in successors(nonleaf))
        if not treechildNode:
            self.log(str(nonleaf) + " is not treechild.")
            self.unset_property('tc')
            if not self.verbose:
              return False
    if 'tc' not in self.properties:
      self.set_property('tc')
      return True
    else:
      return False

  # Output: True if N is nested, False otherwise
  def is_nested(self):
    pass
    #TODO: write me

  # Output: True if N is a galled tree, False otherwise
  def is_galled_tree(self):
    if 'gt' in self.properties:
      return self.properties['gt']
    else:
      result = (self.compute_level() == 1)
      self.set_property('gt', result)
      return result

  # Output: True if N is tree-child, False otherwise
  def is_tree_sibling(self):
    """
    :return:
    """

    if 'ts' in self.properties:
      return self.properties['ts']

    predecessors = self.N.predecessors_iter
    successors = self.N.successors_iter
    for v in self.reticulations:
        # for each reticulation vertex, we build its set of siblings S
        tree_siblings = [c for p in predecessors(v) for c in successors(p) if not self.is_reticulation(c)]
        # we check if S contains at least one tree vertex
        if tree_siblings:
          self.log(str(v) + " has tree siblings: " + str(tree_siblings))
        else:
          self.log(str(v) + " is not tree-sibling.")
          self.unset_property('ts')
          if not self.verbose:
            return False
    if not 'ts' in self.properties:
      self.set_property('ts')
      return True
    else:
      return False

  # Input: a rooted phylogenetic network N
  # Output: True if N is reticulation visible, False otherwise
  def is_reticulation_visible(self):
    """
    :return:
    """

    if 'rv' in self.properties:
      return self.properties['rv']

    for v in self.reticulations:
      if not self.is_stable(v):
        self.log(str(v) + " is not stable.")
        self.unset_property('rv')
        if not self.verbose:
          return False
      else:
        self.log(str(v) + " is stable on " + str(self.stability.get(v, set())))
    if not 'rv' in self.properties:
      self.set_property('rv')
      return True
    else:
      return False



  # Input: a rooted phylogenetic network N
  # Output: True if N is reticulation visible, False otherwise
  def is_nearly_stable(self):
    if 'ns' in self.properties:
      return self.properties['ns']

    # skipping dots may improve speed
    predecessors = self.N.predecessors_iter
    for v in self.N:
      if v != self.root:
        if not self.is_stable(v):
          self.log(str(v) + " is not stable.")

          unstable_parents = [p for p in predecessors(v) if not self.is_stable(p)]

          if unstable_parents:
            self.log("Unstable parents of " + str(v) + ": " + str(unstable_parents))
            self.unset_property('ns')
            if not self.verbose:
              return False
          else:
            self.log("All parents of " + str(v) + " are stable.")
        else:
          self.log(str(v) + " is stable.")
    if not 'ns' in self.properties:
      self.set_property('ns')
      return True
    else:
      return False

  # Input: a rooted phylogenetic network N
  # Output: True if N is component visible, False otherwise
  def is_component_visible(self):
    if 'cv' in self.properties:
      return self.properties['cv']

    for v in self.reticulations:
      s = self.N.successors(v)[0]
      if not self.is_reticulation(s) and not self.is_stable(s):
        self.log(str(s) + " is an unstable component root.")
        self.unset_property('cv')
        if not self.verbose:
          return False
      else:
        self.log(str(s) + " is stable on " + str(self.stability.get(s, [])))
    if not 'cv' in self.properties:
      self.set_property('cv')
      return True
    else:
      return False


  # Input: a rooted phylogenetic network N
  # Output: True if N is compressed, false otherwise
  def is_compressed(self):
    """Returns True if network N is compressed, False otherwise."""

    if 'cp' in self.properties:
      return self.properties['cp']

    # skipping dots may improve speed
    predecessors = self.N.predecessors_iter
    for r in self.reticulations:
      for p in predecessors(r):
        if self.is_reticulation(p):
          self.log("Not compressed: parent " + str(p) + \
                   " of reticulation " + str(r) + \
                   " is also a reticulation.")
          self.unset_property('cp')
          if not self.verbose:
            return False
    if not 'cp' in self.properties:
      self.set_property('cp')
      return True
    else:
      return False

  # input: a vertex and a dict to True/False
  # output: True <=> v has a tree path in N (saving everyone between in the dict)
  def check_TC(self, v, has_TC):
    if not self.is_reticulation(v):
      if v in has_TC:
        return has_TC[v]
      for i in self.N.successors(v):
        has_TC[v] = self.check_TC(i, has_TC)
        if has_TC[v]:
          return True
    return False

  # Input: a rooted phylogenetic network N
  # Output: True if N is nearly tree-child, false otherwise
  def is_nearly_treechild(self):
    if 'ntc' in self.properties:
      return self.properties['ntc']

    has_TC = dict([(self.N.predecessors(v)[0],True) for v in self.leaves])
    
    for r in self.reticulations:
      for p in self.N.predecessors(r):
        if not self.is_reticulation(p):
          if self.check_TC(p, has_TC):
            self.log(str(p) + " has a tree path.")
          else:
            self.log(str(p) + " has no tree path.")
            self.unset_property('ntc')
            if not self.verbose:
              return False
    if not 'ntc' in self.properties:
      self.set_property('ntc')
      return True
    else:
      return False


  # Input: none
  # Output: True if N is genetically stable, False otherwise
  def is_genetically_stable(self):
    if 'gs' in self.properties:
      return self.properties['gs']

    predecessors = self.N.predecessors_iter
    for r in self.reticulations:
      self.log("Testing if " + str(r) + " is stable.")
      if self.is_stable(r):
        self.log(str(r) + " is stable, testing if its parents are stable.")
        stable_parents = [x for x in predecessors(r) if self.is_stable(x)]
        self.log("Stable parents of " + str(r) + ": " + str(stable_parents))
        if len(stable_parents) == 0:
          self.unset_property('gs')
          if not self.verbose:
            return False
      else:
        self.unset_property('gs')
        if not self.verbose:
          return False
    if not 'gs' in self.properties:
      self.set_property('gs')
      return True
    else:
      return False


  # Input: none
  # Output: True if N is tree-based
  # NOTE: we use the characterization of Francis, Semple, Steel, Arxiv'16
  from networkx.algorithms import bipartite
  def is_tree_based(self):
    if 'tb' in self.properties:
      return self.properties['tb']

    # step 1: build G_N
    B = networkx.Graph()
    B.add_nodes_from([str(x)+'a' for x in self.N.nodes()], bipartite=0)
    B.add_nodes_from([str(x)+'b' for x in self.N.nodes()], bipartite=1)
    B.add_edges_from([(str(x)+'a', str(y)+'b') for (x,y) in self.N.edges()])

    # step 2: check if G_N has a matching of size |V| - |X|
    match = networkx.bipartite.maximum_matching(B)
    self.log("found matching of size " + str(len(match)/2) + ": " + str(match) + " (" + str(self.N.number_of_nodes() - len(self.leaves)) + " necessary for tree-based)")

    result = (len(match) / 2 == self.N.number_of_nodes() - len(self.leaves))
    self.set_property('tb', result)
    return result


  # map short descriptions to functions checking classes
  short_to_fn = {'tc': ('tree child', is_tree_child),
                 'ntc': ('nearly tree child', is_nearly_treechild),
                 'gs': ('genetically stable', is_genetically_stable),
                 'ts': ('tree sibling', is_tree_sibling),
                 'rv': ('reticulation visible', is_reticulation_visible),
                 'cv': ('component visible', is_component_visible),
                 'cp': ('compressed', is_compressed),
                 'ns': ('nearly stable', is_nearly_stable),
                 'tb': ('tree based', is_tree_based)
                }

  # report about a property
  def report_property(self, short):
    long = self.short_to_fn[short][0]
    val = self.short_to_fn[short][1](self)
    self.log("=== N is " + vocalize[val] + long + " ===")
    return ";" + vocalize[val] + short



def main():
    """The main part of the program."""
    # build argument parser and parse arguments -------------------------------
    parser = argparse.ArgumentParser(
        prog='RecoPhyNC',
        description='Recognizes phylogenetic network classes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-f','--file', type=str, help='the data file to analyze', dest='file'
    )
    parser.add_argument(
        '-v','--verbose', action='store_true', help="verbose mode", dest='verbose'
    )
    parser.add_argument(
        '-r','--rand', help="create random binary network. --rand x,y creates a binary network\
            with x nodes whose number of reticulations is picked from an exponential distribution \
            with mean y (default: y = x/10)", dest='rand'
    )
    parser.add_argument(
        '-o,--out', help="write graph to file in edgelist format", dest='out'
    )

    arguments = parser.parse_args()

    if not arguments.file and not arguments.rand:
      print "invalid arguments: need either --file or --rand, check --help for information"
      return

    folder = os.path.abspath(os.path.dirname(arguments.file if arguments.file else sys.argv[0]))
    
    # read all data files in specified folder and classify networks -----------
    with open(os.path.join(folder, "results2.csv"),"a") as output:
      PN = PhyloNetwork(arguments.verbose)
      # network initialization and preprocessing

      line = ''
      if arguments.file:
        PN.open(arguments.file)
        line = arguments.file
      elif arguments.rand:
        rand_args = [int(x) for x in arguments.rand.split(',')]
        PN.create_binary(rand_args[0], rand_args[1] if len(rand_args) > 1 else rand_args[0]/10)
        line = "(random-" + str(PN.N.number_of_nodes()) + "-" + str(PN.N.number_of_edges()) + ")"

      if arguments.out:
        with open(os.path.join(folder, arguments.out + ".graph"), "w") as graph_out:
          for uv in PN.N.edges():
            graph_out.write(str(uv[0]) + ' ' + str(uv[1]) + "\n")

      for short in ['ret', 'lvl', 'rcl', 'uc', 'ucb']:
        line += PN.report_value(short)

      for short in network_types:
        line += PN.report_property(short)

      # write information about the network
      output.write(line + "\n")


if __name__ == '__main__':
    main()
