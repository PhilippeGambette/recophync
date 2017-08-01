#!/usr/bin/env python

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


  # report about a property
  def report(self):
    self.log("=== N is " + vocalize[self.get()] + self.long + " ===")
    return ";" + vocalize[self.get()] + self.short



class NumericalNetworkProperty(NetworkProperty):
  # report about a value
  def report(self):
    self.log("=== " + self.long + ": " + str(self.get()) + " ===")
    return ";" + str(self.get())   



# ======================= boolean properties ===========================


class TreeBasedProp(NetworkProperty):
  # Input: none
  # Output: True if N is tree-based
  # NOTE: we use the characterization of Francis, Semple, Steel, Arxiv'16
  from networkx.algorithms import bipartite
  def check(self):
    network = self.network
    N = network.N
    # step 1: build G_N
    B = networkx.Graph()
    B.add_nodes_from([str(x)+'a' for x in N.nodes()], bipartite=0)
    B.add_nodes_from([str(x)+'b' for x in N.nodes()], bipartite=1)
    B.add_edges_from([(str(x)+'a', str(y)+'b') for (x,y) in N.edges()])

    # step 2: check if G_N has a matching of size |V| - |X|
    match = networkx.bipartite.maximum_matching(B)
    self.log("found matching of size " + str(len(match)/2) + ": " + str(match) + " (" + str(N.number_of_nodes() - len(network.leaves)) + " necessary for tree-based)")

    self.set(len(match) / 2 == N.number_of_nodes() - len(network.leaves))


class TreeChildProp(NetworkProperty):
  # Input: rooted network N
  # Output: True if N is tree-child, False otherwise
  def check(self):
    network = self.network
    successors = network.N.successors_iter
    for nonleaf in (v for v, degree in network.N.out_degree_iter() if degree):
        treechildNode = any(network.N.in_degree(c) == 1 for c in successors(nonleaf))
        if not treechildNode:
            self.log(str(nonleaf) + " is not treechild.")
            self.unset()
            if not network.verbose:
              return False
    self.set_if_none()


class NestedProp(NetworkProperty):
  # Output: True if N is nested, False otherwise
  def check(self):
    pass
    #TODO: write me


class GalledTreeProp(NetworkProperty):
  # Output: True if N is a galled tree, False otherwise
  def check(self):
    self.set(self.network.compute_level() == 1)


class TreeSiblingProp(NetworkProperty):
  # Output: True if N is tree-child, False otherwise
  def check(self):
    network = self.network
    N = network.N
    predecessors = N.predecessors_iter
    successors = N.successors_iter
    for v in network.reticulations:
        # for each reticulation vertex, we build its set of siblings S
        tree_siblings = [c for p in predecessors(v) for c in successors(p) if not network.is_reticulation(c)]
        # we check if S contains at least one tree vertex
        if tree_siblings:
          self.log(str(v) + " has tree siblings: " + str(tree_siblings))
        else:
          self.log(str(v) + " is not tree-sibling.")
          self.unset()
          if not network.verbose:
            return False
    self.set_if_none()


class ReticulationVisibleProp(NetworkProperty):
  # Input: a rooted phylogenetic network N
  # Output: True if N is reticulation visible, False otherwise
  def check(self):
    network = self.network
    for v in network.reticulations:
      if not network.is_stable(v):
        self.log(str(v) + " is not stable.")
        self.unset()
        if not network.verbose:
          return False
      else:
        self.log(str(v) + " is stable on " + str(network.stability.get(v, set())))
    self.set_if_none()


class NearlyStableProp(NetworkProperty):
  # Input: a rooted phylogenetic network N
  # Output: True if N is reticulation visible, False otherwise
  def check(self):
    network = self.network
    N = network.N
    predecessors = N.predecessors_iter
    for v in N:
      if v != network.root:
        if not network.is_stable(v):
          self.log(str(v) + " is not stable.")

          unstable_parents = [p for p in predecessors(v) if not network.is_stable(p)]

          if unstable_parents:
            self.log("Unstable parents of " + str(v) + ": " + str(unstable_parents))
            self.unset()
            if not network.verbose:
              return False
          else:
            self.log("All parents of " + str(v) + " are stable.")
        else:
          self.log(str(v) + " is stable.")
    self.set_if_none()


class ComponentVisibleProp(NetworkProperty):
  # Input: a rooted phylogenetic network N
  # Output: True if N is component visible, False otherwise
  def check(self):
    network = self.network
    successors = network.N.successors_iter
    for v in network.reticulations:
      s = successors(v)[0]
      if not network.is_reticulation(s) and not network.is_stable(s):
        self.log(str(s) + " is an unstable component root.")
        self.unset()
        if not network.verbose:
          return False
      else:
        self.log(str(s) + " is stable on " + str(network.stability.get(s, [])))
    self.set_if_none()


class CompressedProp(NetworkProperty):
  # Input: a rooted phylogenetic network N
  # Output: True if N is compressed, false otherwise
  def check(self):
    network = self.network
    predecessors = network.N.predecessors_iter
    for r in network.reticulations:
      for p in predecessors(r):
        if network.is_reticulation(p):
          self.log("Not compressed: parent " + str(p) + \
                   " of reticulation " + str(r) + \
                   " is also a reticulation.")
          self.unset()
          if not network.verbose:
            return False
    self.set_if_none()


class NearlyTreeChildProp(NetworkProperty):
  # input: a vertex and a dict to True/False
  # output: True <=> v has a tree path in N (saving everyone between in the dict)
  def check_TP(self, v, has_TP):
    if not self.network.is_reticulation(v):
      if v in has_TP:
        return has_TP[v]
      for i in self.network.N.successors(v):
        has_TP[v] = self.check_TP(i, has_TP)
        if has_TP[v]:
          return True
    return False

  # Input: a rooted phylogenetic network N
  # Output: True if N is nearly tree-child, false otherwise
  def check(self):
    network = self.network
    predecessors = network.N.predecessors_iter
    has_TP = dict([(next(predecessors(v)),True) for v in network.leaves])
    
    for r in network.reticulations:
      for p in predecessors(r):
        if not network.is_reticulation(p):
          if self.check_TP(p, has_TP):
            self.log(str(p) + " has a tree path.")
          else:
            self.log(str(p) + " has no tree path.")
            self.unset()
            if not network.verbose:
              return False
    self.set_if_none()


class GeneticallyStableProp(NetworkProperty):
  # Input: none
  # Output: True if N is genetically stable, False otherwise
  def check(self):
    network = self.network
    predecessors = network.N.predecessors_iter
    for r in network.reticulations:
      self.log("Testing if " + str(r) + " is stable.")
      if network.is_stable(r):
        self.log(str(r) + " is stable, testing if its parents are stable.")
        stable_parents = [x for x in predecessors(r) if network.is_stable(x)]
        self.log("Stable parents of " + str(r) + ": " + str(stable_parents))
        if len(stable_parents) == 0:
          self.unset()
          if not network.verbose:
            return False
      else:
        self.unset()
        if not network.verbose:
          return False
    self.set_if_none()




# ======================= numerical properties ===========================

class NumReticulations(NumericalNetworkProperty):
  def check(self):
    self.set(len(self.network.reticulations))


class Level(NumericalNetworkProperty):
  def check(self):
    self.set(self.network.max_per_block(self.network.reticulations))


class NumUnstableRoots(NumericalNetworkProperty):
  unstable_roots = None

  # output: the list of unstable component roots in N
  def compute_unstable_roots(self):
    network = self.network
    successors = network.N.successors_iter

    self.unstable_roots = set()
    for u in network.reticulations:
      if not network.is_stable(u):
        v = next(successors(u))
        if not network.is_reticulation(v):
          self.unstable_roots.add(v)
    # if there are unstable roots, we cannot have component visibility
    self.log("unstable roots: " + str(self.unstable_roots))
    network.properties['cv'].set(not self.unstable_roots)
    self.set(not self.unstable_roots)

  # output: the number of unstable component roots in N
  def check(self):
    if self.unstable_roots is None:
      self.compute_unstable_roots()
    self.set(len(self.unstable_roots))


class NumUnstableRootsPerBlock(NumericalNetworkProperty):
  def check(self):
    # compute unstable component roots
    network = self.network
    uc = network.properties['ur']
    
    uc.check()
    self.set(network.max_per_block(uc.unstable_roots))


class NestedDepth(NumericalNetworkProperty):
  # Input: binary rooted network N
  # Output: -1 if not nested, nested depth otherwise
  # TODO: Does not work yet
  def check(self):
      """Returns the nested depth of a rooted binary network N, or -1 if N is not
       nested.
      """
      biconn_comp = networkx.biconnected_components(self.network.N.to_undirected())
      for B in biconn_comp:
          # compute the nested depth of each biconnected component
          for node in B:
              print node
              #...


# the highest number of incoming arcs in N into any connected component of N[R]
# or 0 if R (set of reticulations) is empty
class MaxReticulationSubgraphIndegree(NumericalNetworkProperty):
  def check(self):
    # contract all edges between reticulations
    network = self.network
    successors = network.N.successors_iter

    parents_seen = dict()
    collective_indeg = dict()
    X = Queue.Queue()
    X.put(network.root)

    while not X.empty():
      v = X.get_nowait()
      for s in successors(v):
        if network.is_reticulation(s):
          collective_indeg[s] = collective_indeg.get(s,0) \
                              + (collective_indeg[v] if network.is_reticulation(v) else 1)
        parents_seen[s] = parents_seen.get(s, 0) + 1
        if parents_seen[s] == network.N.in_degree(s):
          X.put(s)
    return max(collective_indeg.values() or [0])







# ======================= main class ===========================


class PhyloNetwork:
  # whether or not verbose output is required
  verbose = None

  # the actual network
  N = None

  # the root of the network
  root = None

  # set of reticulations
  reticulations = set()

  # set of leaves
  leaves = set()

  # a dict mapping each vertex to the set of leaves it is stable on
  stability = dict()

  # network regularity (2 = binary, 0 = irregular)
  regular = None


  def __init__(self, _verbose = False):
    self.verbose = _verbose
    self.N = networkx.DiGraph()
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
                       'nd' : NestedDepth('nested depth', 'nd', self),
                       'ur' : NumUnstableRoots('#unstable component roots', 'ur', self),
                       'urb': NumUnstableRootsPerBlock('#unstable components roots per block', 'urb', self),
                       'rsi': MaxReticulationSubgraphIndegree('max #incoming edges to any reticulation component', 'rsi', self)
                      }


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


  # Input: some edge uv
  # Effect: contract uv in N unto v
  def contract_edge(self, uv, delete_node = True):
    u = uv[0]
    v = uv[1]
    self.N.remove_edge(u,v)
    self.N.add_edges_from([(v,w) for w in self.N.successors(u)])
    self.N.add_edges_from([(w,v) for w in self.N.predecessors(u)])
    if delete_node:
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
    for i in xrange((n - 1)/regularity):
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
      succ = self.N.successors(u)
      for w in succ:
        self.N.remove_edge(u, w)
        self.N.add_edge(v, w)
      self.N.add_edge(u, v)

    self._update_infos()
    self._update_stable()
     


  # ===================== computing values =========================

  # Input: a type t (0,1,2,3)
  # Output: list of vertices of N of type t
  def get_nodes(self, t):
    return [u for u in self.N.nodes() if self.node_type(u) == t]

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
        retis = (m - n + 1) / (self.regular - 1)
        tree = (m - retis) / self.regular
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
    lists_per_block = [ v_set.intersection(B) for B in biconn_comp if len(B) > 2]
    self.log('per block: ' + str(lists_per_block))
    return max([len(x) for x in lists_per_block]) if lists_per_block else 0

  
  def write_dot(self, filename):
    networkx.write_dot(N, filename)




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
        '-o,--out', help="write graph to file in graphviz dot format", dest='out'
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
        PN.write_dot(os.path.join(folder, arguments.out))

      for short in ['r', 'lvl', 'rsi', 'ur', 'urb'] + network_types:
        line += PN.properties[short].report()

      # write information about the network
      output.write(line + "\n")


if __name__ == '__main__':
    main()
