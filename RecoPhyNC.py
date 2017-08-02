#!/usr/bin/env python

#Python3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import os
import re
import networkx
from queue import Queue
import sys
import random
import glob
from collections import defaultdict

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

print("""
RecoPhyNC Copyright (C) 2015 Philippe Gambette
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions (GNU General Public License).
""")

# a dict that gives "not " if the input is False and "" otherwise
vocalize = {False: "not ", True: ""}

numerical_properties = ['r', 'lvl', 'rsi', 'ur', 'urb', 'nd']
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

  # Python3 compat
  __bool__ = __nonzero__

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
    B.add_nodes_from((str(x)+'a' for x in N.nodes_iter()), bipartite=0)
    B.add_nodes_from((str(x)+'b' for x in N.nodes_iter()), bipartite=1)
    B.add_edges_from(((str(x)+'a', str(y)+'b') for (x,y) in N.edges_iter()))

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
        self.log(str(v) + " is stable on " + str(network.stability.get(v, None)))
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
        self.log(str(s) + " is stable on " + str(network.stability.get(s, None)))
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
      for i in self.network.N.successors_iter(v):
        has_TP[v] = self.check_TP(i, has_TP)
        if has_TP[v]:
          return True
    return False

  # Input: a rooted phylogenetic network N
  # Output: True if N is nearly tree-child, false otherwise
  def check(self):
    network = self.network
    predecessors = network.N.predecessors_iter
    has_TP = dict((next(predecessors(v)),True) for v in network.leaves)
    
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


class NestingDepth(NumericalNetworkProperty):
  # a rooted tree of reticulations in N that is the transitive reduction of the graph with
  # (u,v) is an edge <=> v and lowest_stable[v] are strictly between u & lowest_stable[u]
  # NOTE: the depth of this tree (if it exists) is the nesting depth
  nesting_tree = None

  # go upwards in the network until we see a reticulation or a blocker
  # return the blocker/reticulation (or None if the root was reached),
  #     as well as a set of vertices on the path with their distance to u
  def climb_tree_until(self, u, blocker, taboo_nodes):
    while not self.network.is_reticulation(u) and not u == blocker and not u in taboo_nodes:
      taboo_nodes.add(u)
      u = self.network.N.predecessors(u)[0]
    return u

  # travel from a reticulation to its lowest stable vertex w, jumping over nested blocks 
  # if ever we encounter a node that we've seen in a previous climb (taboo_nodes), we
  # can conclude that the network is not nested
  def goto_lowest_stable(self, r, taboo_nodes):
    if r not in self.lowest_stable:
      network = self.network
      predecessors = network.N.predecessors_iter
      # save the lowest stable vertex above r
      stable_on_r = network.stability_tree.predecessors(r)[0]
      for p in predecessors(r):
        # climb the tree from p, jumping over previously evaluated cycles using lowest_stable
        u = p
        while True:
          u = self.climb_tree_until(p, stable_on_r, taboo_nodes)
          self.log(str(r) + ': climbed from ' + str(p) + ' (' + network.node_type(p) + ') to '
                                   + str(u) + ' (' + network.node_type(u) + ')')

          if u == stable_on_r:
            break

          # if we climbed or jumped to a previous jump target we know that N is NOT NESTED
          if u in taboo_nodes:
            self.log('encountered taboo node ' + str(u) + ' between ' + str(stable_on_r) + ' and ' + str(r))
            self.set(-1)
            if not network.verbose:
              return None

          if network.is_reticulation(u):
            # if we arrive at a reticulation, then continue from the high-point of its cycle
            p = network.stability_tree.predecessors(u)[0]
            self.log('jumping from ' + str(u) + ' to ' + str(p))
            # add r -> u to the nesting tree
            self.nesting_tree.add_edge(r, u)
          else:
            # otherwise, we've found stable_on_r
            break


  # compute and return the depth in the nesting_tree
  def depth_in_nesting_tree(self, r, depths):
    if not r in depths:
      p = self.nesting_tree.predecessors(r)
      if p:
        if len(p) > 1:
          raise ValueError('Mathias is bad at programming')
        d = self.depth_in_nesting_tree(p[0], depths) + 1
        depths[r] = d
        return d
      else:
        # by definition of nesting depth, only trees have nesting depth 0,
        # so the reticulation at the root of the nesting tree gets nesting depth 1
        return 1
    else:
      return depths[r]

  # Input: binary rooted network N
  # Output: -1 if not nested, nested depth otherwise
  def check(self):
    """Returns the nested depth of a rooted binary network N, or -1 if N is not nested."""
    network = self.network
    # note that level-0/1 implies nesting depth 0/1
    level = network.properties['lvl'].get()
    if level > 1:
      predecessors = network.N.predecessors_iter
      self.lowest_stable = dict()
      for B in networkx.biconnected_components(network.N.to_undirected()):
        if len(B) > 2:
          BR = B.intersection(network.reticulations)
          self.nesting_tree = networkx.DiGraph()
          self.nesting_tree.add_nodes_from(BR)

          taboo_nodes = set()
          for r in BR:
            self.goto_lowest_stable(r, taboo_nodes)

          if self.val is None:
            self.log('nesting tree: ' + str(self.nesting_tree.edges()))
            depths = dict()
            for r in BR:
              self.depth_in_nesting_tree(r, depths)
            self.log('nesting tree depths: ' + str(depths))
            self.set(max(depths.values() or [1])) # NOTE: if BR is a single cycle, the nesting tree is empty
          elif not network.verbose:
            return
      # if we still didn't set the nesting depths, then |B|=2 for all blocks (N is a tree), so set it to 0
      self.set_if_none(0)
    else:
      self.set(level)




# the highest number of incoming arcs in N into any connected component of N[R]
# or 0 if R (set of reticulations) is empty
class MaxReticulationSubgraphIndegree(NumericalNetworkProperty):
  def check(self):
    # contract all edges between reticulations
    network = self.network
    successors = network.N.successors_iter

    parents_seen = defaultdict(int)
    collective_indeg = defaultdict(int)
    X = Queue()
    X.put(network.root)

    while not X.empty():
      v = X.get_nowait()
      for s in successors(v):
        if network.is_reticulation(s):
          collective_indeg[s] += (collective_indeg[v] if network.is_reticulation(v) else 1)
        parents_seen[s] += 1
        if parents_seen[s] == network.N.in_degree(s):
          X.put(s)
    self.set(max(collective_indeg.values() or [0]))







# ======================= main class ===========================


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

  # stability tree of N: rooted tree in which (u,v) <=> u is the lowest vertex in N that is stable on v
  stability_tree = None

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
    self.stability_tree = networkx.DiGraph()
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
                       'rsi': MaxReticulationSubgraphIndegree('max #incoming edges to any reticulation component', 'rsi', self)
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
    u = uv[0]
    v = uv[1]
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
    X = set(self.N.nodes()) #NOTE: we need a set here to avoid having a node multiple times in the queue
    while X:
      v = X.pop()
      #self.log('considering ' + str(v) + ' with indeg ' + str(self.N.in_degree(v)) + ' & outdeg ' + str(self.N.out_degree(v)))
      if self.N.in_degree(v) == self.N.out_degree(v) == 1:
        u = next(self.N.predecessors_iter(v))
        w = next(self.N.successors_iter(v))
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
    
    # initialize stability with the parents of all leaves
    self.stability.update((next(self.N.predecessors_iter(l)),l) for l in self.leaves)
    # the current bottom-up front
    X = Queue()
    map(X.put, self.leaves)

    # the current number of vertices that can see x
    spread = defaultdict(int)
    # for each vertex, the set of leaves it can reach
    reach = defaultdict(set)
    # for each vertex, how many successors we have processed
    childs_seen = defaultdict(int)

    while not X.empty():
      v = X.get(False)

      for j in reach[v]:
        spread[j] -= 1
      for p in self.N.predecessors_iter(v):
        childs_seen[p] = childs_seen.get(p, 0) + 1
        for j in reach[v].union(set([v])):
          if j not in reach[p]:
            spread[j] += 1
          reach[p].add(j)

      self.log(str(v) + ": spreads: " + str(spread))
      for p in self.N.predecessors_iter(v):
        if childs_seen[p] == self.N.out_degree(p):
          X.put(p)
          #self.stability[p] = set(x for x in reach[p] if spread[x] == 1)
          p_stable = [i for i in reach[p] if spread[i] == 1]
          self.log(str(p) + ' is stable on ' + str(p_stable))
          for x in p_stable:
            self.stability_tree.add_edge(p, x)
            # update stability dict (for leaves) - note that all parents of leaves are initialized
            if x in self.stability:
              self.stability[p] = self.stability[x]

            # forget about x
            reach[p].discard(x)
            del spread[x]

      del reach[v]
    #print "Time Stable vertices: "+str((datetime.datetime.now()-t0).microseconds)+"ms."
    
    if self.verbose:
      self.log('stability:')
      for u in self.stability:
        if self.is_stable(u):
          self.log(str(u) + ': ' + str(self.stability[u]))
      self.log('stability tree: ' + str(self.stability_tree.edges()))

      

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
        '-d','--dir', type=str, help='directory/folder of data files to analyze', dest='dir'
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
        '-o','--out', help="write graph to file in graphviz dot format; if --dir is given, then the data filename and '.dot' are prepended.", dest='out'
    )

    arguments = parser.parse_args()

    if sum(bool(x) for x in [arguments.file, arguments.dir, arguments.rand]) != 1:
      print("Invalid arguments! I need EXACTLY one of --file, --dir, and --rand. Check --help for more information.")
      return

    file_list = ''
    if arguments.dir:
      folder = os.path.abspath(os.path.dirname(arguments.dir))
      file_list = glob.glob(os.path.join(folder, "*"))
    elif arguments.file:
      folder = os.path.abspath(os.path.dirname(arguments.file))
      file_list = [arguments.file]
    else:
      folder = os.path.abspath(os.path.dirname(sys.argv[0]))
      file_list = [0]

    # read all data files in specified folder and classify networks -----------
    with open(os.path.join(folder, "results2.csv"),"a") as output:
      for data_file in file_list:
        PN = PhyloNetwork(arguments.verbose)

        line = ''
        if data_file:
          print("working on ", data_file)
          PN.open(data_file)
          line = os.path.basename(data_file)
        else:
          rand_args = [int(x) for x in arguments.rand.split(',')]
          PN.create_binary(rand_args[0], rand_args[1] if len(rand_args) > 1 else rand_args[0]/10.0)
          line = "(random-" + str(PN.N.number_of_nodes()) + "-" + str(PN.N.number_of_edges()) + ")"

        if arguments.out:
          PN.write_dot(os.path.join(folder, arguments.out + (os.path.basename(data_file) + '.dot' if arguments.dir else '') ))

        for short in numerical_properties + network_types:
          line += PN.properties[short].report()

        # write information about the network
        output.write(line + "\n")



if __name__ == '__main__':
    main()
