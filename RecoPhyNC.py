#!/usr/bin/python

import argparse
import os
import sys
import re
import networkx

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
                         'ns' : ['cv']}
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

    X = set(self.N.nodes()) # TODO: use a FIFO
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
    X = set(self.leaves) # TODO: use FIFO
    # the current number of vertices that can see x
    spread = dict()
    for x in self.leaves:
      spread[x] = 1
    # for each vertex, the set of leaves it can reach
    reach = dict([(x,set()) for x in self.N.nodes()])
    reach.update(dict([(x,set(x)) for x in self.leaves]))
    # for each vertex, how many successors we have processed
    childs_seen = dict()

    while X:
      #TODO: make X a FIFO instead of a set
      i = X.pop()

      for j in reach[i]:
        spread[j] -= 1
      self.log('considering ' + i + ' with reachability ' + str(reach[i]))
      for p in self.N.predecessors(i):
        childs_seen[p] = childs_seen.get(p, 0) + 1
        for j in reach[i]:
          if j not in reach[p]:
            spread[j] = spread.get(j, 0) + 1
          reach[p].add(j)
        if childs_seen[p] == self.N.out_degree(p):
          X.add(p)
          self.stability[p] = set([x for x in self.leaves if x in reach[p] and spread[x] == 1])
      self.log('among ' + str(X) + ', the spread is ' + str(spread))
    #print "Time Stable vertices: "+str((datetime.datetime.now()-t0).microseconds)+"ms."
    
    if self.verbose:
      self.log('Stable vertices: ' + str([x for x in self.N.nodes() if self.is_stable(x) ]))

      

  # Input: none
  # Effect: set root to the first vertex of in-degree 0 and set regularity and reticulations and stable vertices
  def _update_infos(self):
    self.root = None
    self.regular = None
    for i in self.N:
      indeg = self.N.in_degree(i)
      outdeg = self.N.out_degree(i)
      if indeg > 2 or outdeg > 2:
        self.binary = False
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
    
    self.log('N has ' + str(self.N.number_of_nodes()) + " vertices and " \
        + str(self.N.number_of_edges()) + " edges." \
        + " Its root is " + str(self.root) + "."
        + "N is " + ("not " if self.regular == 0 else str(self.regular) + "-") + "regular")
    self.log('reticulations: ' + str(self.reticulations))


  # Input: text file containing a list of arcs
  # Output: network given as a dict associating to each vertex the table of its children
  def open(self, filename):
    with open(filename) as fd:
      for line in fd.readlines():  # match patterns of type "vertex1 vertex2"
        res = re.search("^([^ ]*) ([^ \n\r]*)$", line)
        if res:
          self.N.add_edge(res.group(1), res.group(2))
    self._contract()
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
        n = self.N.number_of_nodes()
        m = self.N.number_of_edges()
        dr = self.N.out_degree(self.root)
        retis = (m - n + 1) / (self.regular - 1)
        tree = (m - retis - dr) / self.regular
        leaf = n - retis - tree - 1
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
    unstable_roots = [u for u in self.reticulations if not self.is_reticulation(self.N.successors(u)[0]) and not self.is_stable(u)]
    self.log("unstable roots: " + str(unstable_roots))
    return unstable_roots

  # output: the number of unstable component roots in N
  def count_unstable_components(self):
    return len(self.list_unstable_components())

  # input: a list of nodes
  # output: the size of the largest sublist that intersects a block in N
  def max_per_block(self, v_set):
    biconn_comp = networkx.biconnected_components(self.N.to_undirected())
    lists_per_block = [ [x for x in self.N.subgraph(B) if x in v_set] for B in biconn_comp ]
    self.log('per block: ' + str(lists_per_block))
    return max([len(x) for x in lists_per_block])

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
    X = set([self.root])

    while X:
      v = X.pop()
      for s in self.N.successors(v):
        if self.is_reticulation(s):
          collective_indeg[s] = collective_indeg.get(s,0) \
                              + (collective_indeg[v] if self.is_reticulation(v) else 1)
        parents_seen[s] = parents_seen.get(s, 0) + 1
        if parents_seen[s] == self.N.in_degree(s):
          X.add(s)
    return max(collective_indeg.values())


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
        self.set_property(impl, True)
    else:
      for impl in negative_implications.get(prop, []):
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


  # Input: rooted network N
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
        self.log("Testing if " + str(v) + " is stable or all its parents are.")
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

    has_TC = dict([(x,True) for v in self.leaves for x in self.N.predecessors(v)[0]])
    
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
        'file', type=str,
        help='the data file to analyze'
    )
    parser.add_argument(
        '--verbose', action='store_true', help="verbose mode"
    )

    arguments = parser.parse_args()

    folder = os.path.abspath(os.path.dirname(arguments.file))

    # read all data files in specified folder and classify networks -----------
    with open(os.path.join(folder, "results2.csv"),"a") as output:
      PN = PhyloNetwork(arguments.verbose)
      # network initialization and preprocessing
      PN.open(arguments.file)

      line = arguments.file

      for short in ['ret', 'lvl', 'rcl', 'uc', 'ucb']:
        line += PN.report_value(short)

      for short in network_types:
        line += PN.report_property(short)

      # write information about the network
      output.write(line + "\n")


if __name__ == '__main__':
    main()
