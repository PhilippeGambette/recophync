
import networkx
from common import *
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



