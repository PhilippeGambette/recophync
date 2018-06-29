
from queue import Queue
import networkx
from collections import defaultdict
from common import *

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
          elif network.is_reticulation(u):
            # if we arrive at a reticulation, then continue from the high-point of its cycle
            p = network.stability_tree.predecessors(u)[0]
            self.log('jumping from ' + str(u) + ' to ' + str(p))
            # add r -> u to the nesting tree
            self.nesting_tree.add_edge(r, u)
          else:
            #u in taboo_nodes:
            # if we climbed or jumped to a previous jump target we know that N is NOT NESTED
            self.log('encountered taboo node ' + str(u) + ' between ' + str(stable_on_r) + ' and ' + str(r))
            self.set(-1)
            return


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


# maximum number of uncommon ancestors in N
class NumUncommonAncestors(NumericalNetworkProperty):

  # output: maximum number of uncommon ancestors of reticulations in N
  def check(self):
      pass




# the max over all reticulations r of the distance of its parents in N-r
class ShortestPath(NumericalNetworkProperty):
    def check(self):
        pass




# the max over all reticulations r of the smallest distance between r and any lowest common ancestor of r's parents
class SmallReticulationHeight(NumericalNetworkProperty):
    def check(self):
        pass
#===================== CONTINUE HERE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                






# the max over all reticulations r of the farthest distance between r and any lowest common ancestor of r's parents
class BigReticulationHeight(NumericalNetworkProperty):
    def check(self):
        pass



