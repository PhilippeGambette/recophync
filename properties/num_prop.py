
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
      stable_on_r = network.dominators[r]
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
            p = network.dominators[u]
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
class UncommonAncestors(NumericalNetworkProperty):

  # return a reverse topological order from r
  def rev_top_order_from(self, r):
    def explore(r):
      if r not in seen:
        seen.add(r)
        order.append(r)
        for j in self.network.N.predecessors(r):
          explore(j)

    order = []
    seen = set()
    explore(r)
    return order

  #Uncommon ancestors of parents of a reticulation node
  #Input  - reticulation node r, network G, reversed network R
  #Output - list of uncommon ancestors of parents of r
  def UA(self, r):
    #return the bitwise OR of all items in the list
    def list_or(L): return reduce(lambda a,b: a|b, L, 0)

    # get the reversed topological order from r upwards
    T = self.rev_top_order_from(r)
    [u1,u2] = self.network.N.predecessors(r)
    if T[0] != u1:
      u1, u2 = u2, u1
    L = dict()
    L[u1] = 1
    T.remove(u1)
    while T:
      q = T.pop(0)
      L[q] = list_or(map(lambda a:L[a], filter(lambda x: x in L, self.network.N.successors(q))))
      if q == u2:
        L[q] |= 2
    result = list(filter(lambda i: L[i] != 3, L.keys()))
    if result:
      self.log("uncommon ancestors for %s: %s" % (str(r), str(result)))
    return result


  # output: maximum number of uncommon ancestors of reticulations in N
  def check(self):
    if self.network.reticulations:
      self.val = max([len(self.UA(u)) for u in self.network.reticulations])
    else:
      self.val = 0




# the max over all reticulations r of the distance of its parents in N-r
class ShortestPath(NumericalNetworkProperty):

  # assign distances from u (offset by "current") in the BFS-tree T
  def assign_distances_tree(self, T, u, distances, current = 0):
    distances[u].append(current)
    for v in T.successors_iter(u):
      self.assign_distances_tree(T, v, distances, current + 1)


  # for a given u, return a dict mapping common ancestors x of u's parents to lists of distances to said parents
  def common_parent_distances(self, u):
    dists = defaultdict(list)
    for v in self.network.N.predecessors_iter(u):
      self.assign_distances_tree(networkx.bfs_tree(self.network.N, v, reverse = True), v, dists)
    # get a list of nodes who have distances from all predecessors
    # finally, filter non-pareto minimal items (which are not lowest ancestors)
    print("filtered: %s" % str(filter(lambda x: len(x) == self.network.N.in_degree(u), dists.values())))
    return pareto_mins(filter(lambda x: len(x) == self.network.N.in_degree(u), dists.values()))


  def check(self):
    # set the small and big reticulation heights at the same time
    srh_prop = self.network.properties['srh']
    brh_prop = self.network.properties['brh']

    self.val = 0
    srh_prop.val = 0
    brh_prop.val = 0

    for u in self.network.N:
      if self.network.N.in_degree(u) > 1:
        dists = self.common_parent_distances(u)
        # update val with the min over all common ancestors of their distances to parents of u
        self.val = max(self.val, min(sum(x) for x in dists))
        # update small reticulation height with the min over all common ancestors of the min of the distances to parents of u
        srh_prop.val = max(srh_prop.val, min(min(x) for x in dists))
        # update big reticulation height with the max over all common ancestors of the max of the distances to parents of u
        brh_prop.val = max(brh_prop.val, max(max(x) for x in dists))
       


# the max over all reticulations r of the smallest distance between r and any lowest common ancestor of r's parents
class SmallReticulationHeight(NumericalNetworkProperty):
  def check(self):
    self.network.properties['sp'].check()


# the max over all reticulations r of the farthest distance between r and any lowest common ancestor of r's parents
class BigReticulationHeight(NumericalNetworkProperty):
    def check(self):
      self.network.properties['sp'].check()



