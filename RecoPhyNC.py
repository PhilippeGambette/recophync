#!/usr/bin/python3
import argparse
import os.path
import re
import networkx
import glob

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

VERBOSE = False

print("""
RecoPhyNC Copyright (C) 2015 Philippe Gambette
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions (GNU General Public License).
""")


# Input: text file containing a list of arcs
# Output: network given as a dict associating to each vertex the table of its children
def open_network(filename):
    """
    :rtype : networkx.DiGraph
    """
    network = networkx.DiGraph()
    with open(filename) as fd:
        for line in fd.readlines():  # match patterns of type "vertex1 vertex2"
            res = re.search("^([^ ]*) ([^ \n\r]*)$", line)
            if res:
                network.add_edge(res.group(1), res.group(2))
    return network


# Input: binary rooted network N
# Output: -1 if not nested, nested depth otherwise
# Does not work yet
def compute_nested_depth(network):
    """Returns the nested depth of a rooted binary network N, or -1 if N is not
     nested.
    :rtype : int
    """
    biconn_comp = networkx.biconnected_components(network.to_undirected())
    for B in biconn_comp:
        # compute the nested depth of each bi-connected component
        for node in B:
            print(node)
    return -1


# Input: binary rooted network N
# Output: level of N
def computeLevel(network):
    """
    :param network:
    :param :
    :return:
    """
    level = 0
    print("Computing level...")
    undirected = network.to_undirected()
    for BN, BE in zip(networkx.biconnected_components(undirected), networkx.biconnected_component_edges(undirected)):
        # compute the level of each bi-connected component
        m = len(BN)
        n = len(BE)
        if m - n + 1 > level:
            level = m - n + 1
    return level


# Input: rooted network N
# Output: True if N is tree-child, False otherwise
def isTreeChild(network):
    """
    :param network:
    :param :
    :return:
    """
    all_nodes_tree_child = True
    for non_leaf in (v for v, degree in network.out_degree() if degree):
        tree_child_node = any(network.in_degree(c) == 1 for c in network.successors(non_leaf))
        if not tree_child_node:
            if VERBOSE:
                print("Vertex", non_leaf, " is not tree-child.")
            all_nodes_tree_child = False
    return all_nodes_tree_child


# Input: rooted network N
# Output: True if N is tree-child, False otherwise
def isTreeSibling(network, reticulations_of_n):
    """
    :param network:
    :param reticulations_of_n:
    :return:
    """
    all_nodes_tree_sibling = True
    for v in reticulations_of_n:
        # for each reticulation vertex, we build its set of siblings S
        siblings = set(c for p in network.predecessors(v) for c in network.successors(p))
        siblings.discard(v)
        # we check if S contains at least one tree vertex
        if VERBOSE:
            print("Siblings of", v, ":", siblings)
        if not any(network.in_degree(s) == 1 for s in siblings):
            if VERBOSE:
                print("Vertex", v, "is not tree-sibling.")
            all_nodes_tree_sibling = False
    return all_nodes_tree_sibling


# Input: a rooted phylogenetic network N
# Output: True if N is reticulation visible, False otherwise
def isReticulationVisible(stable_vertices, reticulations_of_n):
    """
    :param stable_vertices:
    :param reticulations_of_n:
    :param :
    :return:
    """
    ret_visible = True
    for vertex in reticulations_of_n:
        if not ret_visible:
            break
        if VERBOSE:
            print("Testing if vertex", vertex, " is stable.")
        if vertex not in stable_vertices:
            ret_visible = False
            if VERBOSE:
                print("Vertex", vertex, "is not stable.")
        else:
            if VERBOSE:
                print("Vertex", vertex, "is stable.")
    return ret_visible


# Input: a rooted phylogenetic network N
# Output: True if N is reticulation visible, False otherwise
def isNearlyStable(network, stable_vertices):
    """
    :param network:
    :param stable_vertices:
    :return:
    """
    nearly_stable = True
    # skipping dots may improve speed
    for nonroot in (vertex for vertex, degree in network.in_degree() if degree):
        if not nearly_stable:
            break
        if VERBOSE:
            print("Testing if vertex", nonroot, "is stable or all its parents are.")
        if nonroot not in stable_vertices:
            if VERBOSE:
                print("Vertex", nonroot, "is not stable.")
            all_predecessors_stable = True
            for p in network.predecessors(nonroot):
                if p not in stable_vertices:
                    all_predecessors_stable = False
                    if VERBOSE:
                        print("Parent", p, "of vertex", nonroot, "is not stable.")
            if not all_predecessors_stable:
                nearly_stable = False
            else:
                if VERBOSE:
                    print("All parents of vertex", nonroot, "are stable.")
        else:
            if VERBOSE:
                print("Vertex", nonroot, "is stable.")
    return nearly_stable


# Input: a rooted phylogenetic network N rooted in N and a vertex v
# Output: True if v is stable in N, False otherwise
def isStable(v, network, network_root, network_leaves):
    if VERBOSE:
        print("Testing if vertex", v, "is stable.")

    if v == network_root:
        return True

    stable = False
    network_v = network.copy()
    network_v.remove_node(v)
    t = {vertex for vertex in networkx.dfs_tree(network_v, network_root).nodes()}
    for leaf in network_leaves.difference(t):
        if VERBOSE:
            print("Vertex", v, "is stable for", leaf)
        stable = True
    if not stable and VERBOSE:
        print("Vertex", v, "is not stable.")
    return stable


# Input: a rooted phylogenetic network N
# Output: True if N is compressed, false otherwise
def isCompressed(network, network_reticulations):
    """Returns True if network N is compressed, False otherwise."""
    # skipping dots may improve speed
    compressed = True
    for reticulation in network_reticulations:
        for p in (w for w in network.predecessors(reticulation) if network.in_degree(w) > 1):
            if VERBOSE:
                print("Not compressed: the parent", p,
                      "of reticulation vertex", reticulation,
                      "is also a reticulation vertex.")
            compressed = False
    return compressed


# Input: a rooted phylogenetic network N
# Output: True if N is nearly tree-child, false otherwise
def isNearlyTreeChild(network, stable_vertices, reticulations_of_n):
    """
    :param network:
    :param stable_vertices:
    :param reticulations_of_n:
    :return:
    """
    # First check if N is stable
    if not isReticulationVisible(stable_vertices, reticulations_of_n):
        return False

    nearly_tree_child = True
    # skipping dots may improve speed
    # let's filter reticulations directly
    for reticulation in (v for v, d in network.in_degree() if d > 1):
        if not nearly_tree_child:
            break
        # Check if at least one parent of v has the tree path property
        parent_with_tree_path = False
        for p in network.predecessors(reticulation):
            if not parent_with_tree_path:
                visited = set()
                if VERBOSE:
                    print("Testing if vertex " + p + " has the tree path property.")
                if hasTreePath(p, network, visited):
                    parent_with_tree_path = True
        if not parent_with_tree_path:
            nearly_tree_child = False
            if VERBOSE:
                print("Not nearly tree-child: no parent of reticulation vertex", reticulation,
                      "has the tree path property.")
    return nearly_tree_child


# Input: a vertex v and a rooted phylogenetic network N
# Output: True if v has the tree path property in N, else otherwise
def hasTreePath(v, network, visited):
    """
    :param v:
    :param network:
    :param visited:
    :return:
    """
    has_tp = not network.out_degree(v)
    # if v is a tree vertex then we check if one of its children has a tree path
    if v not in visited:
        visited.add(v)  # visited is now a set (see isNearlyTreeChild)
        n_in_degree = network.in_degree
        for tree_successor in (s for s in network.successors(v) if n_in_degree(s) == 1):
            has_tp = hasTreePath(tree_successor, network, visited)
            if has_tp:
                break

    if not has_tp and VERBOSE:
        print("Vertex", v, " does not have the tree path property.")
    return has_tp


# Input: a rooted phylogenetic network N
# Output: True if N is genetically stable, False otherwise
def isGeneticallyStable(network, stable_vertices, reticulations_of_n):
    gen_stable = True
    for reticulation in reticulations_of_n:
        if not gen_stable:
            break
        if VERBOSE:
            print("Testing if vertex", reticulation, "is stable.")
        if reticulation not in stable_vertices:
            gen_stable = False
        else:
            if VERBOSE:
                print("Vertex", reticulation, "is stable, testing if its parents are stable.")
            one_parent_stable = False
            # let's filter stable predecessors directly
            for stable_predecessor in (v for v in network.predecessors(reticulation) if v in stable_vertices):
                one_parent_stable = True
                if VERBOSE:
                    print("Parent", stable_predecessor, "is stable.")
            if not one_parent_stable:
                gen_stable = False
    return gen_stable


# Input: a rooted phylogenetic network N
# Output: N's root
def root(network):
    """
    :param network:
    :return:
    """
    for vertex, degree in network.in_degree():
        if not degree:
            return vertex


# Input: a rooted phylogenetic network N
# Output: N where each in-degree and out-degree 1 vertex is contracted with its parent
def contract(network):
    """
    :param network:
    :return:
    """
    found_candidate = True
    n_in_degree = network.in_degree()
    n_out_degree = network.out_degree()
    n_nodes = list(network.nodes())
    while found_candidate:
        found_candidate = False
        for vertex in (v for v in n_nodes if n_in_degree(v) == 1 and n_out_degree(v) == 1):
            if VERBOSE:
                print(vertex, "was an in-degree 1 out-degree 1 vertex")
            network.add_edge(list(network.predecessors(vertex))[0], list(network.successors(vertex))[0])
            network.remove_node(vertex)
            found_candidate = True


# Input: a rooted phylogenetic network N
# Output: number of leaves of N
def leafNumber(network):
    """
    :param network:
    :return:
    """
    return sum(1 for vertex, degree in network.out_degree() if not degree)


# Input: a rooted phylogenetic network N
# Output: set of leaves of N
def leaves(network):
    """
    :param network:
    :return:
    """
    return {vertex for vertex, degree in network.out_degree() if not degree}


def reticulations(network):
    """Returns the set of reticulations in N."""
    return {vertex for vertex, degree in network.in_degree() if degree > 1}


def main():
    """The main part of the program."""
    # build argument parser and parse arguments -------------------------------
    global VERBOSE
    parser = argparse.ArgumentParser(
        prog='RecoPhyNC',
        description='Recognizes phylogenetic network classes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-f', type=str, dest='folder',
        help='the path to the folder that contains the data files'
    )
    parser.add_argument(
        '--verbose', action='store_true', help="verbose mode"
    )

    arguments = parser.parse_args()
    VERBOSE = arguments.verbose

    folder = os.path.abspath(os.path.dirname(arguments.folder))
    # Filter used for the files contained in the data folder:
    file_filter = "*"

    # read all data files in specified folder and classify networks -----------
    with open(os.path.join(folder, "results2.csv"), "w") as output:
        output.write('file,level,tree-child,nearly tree-child,genetically stable,tree-sibling,reticulation-visible,'
                     'compressed,nearly stable\n')
        for data_file in glob.glob(os.path.join(folder, "data", file_filter)):
            # network initialization and preprocessing
            print("Treating file", data_file)
            network = open_network(data_file)
            contract(network)

            # store identified classes and useful properties in a dictionary,
            # so we can use previously acquired knowledge to speed things up
            # (i.e., if we know that network is tree-child, then we know it's also
            # nearly tree-child so no need to ask a function for the result)
            n_classification = dict()

            # compute stable vertices
            if VERBOSE:
                print("""== Computing the stable vertices ==""")
            # t0=datetime.datetime.now()
            # knowing the root, the leaves and the reticulations is useful at
            # various stages, so let's compute all of that only once
            root_n = root(network)
            leaves_of_n = leaves(network)
            reticulations_of_n = reticulations(network)
            stable_vertices = {
                v for v in network.nodes() if
                not (network.out_degree and network.in_degree and not isStable(v, network, root_n, leaves_of_n))
            }

            if VERBOSE:
                print("Stable vertices of network:", stable_vertices)
                print("""===================================""")
            # print "Time Stable vertices: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            line = data_file
            if VERBOSE:
                print("Edges of network: ", network.edges())

            # level
            # t0 = datetime.datetime.now()
            level = computeLevel(network)
            if VERBOSE:
                print("Level of network:", level)
                print("""""")
            n_classification["level"] = level
            line = ''.join((line, ",", str(level)))  # join is faster
            # print "Time Level: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # tree-child
            # t0 = datetime.datetime.now()
            n_classification["tc"] = n_classification.get("level", -1) == 1 or isTreeChild(network)
            if n_classification["tc"]:
                if VERBOSE:
                    print("""== network is tree-child.==""")
                line = ''.join((line, ",tc"))
            else:
                if VERBOSE:
                    print("""== network is not tree-child. ==""")
                line = ''.join((line, ",not tc"))
            # print "Time Tree-child: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # nearly tree-child
            # t0=datetime.datetime.now()
            ntc = n_classification.get("tc", False) or isNearlyTreeChild(network, stable_vertices, reticulations_of_n)
            n_classification["ntc"] = ntc
            if VERBOSE:
                if ntc:
                    print("""== network is nearly tree-child. ==""")
                else:
                    print("""== network is not nearly tree-child. ==""")
            if ntc:
                line = ''.join((line, ",ntc"))
            else:
                line = ''.join((line, ",not ntc"))
            # print "Time Nearly tree-child: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # genetically stable
            # t0=datetime.datetime.now()
            gen_stab = n_classification.get("ntc", False) or isGeneticallyStable(network, stable_vertices,
                                                                                 reticulations_of_n)
            n_classification["gs"] = gen_stab
            if VERBOSE:
                if gen_stab:
                    print("""== network is genetically stable. ==""")
                else:
                    print("""== network is not genetically stable. ==""")
            if gen_stab:
                line = ''.join((line, ",gs"))
            else:
                line = ''.join((line, ",not gs"))
            # print "Time Genetically stable: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # tree-sibling
            # t0=datetime.datetime.now()

            tr_sib = n_classification.get("gs", False) or isTreeSibling(network, reticulations_of_n)
            n_classification["ts"] = tr_sib
            if VERBOSE:
                if tr_sib:
                    print("""== network is tree-sibling. ==""")
                else:
                    print("""== network is not tree-sibling. ==""")
            if tr_sib:
                line = ''.join((line, ",ts"))
            else:
                line = ''.join((line, ",not ts"))

            # print "Time Tree-sibling: " +str((datetime.datetime.now()-t0).microseconds) + "ms."

            # reticulation-visible
            # t0 = datetime.datetime.now()
            ret_vis = n_classification.get("gs", False) or isReticulationVisible(
                stable_vertices, reticulations_of_n
            )
            n_classification["rv"] = ret_vis
            if VERBOSE:
                if ret_vis:
                    print("""== network is reticulation visible. ==""")
                else:
                    print("""== network is not reticulation visible. ==""")
            if ret_vis:
                line = ''.join((line, ",rv"))
            else:
                line = ''.join((line, ",not rv"))
            # print "Time Reticulation-visible: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # compressed
            # t0 = datetime.datetime.now()
            comp = isCompressed(network, reticulations_of_n)
            n_classification["cp"] = comp
            if VERBOSE:
                if comp:
                    print("""== network is compressed. ==""")
                else:
                    print("""== network is not compressed. ==""")
            if comp:
                line = ''.join((line, ",cp"))
            else:
                line = ''.join((line, ",not cp"))
            # print "Time Compressed: "+str((datetime.datetime.now() - t0).microseconds)+"ms."

            # nearly-stable
            # t0 = datetime.datetime.now()
            ns = n_classification.get("tc", False) or isNearlyStable(network, stable_vertices)
            n_classification["ns"] = ns
            if VERBOSE:
                if ns:
                    print("""== network is nearly stable. ==""")
                else:
                    print("""== network is not nearly stable. ==""")
            if ns:
                line = ''.join((line, ",ns"))
            else:
                line = ''.join((line, ",not ns"))
            # print "Time Compressed: " + str((datetime.datetime.now() - t0).microseconds)+"ms."

            # write information about the network
            output.write(line + "\n")


if __name__ == '__main__':
    main()
