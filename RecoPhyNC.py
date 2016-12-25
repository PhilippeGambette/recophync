#!/usr/sfw/bin/python
import argparse
import os.path, re, networkx, glob

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

# Input: text file containing a list of arcs
# Output: network given as a dict associating to each vertex the table of its children
def open_network(filename):
    """

    :rtype : networkx.DiGraph
    """
    # fd = open(filename,"r")
    # lines = fd.readlines()
    # N=networkx.DiGraph()
    #
    # for line in lines:
    #     # Recognize pattern of type "vertex1 vertex2"
    #     res=re.search("^([^ ]*) ([^ \n\r]*)$",line)
    #     if res:
    #         N.add_edge(res.group(1),res.group(2))
    # return N
    N = networkx.DiGraph()
    with open(filename) as fd:
        for line in fd.readlines():  # match patterns of type "vertex1 vertex2"
            res = re.search("^([^ ]*) ([^ \n\r]*)$", line)
            if res:
                N.add_edge(res.group(1), res.group(2))
    return N


# Input: binary rooted network N
# Output: -1 if not nested, nested depth otherwise
# Does not work yet
def compute_nested_depth(N):
    """Returns the nested depth of a rooted binary network N, or -1 if N is not
     nested.

    :rtype : int
    """
    biconn_comp = networkx.biconnected_components(N.to_undirected())
    for B in biconn_comp:
        # compute the nested depth of each biconnected component
        for node in B:
            print node
            #...
    return -1


# Input: binary rooted network N
# Output: level of N
def computeLevel(N):
    """

    :param N:
    :param :
    :return:
    """
    level = 0
    for BN in networkx.biconnected_component_subgraphs(N.to_undirected()):
        # compute the level of each biconnected component
        m = BN.number_of_edges()
        n = BN.number_of_nodes()
        if m-n+1 > level:
            level = m-n+1
            #print "Biconnected component "+str(i)+": "+str(m)+" edges, "+str(n)+" nodes, cyclomatic number "+str(m-n+1)
    return level


# Input: rooted network N
# Output: True if N is tree-child, False otherwise
def isTreeChild(N):
    """

    :param N:
    :param :
    :return:
    """
    allNodesTreechild = True
    N_in_degree = N.in_degree
    N_successors_iter = N.successors_iter
    for nonleaf in (v for v, degree in N.out_degree_iter() if degree):
        treechildNode = any(N_in_degree(c) == 1 for c in N_successors_iter(nonleaf))
        if not treechildNode:
            if VERBOSE:
                print "Vertex", nonleaf, " is not treechild."
            allNodesTreechild = False
    return allNodesTreechild


# Input: rooted network N
# Output: True if N is tree-child, False otherwise
def isTreeSibling(N, reticulations_of_N):
    """

    :param N:
    :param reticulations_of_N:
    :param verbose:
    :return:
    """
    N_predecessors_iter = N.predecessors_iter
    N_successors_iter = N.successors_iter
    allNodesTreesibling = True
    for v in reticulations_of_N:
        # for each reticulation vertex, we build its set of siblings S
        S = set(
            c for p in N_predecessors_iter(v) for c in N_successors_iter(p)
        )
        S.discard(v)
        # we check if S contains at least one tree vertex
        if VERBOSE:
            print "Siblings of", v, ":", S
        if not any(N.in_degree(s) == 1 for s in S):
            if VERBOSE:
                print "Vertex", v, "is not tree-sibling."
            allNodesTreesibling = False
    return allNodesTreesibling

# Input: a rooted phylogenetic network N
# Output: True if N is reticulation visible, False otherwise
def isReticulationVisible(stableVertices, reticulations_of_N):
    """

    :param N:
    :param stableVertices:
    :param reticulations_of_N:
    :param :
    :return:
    """
    retVisible = True
    for vertex in reticulations_of_N:
        if not retVisible:
            break
        if VERBOSE:
            print "Testing if vertex", vertex, " is stable."
        if vertex not in stableVertices:
            retVisible = False
            if VERBOSE:
                print "Vertex", vertex, "is not stable."
        else:
            if VERBOSE:
                print "Vertex", vertex, "is stable."
    return retVisible



# Input: a rooted phylogenetic network N
# Output: True if N is reticulation visible, False otherwise
def isNearlyStable(N, stableVertices):
    """

    :param N:
    :param stableVertices:
    :param verbose:
    :return:
    """
    nearlyStable = True
    # skipping dots may improve speed
    N_predecessors_iter = N.predecessors_iter
    for nonroot in (vertex for vertex, degree in N.in_degree_iter() if degree):
        if not nearlyStable:
            break
        if VERBOSE:
            print "Testing if vertex", nonroot, "is stable or all its parents are."
        if nonroot not in stableVertices:
            if VERBOSE:
                print "Vertex", nonroot, "is not stable."
            allPredecessorsStable = True
            for p in N_predecessors_iter(nonroot):
                if p not in stableVertices:
                    allPredecessorsStable = False
                    if VERBOSE:
                        print "Parent", p, "of vertex", nonroot, "is not stable."
            if not allPredecessorsStable:
                nearlyStable = False
            else:
                if VERBOSE:
                    print "All parents of vertex", nonroot, "are stable."
        else:
            if VERBOSE:
                print "Vertex", nonroot, "is stable."
    return nearlyStable


# Input: a rooted phylogenetic network N rooted in N and a vertex v
# Output: True if v is stable in N, False otherwise
def isStable(v, N, network_root, network_leaves):
    if VERBOSE:
        print "Testing if vertex",  v, "is stable."

    if v == network_root:
        return True

    stable = False
    Nv = N.copy()
    Nv.remove_node(v)
    T = {vertex for vertex in networkx.dfs_tree(Nv, network_root).nodes_iter()}
    for leaf in network_leaves.difference(T):
        if VERBOSE:
            print "Vertex", v, "is stable for", leaf
        stable = True
    if not stable and VERBOSE:
        print "Vertex", v, "is not stable."
    return stable


# Input: a rooted phylogenetic network N
# Output: True if N is compressed, false otherwise
def isCompressed(N, network_reticulations):
    """Returns True if network N is compressed, False otherwise."""
    # skipping dots may improve speed
    N_in_degree = N.in_degree
    N_predecessors_iter = N.predecessors_iter
    compressed = True
    for reticulation in network_reticulations:
        for p in (w for w in N_predecessors_iter(reticulation) if N_in_degree(w) > 1):
            if VERBOSE:
                print "Not compressed: the parent", p, \
                      "of reticulation vertex", reticulation, \
                      "is also a reticulation vertex."
            compressed = False
    return compressed


# Input: a rooted phylogenetic network N
# Output: True if N is nearly tree-child, false otherwise
def isNearlyTreeChild(N, stableVertices, reticulations_of_N):
    """

    :param N:
    :param stableVertices:
    :param reticulations_of_N:
    :param verbose:
    :return:
    """
    # First check if N is stable
    if not isReticulationVisible(stableVertices, reticulations_of_N):
        return False

    nearlyTreeChild = True
    # skipping dots may improve speed
    N_predecessors_iter = N.predecessors_iter
    # let's filter reticulations directly
    for reticulation in (v for v, d in N.in_degree_iter() if d > 1):
        if not nearlyTreeChild:
            break
        # Check if at least one parent of v has the tree path property
        parentWithTreePath = False
        for p in N_predecessors_iter(reticulation):
            if not parentWithTreePath:
                visited = set()
                if VERBOSE:
                    print "Testing if vertex "+p+" has the tree path property."
                if hasTreePath(p, N, visited):
                    parentWithTreePath = True
        if not parentWithTreePath:
            nearlyTreeChild = False
            if VERBOSE:
                print "Not nearly tree-child: no parent of reticulation vertex", reticulation, "has the tree path property."
    return nearlyTreeChild


# Input: a vertex v and a rooted phylogenetic network N
# Output: True if v has the tree path property in N, else otherwise
def hasTreePath(v, N, visited):
    """

    :param v:
    :param N:
    :param visited:
    :param verbose:
    :return:
    """
    hasTP = not N.out_degree(v)
    # if v is a tree vertex then we check if one of its children has a tree path
    if v not in visited:
        visited.add(v)  # visited is now a set (see isNearlyTreeChild)
        N_in_degree = N.in_degree
        for tree_successor in (s for s in N.successors_iter(v) if N_in_degree(s)==1):
            hasTP = hasTreePath(tree_successor, N, visited)
            if hasTP:
                break

    if not hasTP and VERBOSE:
        print "Vertex", v, " does not have the tree path property."
    return hasTP


# Input: a rooted phylogenetic network N
# Output: True if N is genetically stable, False otherwise
def isGeneticallyStable(N, stableVertices, reticulations_of_N):
    genStable = True
    N_predecessors_iter = N.predecessors_iter
    for reticulation in reticulations_of_N:
        if not genStable:
            break
        if VERBOSE:
            print "Testing if vertex", reticulation, "is stable."
        if reticulation not in stableVertices:
            genStable = False
        else:
            if VERBOSE:
                print "Vertex", reticulation, "is stable, testing if its parents are stable."
            oneParentStable = False
            # let's filter stable predecessors directly
            for stable_predecessor in (v for v in N_predecessors_iter(reticulation) if v in stableVertices):
                oneParentStable = True
                if VERBOSE:
                    print "Parent", stable_predecessor, "is stable."
            if not oneParentStable:
                genStable = False
    return genStable


# Input: a rooted phylogenetic network N
# Output: N's root
def root(N):
    """

    :param N:
    :return:
    """
    for vertex, degree in N.in_degree_iter():
        if not degree:
            return vertex


# Input: a rooted phylogenetic network N
# Output: N where each indegree and outdegre 1 vertex is contracted with its parent
def contract(N):
    """

    :param N:
    :param verbose:
    :return:
    """
    found_candidate = True
    N_in_degree = N.in_degree
    N_out_degree = N.out_degree
    N_nodes = N.nodes
    N_predecessors = N.predecessors
    N_successors = N.successors
    while found_candidate:
        found_candidate = False
        for vertex in (v for v in N_nodes() if N_in_degree(v) == 1 and N_out_degree(v) == 1):
            if VERBOSE:
                print vertex, "was an indegree 1 outdegree 1 vertex"
            N.add_edge(N_predecessors(vertex)[0], N_successors(vertex)[0])
            N.remove_node(vertex)
            found_candidate = True


# Input: a rooted phylogenetic network N
# Output: number of leaves of N
def leafNumber(N):
    """

    :param N:
    :return:
    """
    return sum(1 for vertex, degree in N.out_degree_iter() if not degree)


# Input: a rooted phylogenetic network N
# Output: set of leaves of N
def leaves(N):
    """

    :param N:
    :return:
    """
    return {vertex for vertex, degree in N.out_degree_iter() if not degree}


def reticulations(N):
    """Returns the set of reticulations in N."""
    return {vertex for vertex, degree in N.in_degree_iter() if degree > 1}


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
        'folder', type=str,
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
    with open(os.path.join(folder, "results2.csv"),"w") as output:
        for data_file in glob.glob(os.path.join(folder, "data", file_filter)):
            # network initialization and preprocessing
            print "Treating file", data_file
            N = open_network(data_file)
            contract(N)

            # store identified classes and useful properties in a dictionary,
            # so we can use previously acquired knowledge to speed things up
            # (i.e., if we know that N is tree-child, then we know it's also
            # nearly tree-child so no need to ask a function for the result)
            N_classification = dict()

            # compute stable vertices
            if VERBOSE:
                print """== Computing the stable vertices =="""
            #t0=datetime.datetime.now()
            # knowing the root, the leaves and the reticulations is useful at
            # various stages, so let's compute all of that only once
            root_N = root(N)
            leaves_of_N = leaves(N)
            reticulations_of_N = reticulations(N)
            stableVertices = {
                v for v in N.nodes_iter() if not N.out_degree(v) or not N.in_degree(v)
                or isStable(v, N, root_N, leaves_of_N)
            }

            if VERBOSE:
                print "Stable vertices of N:", stableVertices
                print """===================================
            """
            #print "Time Stable vertices: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            line = data_file
            if VERBOSE:
                print "Edges of N: ", N.edges()

            # level
            #t0=datetime.datetime.now()
            level = computeLevel(N)
            if VERBOSE:
                print "Level of N:", level
                print """
              """
            N_classification["level"] = level
            line = ''.join((line, ";", str(level)))  # join is faster
            #print "Time Level: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # tree-child
            #t0=datetime.datetime.now()
            N_classification["tc"] = N_classification.get("level", -1) == 1 or\
                isTreeChild(N)
            if N_classification["tc"]:
                if VERBOSE:
                    print """== N is treechild.==
                 """
                line = ''.join((line, ";tc"))
            else:
                if VERBOSE:
                    print """== N is not treechild. ==
                 """
                line = ''.join((line, ";not tc"))
            #print "Time Tree-child: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # nearly tree-child
            #t0=datetime.datetime.now()
            ntc = N_classification.get("tc", False) or isNearlyTreeChild(
                N, stableVertices, reticulations_of_N
            )
            N_classification["ntc"] = ntc
            if VERBOSE:
                if ntc:
                    print """== N is nearly tree-child. ==
                 """
                else:
                    print """== N is not nearly tree-child. ==
                 """
            if ntc:
                line = ''.join((line, ";ntc"))
            else:
                line = ''.join((line, ";not ntc"))
            #print "Time Nearly tree-child: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # genetically stable
            #t0=datetime.datetime.now()
            genStab = N_classification.get("ntc", False) or \
                     isGeneticallyStable(N, stableVertices, reticulations_of_N)
            N_classification["gs"] = genStab
            if VERBOSE:
                if genStab:
                    print """== N is genetically stable. ==
                 """
                else:
                    print """== N is not genetically stable. ==
                 """
            if genStab:
                line = ''.join((line, ";gs"))
            else:
                line = ''.join((line, ";not gs"))
            #print "Time Genetically stable: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # tree-sibling
            #t0=datetime.datetime.now()

            trSib = N_classification.get("gs", False) or isTreeSibling(
                N, reticulations_of_N
            )
            N_classification["ts"] = trSib
            if VERBOSE:
                if trSib:
                    print """== N is tree-sibling. ==
                 """
                else:
                    print """== N is not tree-sibling. ==
                 """
            if trSib:
                line = ''.join((line, ";ts"))
            else:
                line = ''.join((line, ";not ts"))

            #print "Time Tree-sibling: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # reticulation-visible
            #t0=datetime.datetime.now()
            retVis = N_classification.get("gs", False) or isReticulationVisible(
                stableVertices, reticulations_of_N
            )
            N_classification["rv"] = retVis
            if VERBOSE:
                if retVis:
                    print """== N is reticulation visible. ==
                 """
                else:
                    print """== N is not reticulation visible. ==
                 """
            if retVis:
                line = ''.join((line, ";rv"))
            else:
                line = ''.join((line, ";not rv"))
            #print "Time Reticulation-visible: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # compressed
            #t0=datetime.datetime.now()
            comp = isCompressed(N, reticulations_of_N)
            N_classification["cp"] = comp
            if VERBOSE:
                if comp:
                    print """== N is compressed. ==
                 """
                else:
                    print """== N is not compressed. ==
                 """
            if comp:
                line = ''.join((line, ";cp"))
            else:
                line = ''.join((line, ";not cp"))
            #print "Time Compressed: "+str((datetime.datetime.now()-t0).microseconds)+"ms."


            # nearly-stable
            #t0=datetime.datetime.now()
            ns = N_classification.get("tc", False) or isNearlyStable(
                N, stableVertices
            )
            N_classification["ns"] = ns
            if VERBOSE:
                if ns:
                    print """== N is nearly stable. ==
                 """
                else:
                    print """== N is not nearly stable. ==
                 """
            if ns:
                line = ''.join((line, ";ns"))
            else:
                line = ''.join((line, ";not ns"))
            #print "Time Compressed: "+str((datetime.datetime.now()-t0).microseconds)+"ms."

            # write information about the network
            output.write(line + "\n")


if __name__ == '__main__':
    main()
