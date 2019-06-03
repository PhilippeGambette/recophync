#!/usr/bin/env python

#Python3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import os
import sys
import glob

from properties.common import numerical_properties, network_types
from Phylo import PhyloNetwork


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
          print("property " + short + ": " + PN.properties[short].report())
          line += ';' + PN.properties[short].report()

        # write information about the network
        output.write(line + "\n")








if __name__ == '__main__':
    main()
