# recophync
[![Build Status](https://travis-ci.com/AndrewQuijano/recophync.svg?branch=master)](https://travis-ci.com/AndrewQuijano/recophync)  
  
[![codecov](https://codecov.io/gh/AndrewQuijano/recophync/branch/main/graph/badge.svg?token=QdrNMDBjw0)](https://codecov.io/gh/AndrewQuijano/recophync)

[RecoPhyNC](http://phylnet.univ-mlv.fr/recophync/) is a Python software which tests 
whether a phylogenetic network, given as a file containing a list of arcs is:
* &#9744; level-k
* &#9745; tree-child
* &#9745; nearly tree-child
* &#9745; genetically stable
* &#9745; reticulation-visible
* &#9745; tree-sibling
* &#9745; compressed
* &#9744; nearly stable 

It has been updated to work with `Python3` and the latest version of `networkx`

Warning: there may be a bug in the currently implemented nearly stable detection algorithm (see https://github.com/PhilippeGambette/recophync/pull/1) as well as in the level computation (see https://github.com/PhilippeGambette/recophync/commit/190651a3943efb479ecaf9d1712969fd3b5fc686).
