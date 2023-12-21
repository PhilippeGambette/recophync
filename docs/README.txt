This README file was generated on 2023-06-13 by Philippe Gambette.
Last updated: 2023-06-13.
 
# GENERAL INFORMATION

## Dataset title: Phylogenetic networks found in scientific publications

## DOI: 10.57745/VIW7B2
 
## Contact email: philippe.gambette@univ-mlv.fr
 
# METHODOLOGICAL INFORMATION 
 
## Methods for processing the data: 

Manual collection of publications containing figures presenting phylogenetic networks, gathered as part of a scientific watch relating to phylogenetic networks.

## Data processing methods:

Manual encoding of the network as a list of aretes:
* in general, by means of an in-depth path
* using leaf labels 
* by taking over the labels of internal nodes if they are labeled, using ad hoc labels (numbers, letters, or letter followed by a number).

## Quality-assurance procedures performed on the data: 

Use of the phylogenetic network visualization tool at https://phylnet.univ-mlv.fr/recophync/networkDraw.php to visually detect possible errors (missing edge, label error), which generally lead to network disconnection.
 
## Other contextual information:

To view or manipulate these networks in Dendroscope 3 software, you can convert each .el extension file to Dendroscope 3's eNewick format, by viewing it using the tool available at https://phylnet.univ-mlv.fr/recophync/networkDraw.php, which will also convert to this format.
 
 # DATA & FILE OVERVIEW

## File naming convention :

n[i].el
where [i] represents the network number, between 1 and the total number of networks.

## File hierarchy convention:

At the root :
* files corresponding to networks, with .el extension (edge list)
* a .tab file containing network metadata:
** nb: network number 
** id: network identifier
** HTML caption: HTML code of a network source reference
** caption: network source reference in text format 