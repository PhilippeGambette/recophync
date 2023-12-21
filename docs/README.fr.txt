Ce fichier README a été généré le 2023-06-13 par Philippe Gambette.

Dernière mise-à-jour le : 2023-06-13.

# INFORMATIONS GÉNÉRALES

## Titre du jeu de données : Phylogenetic networks found in scientific publications

## DOI: 10.57745/VIW7B2
 
## Adresse de contact : philippe.gambette@univ-mlv.fr

# INFORMATIONS MÉTHODOLOGIQUES

## Description des sources et méthodes utilisées pour collecter et générer les données :

Collecte manuelle de publications contenant des figures présentant des réseaux phylogénétiques, réunies dans le cadre d'une veille scientifique relative aux réseaux phylogénétiques.

## Méthodes de traitement des données :

Encodage manuel du réseau comme une liste d'arètes :
* en général par un parcours en profondeur
* en reprenant les étiquettes des feuilles 
* en reprenant les étiquettes de nœuds internes s'ils sont étiquetés, en utilisant des étiquettes ad hoc (nombres, lettres, ou lettre suivie d'un nombre)

## Procédures d’assurance-qualité appliquées sur les données :

Utilisation de l'outil de visualisation de réseaux phylogénétiques à l'adresse https://phylnet.univ-mlv.fr/recophync/networkDraw.php pour détecter visuellement des erreurs éventuelles (oubli d'arête, erreur d'étiquette) qui conduisent généralement à une déconnexion du réseau.

## Autres informations contextuelles :

Pour visualiser ou manipuler ces réseaux dans le logiciel Dendroscope 3, on peut convertir chaque fichier d'extension .el vers le format eNewick de Dendroscope 3, en le visualisant à l'aide de l'outil disponible à l'adresse https://phylnet.univ-mlv.fr/recophync/networkDraw.php qui fera aussi la conversion à ce format.

# APERCU DES DONNEES ET FICHIERS

## Convention de nommage des fichiers :

n[i].el
où [i] représente le numéro du réseau, compris entre 1 et le nombre total de réseaux.

## Arborescence/plan de classement des fichiers :

À la racine :
* les fichiers correspondant aux réseaux, d'extension .el (edge list)
* un fichier .tab contenant les métadonnées relatives aux réseaux :
** nb : numéro du réseau 
** id	: identifiant du réseau
** HTML caption : code HTML d'une référence de la source du réseau
** caption : référence de la source du réseau au format texte