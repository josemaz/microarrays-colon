# Transcriptomic analyses of metformin treatment reveal novel mechanism of action in cervical cancer inhibition
This repository contains code and supplementary materials for paper named "Transcriptomic analyses of metformin treatment reveal novel mechanism of action in cervical cancer inhibition". Karen Griselda De la Cruz-López, Heriberto A Valencia-González, Freddy Omar Beltrán-Anaya, Luis Alfaro-Ruiz, Alfredo Hidalgo-Miranda, Patricio Gariglio-Vidal, Jesús Espinal-Enríquez, Jose Maria Zamora Fuentes, Alejandro García-Carrancá.

# Tree

o
|-- Data
|   `-- .CEL files
|-- Docs
|   `-- Documentation and Notes
`-- R
    `-- Source Code


## Data

*Data availibality is pending approval.*

Put Data files (Data/Tx1-2s.CEL, Data/Vh1-2s.CEL, ...) in :

`$ mkdir Data/`

`$ ls Data/`

Rename explicitly files as:

`mv Tx2-4s.CEL Tx1-4s.CEL`

`mv Tx4-4s.CEL Tx2-4s.CEL`

`mv Tx6-4s.CEL Tx3-4s.CEL`

`mv Vh4-4s.CEL Vh3-4s.CEL`

## Differential Expression

You can make this step with:

`Rscript R/deg-TxVh.R`

#TCGA

TCGA-HM-A3JJ, 	TCGA-FU-A3EO, TCGA-MY-A5BF