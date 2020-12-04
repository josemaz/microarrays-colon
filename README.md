# Transcriptomic analyses of metformin treatment reveal novel mechanism of action in cervical cancer inhibition
This repository contains code and supplementary materials for paper named "Transcriptomic analyses of metformin treatment reveal novel mechanism of action in cervical cancer inhibition
". Karen Griselda De la Cruz-López,1,4 Heriberto A Valencia-González,4 Freddy Omar Beltrán-Anaya,2 Luis Alfaro-Ruiz,2 Alfredo Hidalgo-Miranda,2 Patricio Gariglio-Vidal,4, Jesús Espinal-Enríquez,6,7, Jose Maria Zamora Fuentes Alejandro García-Carrancá.5, *

## Data

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
