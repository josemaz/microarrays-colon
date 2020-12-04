# microarrays-colon
This repository contains code and supplementary materials for paper named "xxx". Karen Apellido, Jose Maria Zamora-Fuentes, Jesus Espinal-Enriquez

## Data

Put Data files (Data/Tx1-2s.CEL, Data/Vh1-2s.CEL, ...) in 
`$ ls Data/`
Rename explicitly files as:
`mv Tx2-4s.CEL Tx1-4s.CEL`
`mv Tx4-4s.CEL Tx2-4s.CEL`
`mv Tx6-4s.CEL Tx3-4s.CEL`
`mv Vh4-4s.CEL Vh3-4s.CEL`

## Diferential Expression

You can make this step with:

`Rscript R/deg-TxVh.R`
