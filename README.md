# T-shaped alignments integrating HIV-1 near full-length genome and partial pol sequences can improve phylodynamic inference of transmission clusters

This is the repository for our work on *T-shaped alignments*, a method to integrate near HIV-1 full-length genome and partial pol sequences for the purposes of phylodynamic inference of transmission clusters. It is under review at PLoS Computational Biology.

## Data

Source data can be downloaded from the [HIV Sequence DB](https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html). Select 2018, HIV-1, full genome, alignment.

First-stage intermediate data, i.e. masked alignments, phylogenetic trees, are available by request.

Intermediate rds data for reproducing figures is directly in this repository.

## To Reproduce

To reproduce first-stage intermediate data directly from the source data, follow the scripts in order in the `scripts` directory. Scripts were written and run on the Slurm scheduler on [OSCAR](https://ccv.brown.edu/services/computing/), and you will likely need access to a HPC cluster yourself to run so many phylogenetic trees.

To reproduce figures, change the `PATH` in `notebooks/figures.qmd` to where this repository is located, and execute. You will need the appropriate R packages as well as Quarto installed, otherwise it should run on your local machine without issue.

## Background

### Genome

![HIV Genome Diagram](img/bi67_0001_1.jpeg)

> The HIV-1 genome encodes nine open reading frames. Three of these encode the Gag, Pol, and Env polyproteins, which are subsequently proteolyzed into individual proteins common to all retroviruses. The four Gag proteins, MA (matrix), CA (capsid), NC (nucleocapsid), and p6, and the two Env proteins, SU (surface or gp120) and TM (transmembrane or gp41), are structural components that make up the core of the virion and outer membrane envelope. The three Pol proteins, PR (protease), RT (reverse transcriptase), and IN (integrase), provide essential enzymatic functions and are also encapsulated within the particle. HIV-1 encodes six additional proteins, often called accessory proteins, three of which (Vif, Vpr, and Nef) are found in the viral particle. Two other accessory proteins, Tat and Rev, provide essential gene regulatory functions, and the last protein, Vpu, indirectly assists in assembly of the virion. The retroviral genome is encoded by an ∼9-kb RNA, and two genomic-length RNA molecules are also packaged in the particle. Thus, in simplistic terms, HIV-1 may be considered as a molecular entity consisting of 15 proteins and one RNA. [https://www.annualreviews.org/doi/10.1146/annurev.biochem.67.1.1](https://www.annualreviews.org/doi/10.1146/annurev.biochem.67.1.1)
