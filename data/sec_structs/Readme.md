What are this RNA data files ?
===============================

## Raw (big) databases
* RNA-Strand 2.0 (secondary_structures_database.dbn) : this file is a dataset supposed to be identical to RNA-Strand 2.0 (actually the file is present on IBISC machines for years now and nobody remembers how it was built). The former RNA Strand website is not online anymore (http://rnasoft.ca/strand).
* bpRNA-1m_90 : this huge database gathers the data from other databases (CRW, PDB, Rfam, RNP, SPR, SRP, ...) and superseeds RNA-Strand (minus the structures that are only in NDB, sadly). Sequences have been prefiltered to have no more than 90% identity. Source : http://bprna.cgrb.oregonstate.edu/
* Pseudobase(++) : A database of biologically validated pseudoknots, from the time discovering a pseudoknot was something unusual. Pseudobase stays famous for its pseudoknot classification scheme. I scraped it myself to build the file. Source : https://www.ekevanbatenburg.nl/PKBASE/PKB.HTML 


## Filtered databases
* verified_secondary_structures.dbn : The subset of RNA-Strand that was experimentally validated (basically, the ones for which a 3D structure was available, so the ones from NDB and PDB).
* The _short.dbn ones : Same as its parent, but filtered using the filter.py script.
* pseudoknots.dbn : Audrey Legendre's scrap of Pseudobase, which, for an unknow reason, does not contain all the available data, but nice descriptions of what the RNAs are.


## Small test databases
* RNA-MoIP dataset : The cherry-picked cases presented in Reinhartz et al. 2012 to show RNA-MoIP's performance.
* applications.dbn : My cherry-picked cases presented in Becquey et al. 2020 to show Biorseo's performance.
* example.dbn : an example database with only one RNA, for testing purposes
* nothing.dbn : an example database with no RNAs, for testing purposes


Enjoy benchmarking RNA structure prediction tools.