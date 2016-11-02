
metaPGAP : metagenomic Pan Genome Analysis Pipeline
===================================================

metaPGAP is a pipeline for building core genome using PAN genome approach. This pipeline takes complete or draft genome assemblies and perform annotation using Prokka. The Prokka predictions then used to build core genome using get_homologues tool. For core genes phylogeny Mafft and RaXML used.

---------------------------------------------------
Full list of required softwares and dependencies:<br/>

metaPGAP steps:<br/>
===============
(1). Download data from https://github.com/mitul-patel/data/archive/master.zip<br/>
(2). Genome Annotation using Prokka<br/>
(3). PAN genome analysis using get_homologues: BDBH, COG and OMCL<br/>
(4). Multiple sequence alignment of CORE genes<br/>
(5). Phylogenetic analysis of CORE genes using RaXML<br/>
(6). Visualization of phylogenetic tree using Newick tools<br/>

---------------------------------------------------
metaPGAP requirements:
======================
Python (http://python.org)<br/>
BioPython (1.5 or higher) with NumPy (http://biopython.org/wiki/Download)<br/>
Prokka (https://github.com/tseemann/prokka)<br/>
get_homologues (https://github.com/eead-csic-compbio/get_homologues)<br/>
Mafft (http://mafft.cbrc.jp/alignment/software/)<br/>
Perl (http://Perl.org)<br/>
AMAS (https://github.com/marekborowiec/AMAS)<br/>
RaXML (https://github.com/stamatak/standard-RAxML)<br/>
Newick tools (http://cegg.unige.ch/newick_utils)<br/>

---------------------------------------------------

Install metaPGAP
================

wget 
unzip 


---------------------------------------------------

Running metaPGAP
================

cd 
python metaPGAP.py

---------------------------------------------------
