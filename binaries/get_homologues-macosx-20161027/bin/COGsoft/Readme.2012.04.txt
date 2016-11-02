############################################################
0. General remarks.
############################################################

This is the April 2012 release of the Eugene Koonin's group COG software, featuring the EdgeSearch algorithm.

#-----------------------------------------------------------
0.1. Disclaimers.

This version of COG-related software is a part of a research project of Eugene Koonin's group at the National Center for Biotechnology Information (NCBI), National Library of Medicine (NLM) at the National Institutes of Health (NIH). This is not an official NCBI software product. The software is distributed "as is" with no stated or implied warranty whatsoever. Although we'll make a reasonable effort to help if problems arise, we are in no way committed to support and maintenance of this software.

#-----------------------------------------------------------
0.2. Credits, contacts etc.

COG software:

Credits:

David M. Kristensen, NCBI
Alexander V. Sorokin, NCBI
Pavel S. Novichkov, NCBI
Yuri I. Wolf, NCBI
Eugene V. Koonin, NCBI

Contact:

Yuri I. Wolf <wolf@ncbi.nlm.nih.gov>

References:

Kristensen DM, Kannan L, Coleman MK, Wolf YI, Sorokin A, Koonin EV, Mushegian A
A low-polynomial algorithm for assembling clusters of orthologous groups from intergenomic symmetric best matches.
Bioinformatics. 2010
http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq229?ijkey=zD7TIWnncGvDGYE&keytype=ref

Other papers that heavily rely on COGs, constructed with this software:

Makarova KS, Sorokin AV, Novichkov PS, Wolf YI, Koonin EV
Clusters of orthologous genes for 41 archaeal genomes and implications for evolutionary genomics of archaea.
Biol Direct. 2, 33, 2007

A blanket citation for the COG approach:

Tatusov RL, Koonin EV, Lipman DJ
A genomic perspective on protein families.
Science 278, 631-637, 1997


EdgeSearch algorithm/implementation:

Credits:

David Kristensen
Lavanya Kannan
Michael Coleman
Arcady Mushegian
all at Stowers Institute for Medical Research, Bioinformatics Department

Contact:

David M. Kristensen <David.Kristensen@nih.gov>


############################################################
1. Software.
############################################################

This distribution package contains the UNIX source files for the following programs (each in its own directory):

COGmakehash
COGreadblast
COGlse
COGtriangles
COGcognitor

Run UNIX make command in each directory to create executable binaries. You might need to modify the makefiles if you are using a C++ compiler other than GNU g++.

All programs will print out a short description of the options when executed with "-h" command line flag.

A note on Mac OSX: the hidden file .DS_Store.tab in each directory interferes with the execution of the COGreadblast program, since it reads all files with the extension '.tab'. One possible work-around is to change the file mask to something other than '.tab' (FMASK variable in bc.h).

############################################################
2. Data.
############################################################

#-----------------------------------------------------------
2.0. Remarks on data.

The description of the data set requirements that will follow assumes that standard NCBI BLAST suite software will be used as a part of the data processing workflow. Thus, the description reflects the limitations imposed by the current (as of May 2010) versions of that software. In many cases the COG-related software has additional flexibility (e.g. less strict conventions on naming proteins and domains) which can be used in alternative workflows. Here we will mostly ignore these options and present the more strict (i.e. most compatible) set of requirements. The internal data formats are described below to allow the possibility to supply these data to the COG software outside of the normal processing workflow.

#-----------------------------------------------------------
2.1. Protein sets.

Each protein is expected to have a unique short identifier - either a string of alphanumeric characters or an integer decimal number not exceeding 2147483647 (NB - longer strings of decimal digits can be used as sequence IDs if they are presented in the context of text IDs rather than numbers). For compatibility with BLAST it is strongly advised that the sequence IDs in FASTA files are constructed as "lcl|<text-id>" or "gi|<num-id>". Use of some non-alphanumeric characters (dot and dash) in names is allowed but discouraged; others (comma, semicolon, colon etc.) are prohibited. Some software does not distinguish upper and lower case letters, so it is strongly recommended to make the names unique in the case-insensitive mode.

#-----------------------------------------------------------
2.2. Proteins and genomes.

Each protein is expected to be assigned to one and only one genome; lists of proteins assigned to a genome are expected to be complete (i.e. include the full set of proteins, encoded by a genome). The protein-genome data are given line-by-line in a comma-separated file as "<prot-id>,<genome-id>". Genome IDs are expected to be alphanumeric strings.

#-----------------------------------------------------------
2.3. BLAST results.

Typically two BLAST search passes are required - one without any low-complexity filtering and composition-based statistics (as it produces scores that, apparently, better reflect the phylogenetic distances); another with low-complexity filtering and/or composition-based statistics (as it produces less false positives). By convention, the COG software expects BLAST results in tabular ("-m 8" or "-m 9") format; each set of BLAST results is expected to reside in it's own directory in *.tab files. Depending on the precise composition of the query and subject deflines, BLAST output may contain either bare protein IDs (e.g. "18978112" or "MTH068") or full sequence IDs from the original fasta files (e.g. "gi|18978112|ref|NP_579469.1|" or "lcl|MTH068"). The BLAST postprocessing software will split the string by delimiter characters and take the specified (by default the 2nd) token for an ID. You need to make sure that BLAST statistics is consistent and comparable between different runs; one useful hint includes forcing the same effective database size in all searches (using "-z" BLAST parameter) when running BLAST searches in batches against different databases.

#-----------------------------------------------------------
2.4. Internal IDs.

Internally all COG software uses numerical IDs for the sequences. Correspondence between the user-supplied IDs and these internal IDs are established by the COGmakehash program (see p. 3.1.1). The correspondence file is called hash.csv and resides in a separate directory together with other processed BLAST data. The file consists of "<num-prot-id>,<user-prot-id>" records.

#-----------------------------------------------------------
2.5. Self-similarity data.

All BLAST similarity scores are measured against the self-similarity of the proteins involved. The self-similarity data is stored in the self.csv file residing in the processed BLAST data directory. The file consists of "<num-prot-id>,<prot-length>,<self-score>" records. This file is prepared by the COGreadblast program (see p. 3.1.2).

#-----------------------------------------------------------
2.6. BLAST hits data.

Processed BLAST hits from the unfiltered BLAST search are stored in the hits.csv file residing in the processed BLAST data directory. The file consists of "<query-num-prot-id>,<subject-num-prot-id>,<query-start>,<query-end>,<subject-start>,<subject-end>,<e-value>,<score>" records. All hits for the same query are stored in contiguous blocks sorted by decreasing score (increasing e-value). This file is prepared by the COGreadblast program (see p. 3.1.2).

#-----------------------------------------------------------
2.7. BLAST neighbors.

Query-subject pairs from the filtered BLAST search are stored in the query2subject.csv file residing in the processed BLAST data directory. The file consists of "<query-num-prot-id>,<subject-num-prot-id>" records. This file is prepared by the COGreadblast program (see p. 3.1.2).

#-----------------------------------------------------------
2.8. LSE job description file.

This file is used by the COGlse program and consists of "<query-genome-id>,<outgroup-genome-id>" records (see p. 3.1.3 for details and examples).

#-----------------------------------------------------------
2.9. LSE data.

The information of lineage-specific (intra-genome) expansions is stored in the form of "<prot-id-1>,<prot-id-2>..." records. Proteins that do not belong to LSEs are listed as single entries. This file is prepared by the COGlse program (see p. 3.1.3).

#-----------------------------------------------------------
2.10. Clustering data.

The output of the COGtriangles program (see p. 3.1.4) is in the following format: "<prot-id>,<genome-id>,<source-prot-id>,<source-prot-length>,<source-prot-start>,<source-prot-end>,<cluster-id>,...". In the basic usage "<prot-id>" and "<source-prot-id>" refer to the same protein record; thus the "<source-prot-length>" refers to the original protein length, "<source-prot-start>" is always 1 and "<source-prot-end>" is always equal to the "<source-prot-length>". In advanced usage the original proteins could be split into domains clustered differently. In this case the "<prot-id>" would refer to the unique identifier assigned to the domain while the rest of the info would refer to the original protein, its length and start-end coordinates relative to it.

As of May 2010, COGtriangles now runs using the EdgeSearch algorithm, which allows proteins to belong to multiple COGs, rather than arbitrarily choosing one when multiple choices could be made, and also produces 2 additional output files, with the hard-coded names all-edges.txt and cog-edges.txt that represent the graph upon which the COGs were built (all of its edges and the edges present in each COG, respectively). A perl script is provided to help deal with this change, for example by reporting all the COGs that a single protein appears to belong to on a single line (with the special prefix "multi:"), or in another mode, to mimic the behavior of the original COGtriangles program for backwards compatability by removing multiple instances of a protein, allowing each protein to only belong to the largest COG that it is seen to belong to (i.e., the one most likely to have been encountered first).

#-----------------------------------------------------------
2.11. COGNITOR data.

The output of the COGcognitor program (see p. 3.2.3) is in the following format: "<prot-id>,<prot-length>,<prot-start>,<prot-end>,<cognitor-score>,<cluster-id>,...". The "<cognitor-score>" indicates the relative confidence in the assignment of the query protein fragment to the cluster.

############################################################
3. Procedure.
############################################################

#-----------------------------------------------------------
3.1. Making COGs.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.1.0. Running BLAST.

Technically, running BLAST searches is outside the scope of COG software. The following sequence of commands is given exclusively for the sake of example. Suppose that you have proteins sets from three genomes in three separate FASTA files: Gen1.fa, Gen2.fa and Gen3.fa.

$ cat Gen1.fa Gen2.fa Gen3.fa > GenThree.fa # pools the sets together

$ makeblastdb -in GenThree.fa -dbtype prot -out GenThree # formats a BLASTable database

$ psiblast -query GenThree.fa -db GenThree -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTno/ThreeByThree.tab # unfiltered BLAST results in the ./BLASTno/ directory

$ psiblast -query GenThree.fa -db GenThree -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out BLASTff/ThreeByThree.tab # filtered BLAST results in the ./BLASTff/ directory

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.1.1. Preparation of the "sequence Universe".

Internally all COG software uses numerical IDs for the sequences. The first step in the data preparation, thus, involves making a table that connects these internal IDs and the IDs in the user-supplied data. Suppose you have a file GenThree.p2o.csv that lists all sequences involved in the project (format "<prot-id>,<genome-id>"). The following command will be used:

$ COGmakehash -i=GenThree.p2o.csv -o=./BLASTconv -s="," -n=1 # makes ./BLASTconv/hash.csv file 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.1.2. Processing of BLAST results.

The program will read the BLAST results from the ./BLASTno/ and ./BLASTff/ directories and will store the pre-processed results in the ./BLASTconv/ directory.

$ COGreadblast -d=./BLASTconv -u=./BLASTno -f=./BLASTff -s=./BLASTno -e=0.1 -q=2 -t=2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.1.3. Lineage specific expansions.

To facilitate the clustering, lineage-specific expansions (a.k.a. in-paralogs) need to be collected first. The program compares all intra-genome similarities to the corresponding outgroup-genome similarities and cluster together proteins from the same genome that are more similar to each other than to any outgroup protein. The comparison is guided by a "job-description" file which consists of "<query-genome-id>,<outgroup-genome-id>" lines. In the simplest case (e.g. three genomes that are roughly equidistant from each other) this file will list each of the query genomes with all other genomes as outgroups:

Gen1,Gen2
Gen1,Gen3
Gen2,Gen1
Gen2,Gen3
Gen3,Gen1
Gen3,Gen2

If, in the above example, Gen1 and Gen2 are closely related strains from one species while Gen3 is relatively distant, you might want to omit those two from the outgroup lists for each other, reducing the job description file to:

Gen1,Gen3
Gen2,Gen3
Gen3,Gen1
Gen3,Gen2

The command would be run as:

$ COGlse -d=./BLASTconv -j=GenThree.job.csv -p=GenThree.p2o.csv -o=GenThree.lse.csv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.1.4. Making clusters from symmetrical best hits.

First the symmetrical best genome-to-genome hits satisfying the reciprocal sequence length coverage criteria (0.75 of both sequences are covered by an HSP by default) are identified. A set of triangles is built from these data; triangles sharing a side are clustered together. After all sequences clustered by triangles are removed, symmetrical best hit pairs are collected from those that are left.

$ COGtriangles -i=./BLASTconv -q=GenThree.p2o.csv -l=GenThree.lse.csv -o=GenThree.cls.csv -t=0.5 -e=0.01 -n="CLS" -s=1

Note that the new COGtriangles program (see also p. 2.10) allows proteins to belong to multiple COGs if they should appear to do so, which can serve as a useful diagnostic indicator of multi-domain proteins or unresolved paralogy that should be more closely investigated. To make the file easier to read, the optional perl script COGtriangles.reformat.pl will change the output to report all the COGs that a single protein belongs to on a single line (with the special prefix "multi:" for such cases). This is the recommended option. The wrapper assumes your perl is in /usr/bin/perl, and that both COGtriangles.reformat.pl and the COGtriangles binary are in the same directory. This script can be run either as a wrapper, where it will run COGtriangles, save the original results (your output file with ".orig" appended to the name), and replace the output file with the reformatted COGs:

$ COGtriangles.reformat.pl -r -i=./BLASTconv -q=GenThree.p2o.csv -l=GenThree.lse.csv -o=GenThree.cls.csv -t=0.5 -e=0.01 -n="CLS" -s=1

If desired, this script could also/instead be run on existing COG results with the alternate syntax:

$ COGtriangles.reformat.pl originalCOGsfile outputCOGsfile

This script can work in several modes, which are explained in the usage statement printed when running the script with no parameters, and also provides additional options such as filtering by size.

#-----------------------------------------------------------
3.2. COGnitor.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.2.0. Running BLAST.

Technically, running BLAST searches is outside the scope of COG software. The following sequence of commands is given exclusively for the sake of example. Suppose that you have proteins sets from two query genomes in two separate FASTA files: Gen1.fa and Gen2.fa and the orthology domains in another FASTA file: COGs.fa.

$ cat Gen1.fa Gen2.fa > GenQuery.fa # pools the query sets together

$ makeblastdb -in GenThree.fa -dbtype prot -out GenThree # formats a BLASTable database for the query set

$ makeblastdb -in COGs.fa -dbtype prot -out COGs # formats a BLASTable database for the target set

$ psiblast -query GenThree.fa -db GenThree -show_gis -outfmt 7 -num_descriptions 10 -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTss/QuerySelf.tab # unfiltered self-hit BLAST results in the ./BLASTss/ directory

$ psiblast -query GenThree.fa -db COGs -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTno/QueryCOGs.tab # unfiltered BLAST results in the ./BLASTno/ directory

$ psiblast -query GenThree.fa -db COGs -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out BLASTff/QueryCOGs.tab # filtered BLAST results in the ./BLASTff/ directory

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.2.1. Preparation of the "sequence Universe".

Internally all COG software uses numerical IDs for the sequences. The first step in the data preparation, thus, involves making a table that connects these internal IDs and the IDs in the user-supplied data. Suppose you need to have a file GenQuery.p2o.csv that lists all sequences involved in the query genomes (format "<prot-id>,<genome-id>") and COG.p2o.csv that lists all orthology domains. The following commands will then be used:

$ cat GenQuery.p2o.csv COG.p2o.csv > tmp.p2o.csv # pools the lists together

$ COGmakehash -i=tmp.p2o.csv -o=./BLASTcogn -s="," -n=1 # makes ./BLASTcogn/hash.csv file 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.2.2. Processing of BLAST results.

The program will read the BLAST results from the ./BLASTno/ and ./BLASTff/ directories and will store the pre-processed results in the ./BLASTconv/ directory.

$ COGreadblast -d=./BLASTcogn -u=./BLASTno -f=./BLASTff -s=./BLASTss -e=0.1 -q=2 -t=2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3.2.3. Running COGNITOR.

To run COGNITOR you need a COG domain assignment file (as described in 2.10.). If your file is called COGs.csv, the following command will be used:

$ COGcognitor -i=./BLASTcogn -t=COGs.csv -q=GenQuery.p2o.csv -o=GenQuery.COG.csv # COGNITOR results in GenQuery.COG.csv

############################################################
4. Release history
############################################################

A short info on what has changed.

#-----------------------------------------------------------
4.2.2. April 2012

Bug fix involving best hit pairs with edges that are not symmetrical. Now every protein in the input files should also appear in the raw output cluster file. Thanks to Natalya Yutin at NCBI.

#-----------------------------------------------------------
4.2.2. March 2012

Combined the two perl scripts introduced in the previous version (COGtriangles.back.pl and COGtriangles.combine.pl) into a single script (COGtriangles.reformat.pl), enhanced further with additional modes of action and filtering options.

#-----------------------------------------------------------
4.2.1. November 2011

Slight changes to accomodate changes in newer C++ compilers, and added two perl scripts to process COGtriangles results (both help introduce better backwards compatability - one by discarding information for strict compatability while the other introduces slight changes in the format to place all COGs that a protein should appear to belong to on a single line).

#-----------------------------------------------------------
4.2.0. May 2010

COGtriangles replaced with a version that uses the EdgeSearch algorithm.

#-----------------------------------------------------------
4.1.2. July 2008

A bug in COGcognitor fixed - thanks to Naama Wald of Hebrew University

#-----------------------------------------------------------
4.1.1. June 2008

A minor bug in COGreadblast fixed - thanks to David Kristensen of Stowers Institute.

#-----------------------------------------------------------
4.1.0. March 2008

COGNITOR added.

#-----------------------------------------------------------
4.0. May 2007

The first public version.
