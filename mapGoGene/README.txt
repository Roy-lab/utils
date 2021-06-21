Usage:
./mapGoGene <gene_association_file> <gene_ontology_terms.obo> biological_process <output_dir>/<prefix_for_output_files>

Example: (for Aspergillus fumigatus, in short, afum)
./mapGoGene afum_gene_association goslim_aspergillus_20151201.obo biological_process GO_test/afum

Restrictions:
The gene association file needs to have a species-specific phrase in its name. 
Currently, the recognized phrases are {fly, yeast, afum, mouse, human, pombe}.
If you wish to add a new species, please make the appropriate changes in "AnnotMgr.C".

