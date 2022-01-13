conda activate cpdb

awk 'NR==1{print "Gene\t",$0}' RNAdata.txt > RNAdata_temp.txt 
awk 'NR!=1' RNAdata.txt >> RNAdata_temp.txt 
mv RNAdata_temp.txt RNAdata.txt 

#Additionally, after the "Gene" description in the first line of RNAdata is a space that leads to an error message and needs to be manually deleted

#INDIR="~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/cellphonedb" directory to run this code

cellphonedb method statistical_analysis --counts-data hgnc_symbol celltype/celltype_annotation.txt RNAdata.txt --project-name celltype --output-path celltype
cellphonedb method statistical_analysis --counts-data hgnc_symbol tumorclusters/tumorsubclusters_annotation.txt RNAdata.txt --project-name tumorclusters --output-path tumorclusters
cellphonedb method statistical_analysis --counts-data hgnc_symbol allclusters/celltypeclusters_annotation.txt RNAdata.txt --project-name allclusters --output-path allclusters
