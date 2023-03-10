##########################################################################################
#																						 #
#																						 #	
#																						 #
#			  	Workflow for analysis of pairwise Fst between localities				 #
#		 		 	  	  in the phenotypic integration study							 #
#																						 #
#																						 #
#																						 #
##########################################################################################

Created 02.23.2023 by Drew Schield

This README contains details processing and analysis steps to summarize population genetic
structure among barn swallow localities used in our phenotypic integration study. The SNP
data will be pulled from the 2021 Mol Ecol global whole genome dataset, including all
localities where there are two or more individuals.

The plan will be to perform pairwise genome-wide Fst analyses between each locality and
the savignii locality in Egypt, so that all comparisons are standardized by a single pair
locality.

General overview of sections in this README:
1. Organizing SNP data input
2. Pairwise Fst analyses


##########################################################################################
#																						 #
#																						 #
#			  				  1. Organizing SNP data input								 #
#		 		 	  	  																 #
#																						 #
##########################################################################################

We have filtered SNPs from the 2021 whole genome dataset for the majority of localities
included in the phenotypic integration dataset. I've since backed these up to VorpalBucket
to save space on Terminator, but will retrieve a reasonable filtered VCF for analysis.

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd /data3/
$mkdir hirundo_phenotypic_integration
$cd hirundo_phenotypic_integration
$mkdir vcf
$mkdir fst

------------------------------------------------------------------------------------------
Retrieved VCF:
------------------------------------------------------------------------------------------

[On desktop]:

$cd /Volumes/VorpalBucket/hirundo_MolEcol_backup_06.23.21/vcf
$scp hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin10kb.ingroup.auto.vcf drewschield@terminator:/data3/hirundo_phenotypic_integration/vcf

------------------------------------------------------------------------------------------
Formatted population lists:
------------------------------------------------------------------------------------------

[On Terminator]:

All formatted as `popmap.<locality>`. There are 20 of them.



##########################################################################################
#																						 #
#																						 #
#			  				  	2. Pairwise Fst analysis								 #
#		 		 	  	  																 #
#																						 #
##########################################################################################

Now we'll perform pairwise Fst analyses.

------------------------------------------------------------------------------------------
Wrote script to run all pairwise against the savignii 'Egypt' population:
------------------------------------------------------------------------------------------

runFst.sh
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for pop in popmap.*; do
	vcftools --vcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin10kb.ingroup.auto.vcf --weir-fst-pop $pop --weir-fst-pop popmap.egypt --out ./fst/$pop
done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

sh runFst.sh

------------------------------------------------------------------------------------------
Look at weighted Fst values:
------------------------------------------------------------------------------------------

$cd fst
$for i in *.log; do echo $i; grep 'weighted Fst estimate' $i; done

Weighted Fst values:

alzamay	0.078706
boatu	0.096982
changchun	0.084985
colorado	0.080094
egypt	-0.12942
harbin	0.10087
israel	0.030276
jiuquan	0.055283
khingui	0.060532
marrakesh	0.026705
moscow	0.031616
narin-talacha	0.071967
qinhuangdao	0.10307
qiqihar	0.082045
takecha	0.079836
tokyo	0.10487
yekatarinburg	0.036518
zakaltoose	0.069245
zhangye	0.039776
zhengzhou	0.095367























