# calls2xls
A set of Perl scripts for intended to simplify variant calling and enable interpretation of NGS data in the clinical genetics setting. The scripts are, in their current state, taylored to be compatible with the the bcbio-nextgen pipeline and the Alamut Visual variant interpretation software.

The scripts enable:

1. Simple generation of bcbio-nextgen config files and execution of the pipeline.

2. Building and populating a local database of variant and coverage depth information.

3. Generation of Excel 2003 files with variant, quality, and other information based on pre-defined subsets of the data (genes) with respect to a provided indication. Links are provided in the excel files that enable simple interaction with Alamut Visual. Variants become annotated with SnpEff as well as local variant frequencies and classifications from a local Alamut Visual variant repository.

4. Syncronizing Excel, Bam, and other quality-related files to a share accessible by the end user (MD / clinical scientist).





