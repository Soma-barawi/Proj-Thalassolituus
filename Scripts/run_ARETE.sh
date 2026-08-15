# Run the ARETE workflow on all Thalassolituus genomes 

#!/usr/bin/bash
nextflow run Application-ARETE/ --input_sample_table '/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/samplesheet.csv' \
	--db_cache '/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Application-ARETE/dbcache' \
	--bakta_db '/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Application-ARETE/dbcache/baktadb/db-light/' \
	--annotation_tools 'mobsuite,rgi,cazy,vfdb,iceberg,bacmet,islandpath,phispy,integronfinder,report' \
	--outdir all_thalasso_output \
  	--run_recombination \
	--skip_poppunk \
	--run_evolccm \
	--run_rspr \
	-entry annotation \
	-profile 'singularity,narval' \
	-resume

