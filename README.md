# DeNovoLOEUF

Pre-requisites: seqr must be installed (uses hail) and will annotate the variants with VEP
Install seqr: https://github.com/broadinstitute/seqr/blob/master/deploy/LOCAL_INSTALL.md

Step 1: Once seqr is installed - run from command line in seqr directory to generate a hail matrix table of annotated variants

cd your/seqr/directory

docker-compose exec --detach pipeline-runner python3 -m seqr_loading SeqrMTToESTask --local-scheduler --reference-ht-path /seqr-reference-data/GRCh38/combined_reference_data_grch38.ht --clinvar-ht-path /seqr-reference-data/GRCh38/clinvar.GRCh38.2020-06-15.ht --vep-config-json-path /vep_configs/vep-GRCh38-loftee.json --es-host elasticsearch --es-index-min-num-shards 1 --sample-type WGS --es-index sopr --genome-version 38 --source-paths /input_vcfs/cohort.targ.vcf.gz --dest-path /input_vcfs/cohort.mt --dont-validate

Step 2: Now run DeNovoLOEUF.py

- Requires pedigree file for trio data (example provided in GitHub)
- list of GenCC genes (provided in GitHub)
- clinical data (example provided in GitHub)
