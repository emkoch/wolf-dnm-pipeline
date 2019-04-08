snakemake --cluster-config config/cluster.json \
	  --cluster "source activate wolf_dnm; sbatch -t {cluster.time} --mem {cluster.mem} -o {cluster.out} -e {cluster.err} --partition=jnovembre" \
	  --printshellcmds --jobs=10 --configfile config/config.json --verbose
