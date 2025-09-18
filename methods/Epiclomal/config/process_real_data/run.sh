SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

conda run -n Epiclomal snakemake --forceall -s ../../snakemake/process_real_data/Snakefile --configfile config.yaml --cores 8 > log 2>&1
conda run -n Epiclomal snakemake -s ../../snakemake/process_real_data/Snakefile --unlock --configfile config.yaml

