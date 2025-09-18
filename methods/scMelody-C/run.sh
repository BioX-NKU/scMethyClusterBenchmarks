SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

n_clusters="6"

conda run -n scMelody Rscript data_generation.r > log 2>&1
conda run -n scMelody python cluster.py $n_clusters >> log 2>&1

