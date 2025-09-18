SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

n_clusters="6"

conda run -n Louvain python build.py
conda run -n Louvain python Leiden-C.py $n_clusters > log 2>&1

