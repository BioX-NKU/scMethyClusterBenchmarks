SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

conda run -n Louvain python build.py
conda run -n Louvain python Leiden.py > log 2>&1

