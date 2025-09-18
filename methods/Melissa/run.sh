SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

conda run -n Melissa Rscript run_melissa.R > log 2>&1
conda run -n Melissa python Performance_Evaluate.py >> log 2>&1
