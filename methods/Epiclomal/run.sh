SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

rm -r ./config/runs_epiclomal/*
rm -r ./config/epiclomal_process/*

bash ./config/process_real_data/run.sh
bash ./config/real_data/run.sh

conda run -n Epiclomal python Performance_Evaluate.py > log 2>&1
