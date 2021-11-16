import glob

import pandas as pd
import os

root = '/home/aaronfishman/mnt/fsx-027/ofat-25/'
branch_paths = sorted(glob.glob(root + './*'))

results = []

for branch_path in branch_paths:
    branch_name = os.path.basename(branch_path)

    # Load csv
    df = pd.read_csv(branch_path+"/stats.csv")

    result = {"name": branch_name}
    for i, row in df.iterrows():
        result[row["name"]] = row["total_errors"]

    results.append(result)

combined = pd.DataFrame(results)

a = 1