import utils
import argparse
import glob
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot(root_path):
    utils.assert_path_exists(root_path)

    dirs = glob.glob(root_path + "/*/stats.csv")

    values = pd.DataFrame(index =[range(1, len(dirs)+1)])
    branches=['branch:']

    for directory in dirs:
        path = Path(directory)
        branch_name = path.parts[-2]
        branches.append(branch_name)
        csv = pd.read_csv(directory)
        values[branch_name] = csv['f_score']

    plt.boxplot(values, showfliers = True)
    plt.xticks(range(0, len(branches)),  branches, rotation=-20)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make box-whisker plot for git branches")

    parser.add_argument("root_path", help="Path to directory containing a directory for each \
        branch where results are saved")

    args = parser.parse_args(sys.argv[1:])

    plot(args.root_path)