import os

pipeline_path = '/home/aaronfishman/repos/btb-seq/'
results_path = '/home/aaronfishman/validation-results/btb-seq/'

def main():
    # Housekeeping
    # make 
    if os.path.isdir(results_path):
        raise Exception("Output results path already exists")

    if not os.path.isdir(pipeline_path):
        raise Exception("Pipeline code repository not found")

    os.makedirs(results_path)

    # simuG - simulate the genome

    # read simulation -- chop up that genome

    # btb-seq

    # performance analysis

    return

if __name__ == '__main__':
    main()