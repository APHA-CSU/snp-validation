import pandas as pd
from sklearn.metrics import confusion_matrix

pipeline = pd.read_csv('blah.tab', delimiter='\t')
simulated = pd.read_csv('AF2122.refseq2simseq.map.txt', delimiter='\t')

simulated_pos = set(simulated['ref_start'].values)
pipeline_pos = set(pipeline['POS'].values)

missed = simulated_pos - pipeline_pos


# TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
tp = len(simulated_pos.intersection(pipeline_pos))

# FP (the pipeline calls a variant that is not in the simulated genome),
fp = len(pipeline_pos - simulated_pos)

# FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
fn =  len(simulated_pos - pipeline_pos)


#### Stats
# precision (positive predictive value) of each pipeline as TP/(TP + FP), 
precision = tp / (tp+fp)

# recall (sensitivity) as TP/(TP + FN)
sensitivity = tp / (tp+fn)

# miss rate as FN/(TP + FN)
miss_rate = fn / (tp+fn)

#  total number of errors (FP + FN) per million sequenced bases
total_errors = fp+fn


### Output

print("tp", tp)
print("fp", fp)
print("fn", fn)


print("precision", precision)
print("sensitivity", sensitivity)
print("miss_rate", miss_rate)
print("total_errors", total_errors)

a = 1