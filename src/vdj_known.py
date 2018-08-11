# explores VDJ_known dataset

#import stuff
import find_files
import sys
import numpy as np

file_path = find_files.find_files('../datasets', 'VDJ')[0]
col = 3
contents = open(file_path).read().splitlines()
colname = contents[0].split('\t')[col]
cdr3s= []
for content in contents[1:]:
  parts = content.split('\t')
  cdr3 = parts[col]
  cdr3s.append(cdr3)

ave_len = sum([len(cdr3) for cdr3 in cdr3s])/len(cdr3s)
max_len = max([len(cdr3) for cdr3 in cdr3s])
min_len = min([len(cdr3) for cdr3 in cdr3s])
residue_freq =dict((i,{}) for i in range(1,max_len+1))

for cdr3 in cdr3s:
  for i, res in enumerate(cdr3):
    i = i+1
    if res not in residue_freq[i]:
      residue_freq[i][res] = 1 
    else:
      residue_freq[i][res] += 1
    
top_residue_fractions = []

for position in residue_freq:
  top_residue = sorted(residue_freq[position].items(), key = lambda item: item[-1])[-1]
  top_residue_fraction = [position,top_residue[0], float(top_residue[-1])/len(cdr3s)]
  print(top_residue_fraction)

print('average length %s' % ave_len)
print('min leght %s, max lenght %s' % (min_len, max_len))
print('source %s, column%s: %s' % (file_path, col, colname))
