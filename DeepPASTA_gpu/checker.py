import pandas as pd
import sys

args = sys.argv

report1 = pd.read_csv(args[1],sep='\t',names=['name','vals'])
report2 = pd.read_csv(args[2],sep='\t',names=['name','vals'])

rep_dict = report1['vals'].to_dict()

for i in report2.index:
	key = i % (len(report1)-1)
	val = report2.loc[i,'vals']
	dict_val = rep_dict[key]
	if dict_val != val :
		print(i, key)
		print(val, dict_val)
print('Reports are the same')
