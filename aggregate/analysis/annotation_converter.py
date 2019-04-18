import pandas as pd


def to_df(filename):
	cols = ['image','aggregates']
	#df = pd.DataFrame(columns = cols)
	with open(filename) as file:
		for line in file:
			linelist = line.split()
			if (linelist[2] == 'None}'):
				#temp = pd.DataFrame([[linelist[0],0]],columns = cols)
				#df.append(temp)
			else:
				#temp = pd.DataFrame([[linelist[0],(len(linelist)-2)/2]],columns = cols)
				#df.append(temp)
			print(df)
	return df