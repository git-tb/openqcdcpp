#%%

import numpy as np
# %%

"""
	Generate string replacement rules for Harmonic Polylogs from the Mathematica output of the exact coefficient functions
	for pasting them into the C++ program
"""

lengths = [1,2,3,4]
idcs = [-1,0,1]

for l in lengths:
	poss = np.arange(l)
	idxtuple = np.zeros(l)
	for j in range(3**l):
		HPLstring = "\"HPL"+str(l)+"("
		for p in poss:
			idxtuple[p] = ((j//3**(l-p-1))%3 - 1)
			HPLstring += "%d,"%(idxtuple[p])
		HPLstring += "x)\""
		HPLrepstring = HPLstring.replace("(","_").replace(")","").replace("-","m").replace(",","").replace("x","")
		print(HPLstring+"->"+HPLrepstring+',')

for l in lengths:
	poss = np.arange(l)
	idxtuple = np.zeros(l)
	for j in range(3**l):
		HPLstring = "HPL"+str(l)+"("
		for p in poss:
			idxtuple[p] = ((j//3**(l-p-1))%3 - 1)
			HPLstring += "%d,"%(idxtuple[p])
		HPLstring += "x)"
		HPLrepstring = HPLstring.replace("(","_").replace(")","").replace("-","m").replace(",","").replace("x","")
		print("double "+HPLrepstring+" = "+HPLstring)
# %%
