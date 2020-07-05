import pandas as pd

def txt_to_df(filename):
	cat = open(filename,'r')
	str = cat.read()
	
	li = list(str.replace('\n', '  ').split('  '))
	#newList = [li[0], li[1]]
	newList = []
	for i in range(len(li)):
		if not li[i] == '':
			newList.append(li[i])
	li = newList
	print(li)
	
	if filename == 'catalog.txt':
		column = 16
	
	if filename == 'visibility.txt':
		column = 4
	
	return makeTable(li,column)
	
def makeTable(li,n):	
	dict = {}
	for i in range(0, n):
		dict[li[i]] = []
	
	x = 0
	for column in dict:
		for i in range(n+x, len(li), n):
			dict[column].append(li[i])
		if ' ' in dict[column]:
			dict[column].remove(' ')
		x+=1
	
	table = pd.DataFrame(dict)
	return table
