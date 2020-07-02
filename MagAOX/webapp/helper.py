import pandas as pd

def txt_to_df(filename):
	cat = open(filename,'r')
	str = cat.read()
	
	li = list(str.replace('\n', '  ').split('  '))
	newList = [li[0], li[1]]
	for i in range(2, len(li)):
		if not li[i] == '':
			newList.append(li[i])
	li = newList
	print(li)
	
	dict = {}
	for i in range(0, 16):
		dict[li[i]] = []
	
	x = 0
	for column in dict:
		for i in range(16+x, len(li), 16):
			dict[column].append(li[i])
		if ' ' in dict[column]:
			dict[column].remove(' ')
		x+=1
	
	table = pd.DataFrame(dict)
	return table.to_html(index = False)