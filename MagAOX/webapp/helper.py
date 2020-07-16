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
	
	if filename == 'static/data/catalog.txt':
		column = 16
	
	if filename == 'static/data/visibility.txt':
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

def getTimeFrame(date): # get starting time and end time for observations on a given day
	#separate date string into year, mont, day (as ints)
	ymd = date.split('-')
	li = []
	for i in ymd:
		li.append(int(i))
	ymd=li
	
	#create datetime object for start time (entered date at 18:00:00) and calculate datetime for end time (7:00:00 next day)
	start=dt.datetime(ymd[0],ymd[1],ymd[2],hour=18)
	end=start + dt.timedelta(hours=13)
	
	#convert datetime objects into strings and put them into a list to return
	return [start.isoformat(' '),end.isoformat(' ')]
