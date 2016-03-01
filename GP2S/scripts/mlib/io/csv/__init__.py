import pylab as PL
import numpy as N
import csv   as CSV
CSV.field_size_limit(1000000000)


class CSVFile:
	"""CSV File class managing column and row headers as well as data type conversion for CSV bulk reading"""
	def __init__(self,file,delimiter=',',dtype='float',rowHeader=True,colHeader=True):
		self._M = None
		self._colH = None
		self._rowH = None
		#1.read raw str datamatrix
		M = readCSV(file,delimiter,typeC='str')
		#2.extrac row/col header if required	
		c0=1*colHeader or 0
		r0=1*rowHeader or 0
		if(rowHeader):
			#extract it
			self._rowH = M[0,c0::]
		if(colHeader):
			self._colH = M[r0::,0]
		#3. inner Matrix
		self._M = N.array(M[r0::,c0::],dtype=dtype)
		pass
		
	def getM(self):
		return self._M
	def getColHeader(self):
		return self._colH
	def getRowHeader(self):
		return self._rowH	

		


def readCSV(file,delimiter=',',typeC='float',toArray=True):
	if isinstance(file,str):
		file = open(file,'rU')
	R = CSV.reader(file,delimiter=delimiter)
	X = []
	for row in R:
		if(typeC=='float'):
			X.append([float(t) for t in row])
		if(typeC=='str'):
			X.append([str(t) for t in row])
	#enforce identical number of columns per row:
	n_col = len(X[0])
	for i in range(len(X)):
		X[i] = X[i][0:n_col]
	if(toArray):
		X = N.array(X)
	return X


def writeCSV(fileName,array,delimiter=','):
	writer = CSV.writer(open(fileName, "wb"))
	writer.writerows(array)
	#PL.save(file, array, fmt='%.6f', delimiter=',')
