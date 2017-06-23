import time,os

def showMatrix(matrix,seq1,seq2,width,high):
	for j in xrange(0,width):	
		if j==0:
			print "|"+"\t",	
		print " | "+seq1[j]+"\t",
	print "|",
	print
	for i in xrange(0,high): #high
		print "| "+seq2[i]+"\t",	
		for j in xrange(0,width): #width
			print " | "+str(matrix[i][j])+"\t",
		print "|"
	print

def showDirection(direction,width,high,seq1,seq2):
	for j in xrange(0,width):
		if j==0:
			print "|"+"\t",	
		print "| "+seq1[j]+"\t",
	print "|",
	print
	for i in xrange(0,high): #high
		print "| "+seq2[i]+"\t",	
		for j in xrange(0,width): #width
			print "| "+str(direction[i][j])+"\t",
		print "|"
	print

def initCalc(M,Ix,Iy,DM,Dx,Dy,h,g,seq1,seq2,width,high):
	M[0][0] = 0
	for i in xrange(high):
		Ix[i][0] = h+(g*i)
		Dx[i][0] = "upIx"
	for j in xrange(width):
		Iy[0][j] = h+(g*j)
		Dy[0][j] = "lIy"
	#---------------------#
	for i in xrange(1,high):
		for j in xrange(1,width):
			calc = []
			if(seq1[j] == seq2[i]):
				calc.append(M[i-1][j-1]+1) #M 0
				calc.append(Ix[i-1][j-1]+1) # Ix 1
				calc.append(Iy[i-1][j-1]+1) #Iy 2
			else:
				calc.append(M[i-1][j-1]-1)
				calc.append(Ix[i-1][j-1]-1)
				calc.append(Iy[i-1][j-1]-1)
			# M MATRIX
			if(calc.index(max(calc)) == 0):
				DM[i][j] = "dM"
			if(calc.index(max(calc)) == 1):
				DM[i][j] = "dIx"
			if(calc.index(max(calc)) == 2):
				DM[i][j] = "dIy"
			M[i][j] = max(calc)
			calc = []
			# Ix MATRIX
			calc.append(M[i-1][j]+h+g)
			calc.append(Ix[i-1][j]+g)
			if(calc.index(max(calc)) == 0):
				Dx[i][j] = "upM"
			if(calc.index(max(calc)) == 1):
				Dx[i][j] = "upIx"
			Ix[i][j] = max(calc)
			calc = []
			# Iy MATRIX
			calc.append(M[i][j-1]+h+g)
			calc.append(Iy[i][j-1]+g)
			if(calc.index(max(calc)) == 0):
				Dy[i][j] = "dM"
			if(calc.index(max(calc)) == 1):
				Dy[i][j] = "lIy"
			Iy[i][j] = max(calc)

def main():
	os.system('cls' if os.name == 'nt' else 'clear')
	start = time.clock()
	print "\t\t"+"-------------------------------"			
	print "\t\t"+"GLOBAL ALIGNMENT AFFINE PENALTY"
	print "\t\t"+"-------------------------------"
	seq1 = '-' + raw_input("Insert the first sequence: ")
	seq2 = '-' + raw_input("Insert the second sequence: ")
	h = input("Gap Opening Penalty: ")
	g = input("Gap Extension Penalty: ")
	print
	width = len(seq1)
	high = len(seq2)
	infinity = float("-inf")
	M = [[infinity for x in xrange(width)] for x in xrange(high)]
	Ix = [[infinity for x in xrange(width)] for x in xrange(high)]
	Iy = [[infinity for x in xrange(width)] for x in xrange(high)]
	DM = [[" " for x in xrange(width)] for x in xrange(high)]
	Dx = [[" " for x in xrange(width)] for x in xrange(high)]
	Dy = [[" " for x in xrange(width)] for x in xrange(high)]
	# INIT CALC
	initCalc(M,Ix,Iy,DM,Dx,Dy,h,g,seq1,seq2,width,high)
	print "M"
	showMatrix(M,seq1,seq2,width,high)
	showDirection(DM,width,high,seq1,seq2)
	print
	print "Ix"
	showMatrix(Ix,seq1,seq2,width,high)
	showDirection(Dx,width,high,seq1,seq2)
	print
	print "Iy"
	showMatrix(Iy,seq1,seq2,width,high)
	showDirection(Dy,width,high,seq1,seq2)

# MAIN #
main()