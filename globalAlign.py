
import time, os
from Bio.SubsMat.MatrixInfo import blosum62 as blosum62
import mechanize
import re

# init calc
def initCalc(matrix,width,high,gap):
	aux = 0
	for i in xrange(high):
		matrix[i][0] = aux
		aux -= gap
	aux = 0
	for i in xrange(width):
		matrix[0][i] = aux
		aux -= gap

# init direction
def initDirection(direction,width,high):
	direction[0][0] = "\033[1;32mnone\033[1;m"
	for i in xrange(1,high):
		direction[i][0] = "up"
	for i in xrange(1,width):
		direction[0][i] = "left"

# cell/direction calc
def cellCalc(matrix,direction,width,high,seq1,seq2,opt,gap):
	calc = []
	for i in xrange(1,high): #high
		for j in xrange(1,width): #width
			if opt==1: # ALTO
				pair = (seq1[j],seq2[i])
				calc.append(matrix[i-1][j]-gap) # UP 0
				calc.append(matrix[i-1][j-1]+scoreMatch(pair,blosum62)) # DIAGONAL 1
				calc.append(matrix[i][j-1]-gap) # LEFT 2
				matrix[i][j] = max(calc)
				if(calc.index(max(calc)) == 0):
					direction[i][j] = "up"
				if (calc.index(max(calc)) == 1):
					direction[i][j] = "diag"
				if(calc.index(max(calc)) == 2):
					direction[i][j] = "left"
				calc = []
			if opt==2: # BAIXO
				pair = (seq1[j],seq2[i])
				calc.append(matrix[i][j-1]-gap) # LEFT 0
				calc.append(matrix[i-1][j-1]+scoreMatch(pair,blosum62)) # DIAGONAL 1
				calc.append(matrix[i-1][j]-gap) # UP 2
				matrix[i][j] = max(calc)
				if(calc.index(max(calc)) == 0):
					direction[i][j] = "left"
				if (calc.index(max(calc)) == 1):
					direction[i][j] = "diag"
				if(calc.index(max(calc)) == 2):
					direction[i][j] = "up"
				calc = []

def scoreMatch(pair,blosum62):
	if pair not in blosum62:
		return blosum62[(tuple(reversed(pair)))]
	else:
		return blosum62[pair]

def scorePairwise(alignX,alignY,blosum62,gap):
    score = 0
    for i in range(len(alignX)):
        pair = (alignX[i], alignY[i])
        if '-' in pair:
        	score -= gap
        else:
        	score += scoreMatch(pair,blosum62)
    return score

def globalAlign(matrix,direction,width,high,seq1,seq2,gap):
	calc = []
	path = "start"
	alignX = ""
	alignY = ""
	i = high-1
	j = width-1
	while 1:
		if(direction[i][j] == "up"):
			path +=" up"
			alignX += "-"
			alignY += seq2[i]
			direction[i][j] = "\033[1;32mup\033[1;m"
			i = i-1
		elif(direction[i][j] == "left"):
			path += " left"
			alignX += seq1[j]
			alignY += "-"
			direction[i][j] = "\033[1;32mleft\033[1;m"
			j = j-1
		elif(direction[i][j] == "diag"):
			path += " diagonal"
			alignX += seq1[j]
			alignY += seq2[i]
			direction[i][j] = "\033[1;32mdiag\033[1;m"
			i = i-1
			j = j-1
		elif(i==0 or j==0):
			break
	print
	alignX = alignX[::-1]
	alignY = alignY[::-1]
	lstX = [alignX[i:i+50] for i in range(0, len(alignX), 50)]
	lstY = [alignY[i:i+50] for i in range(0, len(alignY), 50)]
	score = scorePairwise(alignX, alignY, blosum62, gap)
	print
	print "Matrix: BLOSUM62"
	print "Gap: "+str(gap)
	print "Length: "+str(len(alignX))
	print "Score: "+str(score)
	print
	for (lX,lY) in zip(lstX,lstY):
			print "\t\t"+lX
			print "\t\t"+lY
			print
	print 

def needleAlign():
	print "\t\t\t"+"--------------------"			
	print "\t\t\t"+"EBI GLOBAL ALIGNMENT"
	print "\t\t\t"+"--------------------"
	file1 = "protein1"
	file2 = "protein2"
	url = "http://www.ebi.ac.uk/Tools/services/web_emboss_needle/toolform.ebi"
	br = mechanize.Browser()
	br.set_handle_robots(False) # ignore robots
	br.open(url)
	br.select_form(predicate=lambda f: f.attrs.get('id', None) == 'jd_toolSubmissionForm')
	br.set_all_readonly(False)
	br.form.add_file(open(file1), 'text/plain', file1, name='aupfile')
	br.form.add_file(open(file2), 'text/plain', file2, name='bupfile')
	response = br.submit()
	content = response.read()
	with open("result", "w") as f:
	    f.write(content)
	textfile = open("result", 'r')
	filetext = textfile.read()
	textfile.close()
	lst = re.findall("(emboss_needle-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w{2}))", filetext)
	jobID = str(lst[0][0])
	os.remove("result")
	time.sleep(8)
	url2 = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"+jobID+"/aln"
	br.retrieve(url2,'data')[0]
	with open("data", 'r') as fileIn:
		print fileIn.read()
	os.remove("data")

def stretcherAlign():		
	print "\t\t\t"+"--------------------"			
	print "\t\t\t"+"EBI GLOBAL ALIGNMENT"
	print "\t\t\t"+"--------------------"
	file1 = "protein1"
	file2 = "protein2"
	url = "http://www.ebi.ac.uk/Tools/services/web_emboss_stretcher/toolform.ebi"
	br = mechanize.Browser()
	br.set_handle_robots(False) # ignore robots
	br.open(url)
	br.select_form(predicate=lambda f: f.attrs.get('id', None) == 'jd_toolSubmissionForm')
	br.set_all_readonly(False)
	br.form.add_file(open(file1), 'text/plain', file1, name='aupfile')
	br.form.add_file(open(file2), 'text/plain', file2, name='bupfile')
	response = br.submit()
	content = response.read()
	with open("result", "w") as f:
	    f.write(content)
	textfile = open("result", 'r')
	filetext = textfile.read()
	textfile.close()
	lst = re.findall("(emboss_stretcher-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w{2}))", filetext)
	jobID = str(lst[0][0])
	os.remove("result")
	time.sleep(5)
	url2 = "http://www.ebi.ac.uk/Tools/services/rest/emboss_stretcher/result/"+jobID+"/aln"
	br.retrieve(url2,'data')[0]
	with open("data", 'r') as fileIn:
		print fileIn.read()
	os.remove("data")

def main():
	start = time.clock()
	print "\t\t\t"+"----------------"			
	print "\t\t\t"+"GLOBAL ALIGNMENT"
	print "\t\t\t"+"----------------"
	seq1 = '-'
	seq2 = '-'
	gapOpen = 0
	gapExtend = 0
	gap = 0
	try:
		with open("protein1", "r") as f1:
			seq1 += f1.read()
		with open("protein2", "r") as f2:
			seq2 += f2.read()
	except IOError:
		print "File doesn't exist !"
	opt = input("Highroad = 1 or Lowroad = 2: ")
	while 1:
		if opt==1 or opt==2:
			break
		else:
			opt = input("Wrong choice, please try again: ")
	print
	print "1 - EMBOSS Needle\n2 - EMBOSS Stretcher"
	opt = input("Choose an EBI Alignment tool: ")
	while 1:
		if opt == 1:
			gapOpen = 10
			gapExtend = 0.5
			gap = round((gapOpen*1+gapExtend*6)/7,0)+1
			width = len(seq1)
			high = len(seq2)
			matrix = [[0 for x in xrange(width)] for x in xrange(high)]
			direction = [[" " for x in xrange(width)] for x in xrange(high)]
			# INIT & CALC
			initCalc(matrix,width,high,gap)
			initDirection(direction,width,high)
			cellCalc(matrix,direction,width,high,seq1,seq2,opt,gap)
			# GLOBAL ALIGN
			globalAlign(matrix,direction,width,high,seq1,seq2,gap)
			# EMBOSS NEEDLE
			needleAlign()
			end = time.clock()
			print "Program Running Time: "+str(abs(start-end))+" seconds"
			print
			return
		elif opt == 2:
			gapOpen = 12
			gapExtend = 2
			gap = round((gapOpen*1+gapExtend*6)/7,0)+1
			width = len(seq1)
			high = len(seq2)
			matrix = [[0 for x in xrange(width)] for x in xrange(high)]
			direction = [[" " for x in xrange(width)] for x in xrange(high)]
			# INIT & CALC
			initCalc(matrix,width,high,gap)
			initDirection(direction,width,high)
			cellCalc(matrix,direction,width,high,seq1,seq2,opt,gap)
			# GLOBAL ALIGN
			globalAlign(matrix,direction,width,high,seq1,seq2,gap)
			# EMBOSS STRETCHER
			stretcherAlign()
			end = time.clock()
			print "Program Running Time: "+str(abs(start-end))+" seconds"
			print
			return
		else:
			print "Choose an EBI Alignment tool: ",
			opt = input()
