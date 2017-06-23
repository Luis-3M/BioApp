import time, os
from Bio.SubsMat.MatrixInfo import blosum62 as blosum62
from Bio.SubsMat.MatrixInfo import blosum50 as blosum50
import mechanize
import re

# init calc
def initCalc(matrix,width,high):
	for i in xrange(high):
		matrix[i][0] = 0
	for i in xrange(width):
		matrix[0][i] = 0

# init direction
def initDirection(direction,width,high):
	direction[0][0] = "\033[1;32mnone\033[1;m"
	for i in xrange(1,high):
		direction[i][0] = "up"
	for i in xrange(1,width):
		direction[0][i] = "left"

# cell/direction calc
def cellCalc(matrix,direction,width,high,seq1,seq2,gap,opt):
	calc = []
	for i in xrange(1,high): #high
		for j in xrange(1,width): #width
				pair = (seq1[j],seq2[i])
				if opt==1 or opt == 2:
					calc.append(matrix[i-1][j-1]+scoreMatch(pair,blosum62)) # DIAGONAL 0
				elif opt==3:
					calc.append(matrix[i-1][j-1]+scoreMatch(pair,blosum50)) # DIAGONAL 0
				calc.append(matrix[i-1][j]-gap) # UP 1
				calc.append(matrix[i][j-1]-gap) # LEFT 2
				calc.append(0)
				matrix[i][j] = max(calc)
				if(calc.index(max(calc)) == 0):
					direction[i][j] = "diag"
				if (calc.index(max(calc)) == 1):
					direction[i][j] = "up"
				if(calc.index(max(calc)) == 2):
					direction[i][j] = "left"
				calc = []

def getMaxValue(matrix,width,high):
	maxVal = -1
	posI = 0
	posJ = 0
	for i in xrange(0,high): #high
		for j in xrange(0,width): #width
			if(matrix[i][j] > maxVal):
				maxVal = matrix[i][j]
				posI = i
				posJ = j
	return posI, posJ

def scoreMatch(pair,blosum):
	if pair not in blosum:
		return blosum[(tuple(reversed(pair)))]
	else:
		return blosum[pair]

def scorePairwise(alignX,alignY,blosum,gap):
    score = 0
    for i in range(len(alignX)):
        pair = (alignX[i], alignY[i])
        if '-' in pair:
        	score -= gap
        else:
        	score += scoreMatch(pair,blosum)
    return score

def localAlign(matrix,direction,width,high,seq1,seq2,gap,opt):
	calc = []
	path = "start"
	alignX = ""
	alignY = ""
	#FINDING MAX ON MATRIX
	i,j = getMaxValue(matrix,width,high)
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
		if(matrix[i][j]==0):
			break
	print
	alignX = alignX[::-1]
	alignY = alignY[::-1]
	lstX = [alignX[i:i+50] for i in range(0, len(alignX), 50)]
	lstY = [alignY[i:i+50] for i in range(0, len(alignY), 50)]
	if opt==1 or opt == 2:
		score = scorePairwise(alignX, alignY, blosum62, gap)
		print "Matrix: BLOSUM62"
	elif opt==3:
		score = scorePairwise(alignX, alignY, blosum50, gap)
		print "Matrix: BLOSUM50"
	print "Gap: "+str(gap)
	print "Length: "+str(len(alignX))
	print "Score: "+str(score)
	print
	for (lX,lY) in zip(lstX,lstY):
			print "\t\t"+lX
			print "\t\t"+lY
			print
	print 

def waterAlign():
	print "\t\t\t"+"--------------------"			
	print "\t\t\t"+"EBI LOCAL ALIGNMENT"
	print "\t\t\t"+"--------------------"
	file1 = "protein1"
	file2 = "protein2"
	url = "http://www.ebi.ac.uk/Tools/services/web_emboss_water/toolform.ebi"
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
	lst = re.findall("(emboss_water-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w{2}))", filetext)
	jobID = str(lst[0][0])
	os.remove("result")
	time.sleep(5)
	url2 = "http://www.ebi.ac.uk/Tools/services/rest/emboss_water/result/"+jobID+"/aln"
	br.retrieve(url2,'data')[0]
	with open("data", 'r') as fileIn:
		print fileIn.read()
	os.remove("data")

def matcherAlign():
	print "\t\t\t"+"--------------------"			
	print "\t\t\t"+"EBI LOCAL ALIGNMENT"
	print "\t\t\t"+"--------------------"
	file1 = "protein1"
	file2 = "protein2"
	url = "http://www.ebi.ac.uk/Tools/services/web_emboss_matcher/toolform.ebi"
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
	lst = re.findall("(emboss_matcher-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w{2}))", filetext)
	jobID = str(lst[0][0])
	os.remove("result")
	time.sleep(5)
	url2 = "http://www.ebi.ac.uk/Tools/services/rest/emboss_matcher/result/"+jobID+"/aln"
	br.retrieve(url2,'data')[0]
	with open("data", 'r') as fileIn:
		print fileIn.read()
	os.remove("data")

def LAlign():
	print "\t\t\t"+"--------------------"			
	print "\t\t\t"+"EBI LOCAL ALIGNMENT"
	print "\t\t\t"+"--------------------"
	file1 = "protein1"
	file2 = "protein2"
	url = "http://www.ebi.ac.uk/Tools/services/web_lalign/toolform.ebi"
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
	lst = re.findall("(lalign-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w*\d*)-(\w{2}))", filetext)
	jobID = str(lst[0][0])
	os.remove("result")
	time.sleep(5)
	url2 = "http://www.ebi.ac.uk/Tools/services/rest/lalign/result/"+jobID+"/aln"
	br.retrieve(url2,'data')[0]
	with open("data", 'r') as fileIn:
		print fileIn.read()
	os.remove("data")

def main():
	start = time.clock()
	print "\t\t\t"+"---------------"			
	print "\t\t\t"+"LOCAL ALIGNMENT"
	print "\t\t\t"+"---------------"
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
	print "1 - EMBOSS Water\n2 - EMBOSS Matcher\n3 - EMBOSS LALIGN"
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
			initCalc(matrix,width,high)
			initDirection(direction,width,high)
			cellCalc(matrix,direction,width,high,seq1,seq2,gap,opt)
			# LOCAL ALIGN
			localAlign(matrix,direction,width,high,seq1,seq2,gap,opt)
			# EMBOSS WATER
			waterAlign()
			end = time.clock()
			print "Program Running Time: "+str(abs(start-end))+" seconds"
			print
			return
		elif opt == 2:
			gapOpen = 14
			gapExtend = 4
			gap = round((gapOpen*1+gapExtend*6)/7,0)+1
			width = len(seq1)
			high = len(seq2)
			matrix = [[0 for x in xrange(width)] for x in xrange(high)]
			direction = [[" " for x in xrange(width)] for x in xrange(high)]
			# INIT & CALC
			initCalc(matrix,width,high)
			initDirection(direction,width,high)
			cellCalc(matrix,direction,width,high,seq1,seq2,gap,opt)
			# LOCAL ALIGN
			localAlign(matrix,direction,width,high,seq1,seq2,gap,opt)
			# EMBOSS MATCHER
			matcherAlign()
			end = time.clock()
			print "Program Running Time: "+str(abs(start-end))+" seconds"
			print
			return
		elif opt == 3:
			gapOpen = 12
			gapExtend = 2
			gap = round((gapOpen*1+gapExtend*6)/7,0)+1
			width = len(seq1)
			high = len(seq2)
			matrix = [[0 for x in xrange(width)] for x in xrange(high)]
			direction = [[" " for x in xrange(width)] for x in xrange(high)]
			# INIT & CALC
			initCalc(matrix,width,high)
			initDirection(direction,width,high)
			cellCalc(matrix,direction,width,high,seq1,seq2,gap,opt)
			# LOCAL ALIGN
			localAlign(matrix,direction,width,high,seq1,seq2,gap,opt)
			# EMBOSS LALIGN
			LAlign()
			end = time.clock()
			print "Program Running Time: "+str(abs(start-end))+" seconds"
			print
			return
		else:
			print "Choose an EBI Alignment tool: ",
			opt = input()