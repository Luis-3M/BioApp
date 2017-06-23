import time

matrix=[[0, 19, 27, 8, 33, 18, 13],
        [19, 0, 31, 18, 36, 1, 13],
        [27, 31, 0, 26, 41, 32, 29],
        [8, 18, 26, 0, 31, 17, 14],
        [33, 36, 41, 31, 0, 35, 28],
        [18, 1, 32, 17, 35, 0, 12],
        [13, 13, 29, 14, 28, 12, 0]]

def showResult(result):
    for i in result:
        print "\t-> Append: " + i[0] + " " + i[1]
    print

def getMinValue(v):
    minVal = float('inf')
    length = len(v)
    for i in xrange(length):
        for j in xrange(length):
            if v[i][j] < minVal and v[i][j] != 0:
                minVal = v[i][j]
                pos = (i,j)
    return minVal, pos

def showMatrix(vec,matrix):
    p='\t'
    for i in vec:
        p+=i
        p+='\t'
    p+='\n'
    length = len(matrix)
    for j in xrange(length):
        p+= vec[j]
        p+='\t'
        for k in matrix[j]:
            p+=str(k)
            p+='\t'
        p+='\n'
    return p

def main():
    start = time.clock()
    result = []
    vec = []
    codID = ord('A')
    length = len(matrix)
    for i in xrange(length):
        vec.append(chr(codID))
        codID+=1
    print showMatrix(vec,matrix)
    while length > 1:
        value,pos = getMinValue(matrix)
        for i in xrange(length):
            if i!=pos[0]:
                matrix[i][pos[0]] = round((matrix[i][pos[0]]+matrix[i][pos[1]])/2.0,2)
                matrix[pos[0]][i] = matrix[i][pos[0]]
        matrix.pop(pos[1])
        for coluna in matrix:
            coluna.pop(pos[1])
        result.append((vec[pos[0]], vec[pos[1]]))
        vec[pos[0]] += vec[pos[1]]
        vec.pop(pos[1])
        print showMatrix(vec,matrix)
        showResult(result)
        print "---------------------------------------------------------"
        result = []
        length-=1
    print
    end = time.clock()
    print "Program Running Time: "+str(abs(start-end))+" seconds"
    print
    return
