from globalAlign import main as globalAlign
from localAlign import main as localAlign
from viterbi import main as viterbi
from upgma import main as upgma
import os, sys

def welcome():
	print "\t\t\t\t\t\t#################################################"
	print "\t\t\t\t\t\t#\t\t\t\t\t\t#"
	print "\t\t\t\t\t\t#\t Welcome to Bio Application V.1.0.\t#"
	print "\t\t\t\t\t\t#\t\t\t\t\t\t#"
	print "\t\t\t\t\t\t#################################################"

def menu():
	print "1 - EBI Web Service Access"
	print "2 - Hidden Markov Model"
	print "3 - Phylogenetic Tree"
	option = input("Select an option: ")
	if option == 1:
		print "1a - Global Alignment Tool"
		print "1b - Local Alignment Tool"
		option2 = raw_input("Select an option: ")
		if option2 == "1a":
			globalAlign()
		elif option2 == "1b":
			localAlign()
	if option == 2:
		viterbi()
	if option == 3:
		print "3a - UPGMA"
		option2 = raw_input("Select an option: ")
		if option2 == "3a":
			upgma()

if __name__ == "__main__":
	os.system('cls' if os.name == 'nt' else 'clear')
	welcome()
	menu()
	sys.exit(0)
