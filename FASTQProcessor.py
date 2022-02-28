#Louie Alexander
#This program reads a FASTQ file and stores all sequence and quality data.  User is queried for a sequence of interest, then program reports the amino acids (and counts) that sequence encodes.
#User is queried for a minimal acceptable quality score, then program trims all sequences so that low quality segments are removed.  Trimmed sequences are then output in a new FASTQ file.

def parser(file):
	header = ""
	sequence = ""
	qualityLine = False
	#Accept file name as parameter, then open and read line by line
	with open(file) as fh:
		for line in fh:
			#Remove new line and blank space at end
			line = line.strip()

			if line.startswith("@"):
				#Store "@" line as the sequence header 
				header = line.lstrip("@")

			elif line.startswith("+"):
				#When a quality header is reached, do nothing with that line and flag the next line
				qualityLine= True

			elif qualityLine:
				#After reaching quality data, yield header, sequence, and quality data as tuple and reset
				yield(header, sequence, line)
				sequence = ""
				header = ""
				qualityLine = False

			else:
				#If line doesn't start with '@', '+', and is not quality data, then store as sequence data
				sequence += line

def translationCount(sequence, frameShift):
	#sequence is the DNA sequence, and frameShift is the desired adjustment of the reading frame (0, 1, 2).  
	#Need to convert to RNA sequence by replacing T with U
	rnaSeq = sequence.replace("T", "U")

	#Set up dictionary of all amino acids to keeps counts of amino acid content in sequence
	aminoAcids = {"Methionine":0, "Alanine":0, "Arginine":0, "Asparagine":0, "Aspartate":0, "Cysteine":0, "Glutamate":0, "Glutamine":0, "Glycine":0, "Histidine":0, "Isoleucine":0, "Leucine":0, "Lysine":0, "Lysine":0, "Phenylalanine":0, "Proline":0, "Serine":0, "Threonine":0, "Tryptophan":0, "Tyrosine":0, "Valine":0}
	
	#Shift reading frame based on user input, then read to end of sequence
	for n in range(frameShift, len(rnaSeq), 3):
		#Ensure that there are at least 3 nucleotides left in the sequence to determine a codon
		if n + 2 < len(rnaSeq):
			#Each codon is 3 nucleotides long
			codon = rnaSeq[n] + rnaSeq[n+1] + rnaSeq[n+2]

			if codon in ("AUG"):
				aminoAcids["Methionine"] += 1
			elif codon in ("AUU", "AUC", "AUA"):
				aminoAcids["Isoleucine"] += 1
			elif codon in ("CGG", "CGA", "CGC", "CGU", "AGG", "AGA"):
				aminoAcids["Arginine"] += 1
			elif codon in ("CAA", "CAG"):
				aminoAcids["Glutamine"] += 1
			elif codon in ("CAU", "CAC"):
				aminoAcids["Histidine"] += 1
			elif codon in ("CCC", "CCG", "CCA", "CCU"):
				aminoAcids["Proline"] += 1
			elif codon in ("CUU", "CUG", "CUC", "CUA", "UUA", "UUG"):
				aminoAcids["Leucine"] += 1
			elif codon in ("UGG"):
				aminoAcids["Tryptophan"] += 1
			elif codon in ("UGU", "UGC"):
				aminoAcids["Cysteine"] += 1
			elif codon in ("UAU", "UAC"):
				aminoAcids["Tyrosine"] += 1
			elif codon in ("UCU", "UCC", "UCA", "UCG", "AGU", "AGC"):
				aminoAcids["Serine"] += 1
			elif codon in ("UUU", "UUC"):
				aminoAcids["Phenylalanine"] += 1
			elif codon in ("GGG", "GGC", "GGA", "GGU"):
				aminoAcids["Glycine"] += 1
			elif codon in ("GAA", "GAG"):
				aminoAcids["Glutamate"] += 1
			elif codon in ("GAC", "GAU"):
				aminoAcids["Aspartate"] += 1
			elif codon in ("GCU", "GCC", "GCA", "GCG"):
				aminoAcids["Alanine"] += 1
			elif codon in ("GUU", "GUC", "GUA", "GUG"):
				aminoAcids["Valine"] += 1
			elif codon in ("AAA", "AAG"):
				aminoAcids["Lysine"] += 1
			elif codon in ("AAU", "AAC"):
				aminoAcids["Asparagine"] += 1				
			elif codon in ("ACU", "ACC", "ACA", "ACG"):
				aminoAcids["Threonine"] += 1

	return aminoAcids

def trimmer(listSequenceTuple, quality, fileName):
	#Begin writing new file for trimmed sequences
	with open(fileName, "w") as file:
		#"gene" represents the set of: sequence header, sequence data, and quality data
		for gene in listSequenceTuple:
			#Pull quality data out of each gene tuple
			qualityScores = gene[2]
			avgScore = 0
			cutoffIndex = 0

			#Use sliding-window of 3 bases long to test for the average quality score
			for n in range(0, len(qualityScores), 3):
				#Read through the sequence until you reach the end
				if n+2 <= len(qualityScores):
					quality1 = ord(qualityScores[n]) - 64
					quality2 = ord(qualityScores[n+1]) - 64
					quality3 = ord(qualityScores[n+2]) - 64
					avgScore = (quality1+quality2+quality3)/3

				#Calculate average score for 3 bases, and if average is less than the cutoff, save current index to trim later and stop looking
				if avgScore < quality:
					cutoffIndex = n
					break

			#Recreate fastq file, but use the index [0:n] to only write the high quality part of the sequence.  Bases after the index 'n' are trimmed
			file.write("@" + gene[0] + "\n")
			file.write(gene[1][0:n] + "\n")
			file.write("+" + gene[0] + "\n")
			file.write(qualityScores[0:n] + "\n")

def main():
	while True:
		try:
			#Store the fastq file name for use later
			fileName = input("Please type the name of a FASTQ file you wish to process:\n")
			
			#Make sure user gave valid file name
			if not (fileName.endswith(".fastq") or fileName.endswith(".fq")):
				raise ValueError

		except ValueError:
			print("Please name a file with the .fastq extension.\n")

	#Store all sequence data (tuples) in a list from parser function
	seqList = []

	for n in parser(fileName):
		seqList.append(n)

	print("There are", len(seqList), "total sequences.\n")

	#Exception handling loops through user queries until all inputs are accepted
	while True:
		try:
			#Ask user which sequence they want to look at, then convert that (-1) to an index to use in seqList
			sequenceIndex = int(input("Let's see which amino acids are coded for by a given sequence. Please select your desired sequence by typing a number between 1 and the total number of sequences.\n")) - 1
			readingFrame = input("Please select a reading frame adjustment ('0', '+1', or '+2') for this sequence.\n")

			#If user doesn't want any sequences trimmed, they may select 0
			qualityCutoff = int(input("Please select a quality score cutoff for trimming sequences.  This score must be between 0 and 42; please select 0 if you don't want sequences trimmed.\n"))

			#Flag whether user wants sequences trimmed
			if qualityCutoff == 0:
				trimSeq = False
			else: 
				trimSeq = True
				trimmedName = input("Please type a file name for the trimmed sequences:\n")

			#Raise value errors if incorrect input
			if readingFrame not in ("0", "+1", "+2") or (sequenceIndex < 0 or sequenceIndex >= len(seqList)) or (qualityCutoff > 42 or qualityCutoff < 0):
				raise ValueError

		except ValueError:
			print("When selecting a sequence, please type an integer between 1 and", len(seqList), ". When selecting a reading frame, please select '0', '+1', or '+2'. When selecting a quality score cutoff, please choose an integer between 0 and 42.")
		else:
			break

	#Strip non-integer characters from readingFrame
	if readingFrame in ("+1", "+2"):
		readingFrame.lstrip("+")

	#Make readingFrame an integer for use later
	readingFrame = int(readingFrame)

	#Use sequenceIndex to select the correct tuple out of seqList, then use index [1] to select the string of sequence nucleotides
	aminoAcidDict = translationCount(seqList[sequenceIndex][1], readingFrame)

	print("When reading the sequence '", seqList[sequenceIndex][0], "' with a +", readingFrame, "reading frame, codons encode the following amino acids:")
	
	#Parse through the amino acid dictionary and report each amino acid and their count
	for aminoAcid, count in aminoAcidDict.items():

		#If amino acid is not encoded, then skip it
		if count == 0:
			continue		
		else:
			print(aminoAcid, ":", count, "times")

	if trimSeq:
		#Trim all sequences based on the quality cutoff and write into new file using the requested file name
		if not (trimmedName.endswith(".fastq") or trimmedName.endswith(".fq")):
			trimmedName += ".fastq"

		trimmer(seqList, qualityCutoff, trimmedName)
		print("Your trimmed sequences are available in the file", trimmedName)


if __name__ == '__main__':
	main()