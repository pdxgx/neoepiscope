#(chrom, pos, ref, alt)
def reqReadLength(mutList):
	infoList = []
	mutLocs = set()
	lastPos = int(mutList[0][1])
	stInd = lastPos - 32
	for mutInfo in mutList:
		mutPos = int(mutInfo[1])
		if((mutPos-lastPos)<33): endInd = mutPos+32
		else:
			infoList.append((stInd, endInd, mutLocs))
			mutLocs = set()
			stInd = mutPos-32
			endInd = mutPos+32
		mutLocs.add((mutPos-stInd, mutInfo[3]))
		lastPos = mutPos
	infoList.append((stInd, endInd, mutLocs))
	return infoList


