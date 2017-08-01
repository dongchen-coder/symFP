
#Debug info

iter_space = [[1, 1024], [1, 1024]]
checking_space = [[1, 1024], [4,1021]]

#iter_space = [[1, 10240], [1, 10240]]
#checking_space = [[1, 10240], [4,10237]]

DebugRT = []


def CalDistanceVector(source, sink):
	#print 'Distance Vector Calculation:'

	# distance vector calculation
	distVec = [source[0] - sink[0], source[1] - sink[1]]
	
	# Check valid
	validFlag = 1
	for item in distVec:
		if (item > 0):
			break
		elif (item == 0):
			continue
		else:
			validFlag = 0
			break

	if (validFlag == 0):
		distVec = []
		
	return distVec

def ValidIterationRange(distVec, iterSpace):
	#print 'Iteration Range Calculation:'
	validRange = []

	if (distVec != []):
		validRange = [[max(iterSpace[0][0], iterSpace[0][0]-distVec[0]), min(iterSpace[0][1], iterSpace[0][1]-distVec[0])], [max(iterSpace[1][0], iterSpace[1][0]-distVec[1]), min(iterSpace[1][1],iterSpace[1][1]-distVec[1])]]
		

	return validRange


def CalDistanceVectorAlign(source, sink, sourceAlignPos, sinkAlignPos):
	
	if (source == sink and sourceAlignPos == sinkAlignPos):
		return []

	distVec = [source[0] - sink[0], source[1] - sink[1] + sinkAlignPos - sourceAlignPos]

	# Check valid
	validFlag = 1
	for item in distVec:
		if (item > 0):
			break
		elif (item == 0):
			continue
		else:
			validFlag = 0
			break

	if (validFlag == 0):
		distVec = []

	return distVec

def sortDistVecList(distVecList):
	
	sortedDistVecList = []

	for item in distVecList:
		idx = 0
		for i in range(len(sortedDistVecList)):
			if (sortedDistVecList[idx][0][0] > item[0][0]):
				break
			if (sortedDistVecList[idx][0][0] == item[0][0] and sortedDistVecList[idx][0][1] > item[0][1]):
				break
			idx += 1
		sortedDistVecList.insert(idx, item)

	return sortedDistVecList


def CalRT(distVec, source, sink, arrayIdx):

	M = 6
	#rt = (distVec[0] * 1024 + distVec[1]) * M - arrayIdx.index(source) + arrayIdx.index(sink)
	rt = (distVec[0] * (iter_space[1][1]) + distVec[1]) * M - arrayIdx.index(source) + arrayIdx.index(sink)

	if (rt <= 0):
		return rt + M

	return rt

def CalIntersection(range1, range2):

	if (range1[0] > range2[1] or range1[1] < range2[0]):
		return []

	intersect = []
	intersect.extend(range1)
	intersect.extend(range2)
	intersect.sort()
	intersect.pop(3)
	intersect.pop(0)

	return intersect

# calculate rect1 - rect2
def CalRectanglarSub(rect1, rect2):

	if ([] in rect2):
		return rect1

	if (rect2 == rect1):
		return []

	remainRec = []

	iList = []
	jList = []

	#print rect1, rect2
	
	iList.append(rect2[0])
	if (rect1[0][0] < rect2[0][0]):
		iList.append([rect1[0][0], rect2[0][0]-1])
	if (rect1[0][1] > rect2[0][1]):
		iList.append([rect2[0][1]+1, rect1[0][1]])

	jList.append(rect2[1])
	if (rect2[1][0] < rect2[1][0]):
		jList.append([rect1[1][0], rect2[1][0]-1])
	if (rect1[1][1] > rect2[1][1]):
		jList.append([rect2[1][1]+1, rect1[1][1]])

	for i_item in iList:
		for j_item in jList:
			remainRec.append([i_item, j_item])

	remainRec.remove(rect2)

	#print remainRec

	return remainRec


def CalShortestDistVec(sortedCandiateList, source, arrayIdx):

	#print 'sorted Candidate List'
	#print sortedCandiateList


	#SDVstack = [[float('inf'), [1, 1024], [1, 1024]]]
	SDVstack = [[float('inf'), iter_space[0], iter_space[1]]]
	SDVstackResult = []

	for item in sortedCandiateList:
		
		print item[0]
		
		distVec = item[0][0]
		validRangei = item[0][1][0]
		validRangej = item[0][1][1]
		sink = item[1]		
		rt = CalRT(distVec, source, sink, arrayIdx)

		DebugRT.append(rt)
		#print 'rt', rt

		#print 'Source ', source, ' sink ', sink
		#print distVec, rt

		while (SDVstack):
				
			sdv = SDVstack.pop()
			
			#print 'check', sdv
			#print 'here', sdv, validRangei
			intersecti = CalIntersection(sdv[1], validRangei)
			intersectj = CalIntersection(sdv[2], validRangej)
			if (intersecti != [] and intersectj != []):
				SDVstackResult.append([min(rt, sdv[0]), intersecti, intersectj])
				# add the rest of sdv to sdv stack result
				remainRects = CalRectanglarSub([sdv[1], sdv[2]], [intersecti, intersectj])
				for remainRectRange in remainRects:
					remainRect = [sdv[0]]
					remainRect.append(remainRectRange[0])
					remainRect.append(remainRectRange[1])
					SDVstackResult.append(remainRect)
			else:
				SDVstackResult.append(sdv)

		#print 'temp', SDVstackResult
		#print ''			

		SDVstack = list(SDVstackResult)
		SDVstackResult = []
			
	return SDVstack

def CalRTHistogram(SDVstack):

	rtHisto = {}
	for sdv in SDVstack:
		if (sdv[0] in rtHisto.keys()):
			rtHisto[sdv[0]] += (sdv[1][1]-sdv[1][0] + 1) * (sdv[2][1]-sdv[2][0] + 1) 
		else:
			rtHisto[sdv[0]] = (sdv[1][1]-sdv[1][0] + 1) * (sdv[2][1]-sdv[2][0] + 1)
	return rtHisto


def DumpRTHistogram(rtHisto):

	count = 0
	for item in rtHisto:
		# take B into account
		rtHisto[item] = rtHisto[item] / 4	

		if (item == 6):
			#rtHisto[item] += 1024*1024 * 3 / 4
			rtHisto[item] += iter_space[0][1] * iter_space[1][1] * 3 / 4;

		count += rtHisto[item]

	for item in rtHisto:
		print item, rtHisto[item], round(float(rtHisto[item]) / count, 6)

	return


def Debug_checkIntersection(SDVstack):

	for i in range(len(SDVstack)):
		for j in range(i+1, len(SDVstack)):
			itemi = SDVstack[i]
			itemj = SDVstack[j]
			intersecti = CalIntersection(itemi[1], itemj[1])
			intersectj = CalIntersection(itemi[2], itemj[2])
		
			if (intersecti != [] and intersectj != []):
				print 'Error here, have intersection'
				return
	return

arrayIdx = [[0, 0], [0, 1], [0, -1], [-1, 0], [1, 0]]

distVecList = {}

print 'Pre-processing, calculate distance vector and valid range'

for source in arrayIdx:
	for sink in arrayIdx:

		distVec = CalDistanceVector(source, sink)
		#validRange = ValidIterationRange(distVec, [[1, 1024], [1, 1024]])
		#validRange = ValidIterationRange(distVec, [[2, 1023], [4, 1021]])
		#validRange = ValidIterationRange(distVec, [[1, 1024], [4,1021]])
		validRange = ValidIterationRange(distVec, checking_space)

		for sourceAlignPos in [0, 1, 2, 3]:
			for sinkAlignPos in [0, 1, 2, 3]:
				distVecAlign = CalDistanceVectorAlign(source, sink, sourceAlignPos, sinkAlignPos)
				#validRangeAlign = ValidIterationRange(distVecAlign, [[1, 1024], [1, 1024]])
				#validRangeAlign = ValidIterationRange(distVecAlign, [[2, 1023], [4, 1021]])
				#validRangeAlign = ValidIterationRange(distVecAlign, [[1, 1024], [4, 1021]])
				validRangeAlign = ValidIterationRange(distVecAlign, checking_space)

				#if (sourceAlignPos == 0 and source == [0, 0]):
				#	print 'Process array index: Source: '+ str(source) + ' Sink ' + str(sink)
				#	print '   Alignment ', sourceAlignPos, sinkAlignPos
				#	print distVecAlign

				#if ([source, sink, sourceAlignPos, sinkAlignPos] not in distVecList.keys()):
				distVecList[str(source) + ' ' + str(sink) + ' ' + str(sourceAlignPos) + ' ' + str(sinkAlignPos)] = [distVecAlign, validRangeAlign]


print 'Find the shortest distance vector'
ShortestDistVecAll = []

for source in arrayIdx:
	for sourceAlignPos in [0, 1, 2, 3]:

		print 'Process array index: Source: '+str(source) + ' Alignment Position ' + str(sourceAlignPos)

		candidateList = []

		for sink in arrayIdx:
			for sinkAlignPos in [0, 1, 2, 3]:
				if (distVecList[str(source) + ' ' + str(sink) + ' ' + str(sourceAlignPos) + ' ' + str(sinkAlignPos)] != [[], []]): 			
					candidateList.append([distVecList[str(source) + ' ' + str(sink) + ' ' + str(sourceAlignPos) + ' ' + str(sinkAlignPos)], sink])

		sortedCandidateList = sortDistVecList(candidateList)
		#print source, sourceAlignPos		
	
		ShortestDistVec = []	
		ShortestDistVec = CalShortestDistVec(sortedCandidateList, source, arrayIdx)

		print 'Shortest distance result', ShortestDistVec
		
		ShortestDistVecAll.extend(ShortestDistVec)
#		break
#	break

		Debug_checkIntersection(ShortestDistVec)

rtHisto = CalRTHistogram(ShortestDistVecAll)

DumpRTHistogram(rtHisto)
			
DebugRT.sort()

print DebugRT
