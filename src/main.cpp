/** 
 * This code is implementation of the following paper:
 *
 * S. Min, S. G. Park, K. Park, D. Giammarresi, G. F. Italiano, and W.-S. Han.
 * "Symmetric Continuous Subgraph Matching with Bidirectional Dynamic Programming"
 *
 *
 * @date 2021.02.02
 * @author Seunghwan Min, Sung Gwan Park
 */

#include <cstdio>
#include <iostream>
#include "DCS.hpp"

using namespace std;


/**
 * argv0: exe
 * argv1: initial data graph
 * argv2: data graph update stream
 * argv3: query graph
 * argv4: the number of update operations to process(defualt: 1000000007)
 */

int main(int argc, char* argv[])
{
	char* dataGraphPath = argv[1];
	char* streamPath = argv[2];
	char* queryGraphPath = argv[3];
	
	int numQueries = 1000000007;
	if(argc >= 5) numQueries = atoi(argv[4]);


	Timer readDataGraphTimer, readQueryGraphTimer, memoryAllocationTimer, initTimer, queryProcessTimer;
	DCSStructure DCS;

	readQueryGraphTimer.start();
	readQueryGraph(DCS.Q, queryGraphPath);
	readQueryGraphTimer.end();

	readDataGraphTimer.start();
	char ch;
	int leftVertex, rightVertex, eLabel;

	char* inFileBufferPtr = initFile(streamPath);
	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);
			maxDataVertex = max(maxDataVertex, max(abs(leftVertex), abs(rightVertex)));
		}
	}
	clearFile();

	readDataGraph(DCS.G, dataGraphPath);
	readDataGraphTimer.end();

	DCS.rootVertex = DCS.Q.buildDAG();

	memoryAllocationTimer.start();
	DCS.allocate();
	memoryAllocationTimer.end();

	initTimer.start();
	DCS.Q.buildNLF();
	DCS.G.buildNLF();
	DCS.Q.constructLeafNEC();
	DCS.init();
	initTimer.end();
	
	cout << "readDataGraph: " << readDataGraphTimer << "ms" << endl;
	cout << "readQueryGraph: " << readQueryGraphTimer << "ms" << endl;
	cout << "memory allocation: " << memoryAllocationTimer << "ms" << endl;
	cout << "initialize DCS: " << initTimer << "ms" << endl;
	
	inFileBufferPtr = initFile(streamPath);
	queryProcessTimer.start();
	int numUpdate = 0;

	while(*inFileBufferPtr)
	{
		//cout << "numUpdate: " << numUpdate << endl;
		ch = parseChar(&inFileBufferPtr);
		if(!ch) break;

		numUpdate++;
		switch(ch)
		{
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);

			bool insert = (leftVertex >= 0 && rightVertex >= 0);
			if(!insert)
			{
				leftVertex = -leftVertex - 1;
				rightVertex = -rightVertex - 1;
			}

			if(queryEdgeHash.find(make_tuple(DCS.G.vLabels[leftVertex], DCS.G.vLabels[rightVertex], eLabel)) == queryEdgeHash.end()) break;

			if(insert)
			{
				DCS.G.insertEdge(leftVertex, rightVertex, eLabel);
				DCS.G.updateNLF(leftVertex, rightVertex, eLabel, true);
				auto DCSEdges = DCS.findDCSChangedEdges(leftVertex, rightVertex, eLabel);
				DCS.insertEdge(DCSEdges);
				DCS.findMatches(DCSEdges);
			}
			else
			{
				auto DCSEdges = DCS.findDCSChangedEdges(leftVertex, rightVertex, eLabel);
				DCS.findMatches(DCSEdges);
				DCS.G.deleteEdge(leftVertex, rightVertex, eLabel);
				DCS.G.updateNLF(leftVertex, rightVertex, eLabel, false);
				DCS.deleteEdge(DCSEdges);
			}
			break;
		}
		if(numUpdate >= numQueries) break;
	}
	queryProcessTimer.end();

	cout << "Query Process: " << queryProcessTimer << "ms" << endl;
	cout << "numMatches: " << numMatches << endl;
	
	return 0;
}
