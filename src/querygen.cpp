#include <stdio.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>

#include "util.hpp"


using namespace std;

class Edge {
public:
	int v1;
	int v2;
	int eLabel;
	bool isStreamEdge;

	Edge(const int _v1 = -1, const int _v2 = -1, const int _eLabel = -1, const bool _isStreamEdge = false) : v1(_v1), v2(_v2), eLabel(_eLabel), isStreamEdge(_isStreamEdge) {}

	bool operator < (const Edge& other) const
	{
		if(v1 != other.v1) return v1 < other.v1;
		else if(v2 != other.v2) return v2 < other.v2;
		else return eLabel < other.eLabel;
	}

	bool operator != (const Edge& other) const
	{
		return (v1 != other.v1 || v2 != other.v2 || eLabel != other.eLabel);
	}
	bool operator == (const Edge& other) const
	{
		return (v1 == other.v1 && v2 == other.v2 && eLabel == other.eLabel);
	}
};

int numELabels = 0;
int numVertices = 0;
vector<int> dataVLabels;
vector<vector<Edge>> dataEdges;
vector<vector<int>> dataEdgesLabels;
vector<Edge> dataEdgesMerged;
void readGraph(FILE* graphFile, bool isStream)
{
	int leftVertex, rightVertex, vLabel, eLabel, graphID;
	char ch;

	int fileSize = getFileSize(graphFile);
	char* inFileBuffer = new char[fileSize + 1];
	fread(inFileBuffer, sizeof(char), fileSize, graphFile);
	inFileBuffer[fileSize] = 0;
	char* inFileBufferPtr = inFileBuffer;

	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'v':
			leftVertex = parseInt(&inFileBufferPtr);
			if(leftVertex >= numVertices)
			{
				dataEdges.resize(leftVertex+1);
				dataEdgesLabels.resize(leftVertex+1);
				dataVLabels.resize(leftVertex+1, 0);
				numVertices = leftVertex+1;
			}

			vLabel = parseInt(&inFileBufferPtr);
			dataVLabels[leftVertex] = vLabel;
			break;

		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			if(leftVertex < 0)
			{
				leftVertex = -leftVertex-1;
				rightVertex = -rightVertex-1;
			}
			
			eLabel = parseInt(&inFileBufferPtr);
			if(max(leftVertex, rightVertex) >= numVertices)
			{
				int newSize = max(leftVertex, rightVertex) + 1;
				dataEdges.resize(newSize);
				dataEdgesLabels.resize(newSize);
				dataVLabels.resize(newSize, 0);
				numVertices = newSize;
			}
			if(eLabel >= numELabels)
			{
				numELabels = eLabel + 1;	
			}

			dataEdgesLabels[leftVertex].emplace_back(eLabel);
			dataEdgesLabels[rightVertex].emplace_back(eLabel);
			dataEdges[leftVertex].emplace_back(leftVertex, rightVertex, eLabel, isStream);
			dataEdges[rightVertex].emplace_back(leftVertex, rightVertex, eLabel, isStream);
			dataEdgesMerged.emplace_back(leftVertex, rightVertex, eLabel, isStream);
			break;

		case 't':
			ch = parseChar(&inFileBufferPtr);
			graphID = parseInt(&inFileBufferPtr);
			break;
		}
	}
}

bool isDuplicate(vector<Edge> &edges, Edge e)
{
	for(auto other: edges)
	{
		if(e.v1 == other.v1 && e.v2 == other.v2) return true;
		if(e.v1 == other.v2 && e.v2 == other.v1) return true;
	}
	return false;
}
int main(int argc, char* argv[])
{
	char* dataGraphPath = argv[1];
	char* streamPath = argv[2];
	char* queryFolderPath = argv[3];
	int queryCount = atoi(argv[4]);
	int querySize = atoi(argv[5]); 

	srand(time(0));

	if(strcmp(streamPath, "None") != 0) {
		FILE* dataGraphFile = fopen(dataGraphPath, "rb");
		readGraph(dataGraphFile, false);
	
		FILE* streamGraphFile = fopen(streamPath, "rb");
		readGraph(streamGraphFile, true);
	}
	else {
		FILE* dataGraphFile = fopen(dataGraphPath, "rb");
		readGraph(dataGraphFile, true);
	}

	for(auto &V : dataEdgesLabels)
	{
		sort(V.begin(), V.end());
		V.erase(unique(V.begin(), V.end()), V.end());
	}

	int edgeCountSum = 0;
	for(int queryNo = 1; queryNo <= queryCount; queryNo++)
	{
		vector <int> queryVertices;
		vector <Edge> queryEdges;
		Edge beginEdge;
		do 
		{
			int idx = rand() % dataEdgesMerged.size();
			beginEdge = dataEdgesMerged[idx];
		} while(!beginEdge.isStreamEdge);

		queryVertices.push_back(beginEdge.v1);
		queryVertices.push_back(beginEdge.v2);
		queryEdges.push_back(beginEdge);

		int currVertex = beginEdge.v2;
		int iterCount = 0;
		while(queryEdges.size() < querySize && iterCount < 10000)
		{
			currVertex = queryVertices[rand() % queryVertices.size()];

			int idx;
			int eLabel = dataEdgesLabels[currVertex][rand() % dataEdgesLabels[currVertex].size()];
			do {
				idx = rand() % dataEdges[currVertex].size();
				break;
			} while(dataEdges[currVertex][idx].eLabel != eLabel);
			
			Edge nextEdge = dataEdges[currVertex][idx];
			int nextVertex = nextEdge.v1 ^ nextEdge.v2 ^ currVertex;

			for(idx = 0; idx < queryEdges.size(); idx++)
			{
				if(queryEdges[idx].v1 == nextEdge.v1 && queryEdges[idx].v2 == nextEdge.v2) break;
				if(queryEdges[idx].v2 == nextEdge.v1 && queryEdges[idx].v1 == nextEdge.v2) break;
			}
			if(idx == queryEdges.size() && nextEdge.v1 != nextEdge.v2)
			{
				queryEdges.push_back(nextEdge);
			}

			if(find(queryVertices.begin(), queryVertices.end(), nextVertex) == queryVertices.end())
			{
				queryVertices.push_back(nextVertex);

				for(auto e: dataEdges[nextVertex])
				{
					if(find(queryVertices.begin(), queryVertices.end(), e.v1) == queryVertices.end()) continue;
					if(find(queryVertices.begin(), queryVertices.end(), e.v2) == queryVertices.end()) continue;
					if(isDuplicate(queryEdges, e)) continue;
					if(e.v1 == e.v2) continue;

					if(rand() % 10 < 9) queryEdges.push_back(e);
				}
			}
			
			currVertex = -1;
			iterCount++;
		}
		random_shuffle(queryVertices.begin(), queryVertices.end());
			
		for(auto &it: queryEdges)
		{
			it.v1 = find(queryVertices.begin(), queryVertices.end(), it.v1) - queryVertices.begin();
			it.v2 = find(queryVertices.begin(), queryVertices.end(), it.v2) - queryVertices.begin();
		}
		if(queryEdges.size() != querySize || queryEdges.size() + 1 <= queryVertices.size()) 
		{
			queryNo--;
			continue;
		}

		char queryGraphPath[105];
		sprintf(queryGraphPath, "%s/Q_%d", queryFolderPath, queryNo);
		FILE* queryGraphFile = fopen(queryGraphPath, "w");
		fprintf(queryGraphFile, "t # s 1\n");
		for(int i = 0; i < queryVertices.size(); i++)
		{
			fprintf(queryGraphFile, "v %d %d -1\n", i, dataVLabels[queryVertices[i]]);
		}
		for(auto &it: queryEdges)
		{
			fprintf(queryGraphFile, "e %d %d %d\n", it.v1, it.v2, it.eLabel);
		}

		printf("G_%d done, queryVertices.size() = %d, queryEdges.size() == %d\n", queryNo, queryVertices.size(), queryEdges.size());
		for(auto it : queryVertices) printf("%d ", it);
		printf("\n");

		edgeCountSum += queryEdges.size();
		fclose(queryGraphFile);
	}
	printf("Average number of edges: %lf\n", (double)edgeCountSum / queryCount);
	return 0;
}
