/**
 * This code is implementation of the following paper:
 *
 * S. Min, S. G. Park, K. Park, D. Giammarresi, G. F. Italiano, and W.-S. Han.
 * "Symmetric Continuous Subgraph Matching with Bidirectional Dynamic Programming"
 *
 *
 * @date 2020.07.14
 * @author Seunghwan Min, Sung Gwan Park
 */

#ifndef _GRAPH_H
#define _GRAPH_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>
#include <cstring>
#include "util.hpp"

#define SORT_THRESHOLD 16
using namespace std;

int maxDataVertex = -1, maxQueryVertex = -1;
map<tuple<int, int, int>, int> queryEdgeHash;
map<tuple<int, int, bool>, int> pairLabelHash;
int NLFBitSize, NLFArraySize;

struct Edge
{
	int v; // opposite vertex
	int vLabel; // vertex label
	int eLabel; // edge label
	int duplicateCount; // number of eqaul edges

	Edge(const int _v = -1, const int _vLabel = -1, const int _eLabel = -1, const int _pairLabel = -1, const int _duplicateCount = 1) : v(_v), vLabel(_vLabel), eLabel(_eLabel), duplicateCount(_duplicateCount) {}

	Edge(const Edge& other) : v(other.v), vLabel(other.vLabel), eLabel(other.eLabel), duplicateCount(other.duplicateCount) {}

	bool operator<(const Edge& other) const
	{
		if(eLabel != other.eLabel) return eLabel < other.eLabel;
		else if(vLabel != other.vLabel) return vLabel < other.vLabel;
		else return v < other.v;
	}

	bool operator==(const Edge& other) const
	{
		return (eLabel == other.eLabel && vLabel == other.vLabel && v == other.v);
	}

	bool operator!=(const Edge& other) const
	{
		return !(*this == other);
	}

	Edge& operator=(const Edge& other)
	{
		this->v = other.v;
		this->vLabel = other.vLabel;
		this->eLabel = other.eLabel;
		this->duplicateCount = other.duplicateCount;
		return *this;
	}

	inline int& getVertexID()
	{
		return v;
	}

	inline int& getEdgeLabel()
	{
		return eLabel;
	}

	inline int& getVertexLabel()
	{
		return vLabel;
	}

	inline int& getDuplicateCount()
	{
		return duplicateCount;
	}
};

class Graph
{
public:
	size_t numVertices; // number of vertices in graph
	vector<int> vLabels; // label of vertices

	typedef PODVector<Edge> adjListType; // vector of Edge
	vector<adjListType> outEdges; // outEdges[v] keeps the outgoing edges of v
	vector<adjListType> inEdges; // inEdges[v] keeps the incoming edges of v
	vector<size_t> outSortedSizes; // outSortedSizes[v] keeps the size of sorted region of outEdges[v]
	vector<size_t> inSortedSizes; // inSortedSizes[v] keeps the size of sorted region of inEdges[v]
	
	class EdgeIterator : public iterator<forward_iterator_tag, Edge>
	{
		friend class vector<Edge>;
	public:
		size_t* sortedSize; // points outSortedSizes[v]/inSortedSizes[v] for some v
		adjListType* adjList; // points outEdges[v]/inEdges[v] for some v
		
		typename adjListType::iterator pointer;
		Edge e;

		EdgeIterator() : sortedSize(NULL), adjList(NULL), pointer(NULL) {}

		EdgeIterator(const EdgeIterator& other) : 
			sortedSize(other.sortedSize), adjList(other.adjList), pointer(other.pointer), e(other.e) {}

		EdgeIterator(
				size_t* _size,
				adjListType* _adjList,
				typename adjListType::iterator _pointer,
				Edge& _e) :
			sortedSize(_size), adjList(_adjList), pointer(_pointer), e(_e) {}

		EdgeIterator& operator=(const EdgeIterator& other)
		{
			if(this != &other)
			{
				this->sortedSize = other.sortedSize;
				this->adjList = other.adjList;
				this->pointer = other.pointer;
				this->e = other.e;
			}
			return *this;
		}

		int& operator*()
		{
			return pointer->getVertexID();
		}

		int getDuplicateCount()
		{
			return pointer->getDuplicateCount();
		}

		EdgeIterator& operator++()
		{
			Graph::getNextEdge(++pointer, *sortedSize, *adjList, e);
			return *this;
		}

		EdgeIterator& operator++(int)
		{
			Graph::getNextEdge(++pointer, *sortedSize, *adjList, e);
			return *this;
		}

		EdgeIterator& operator+=(int rhs)
		{
			while(rhs > 0)
			{
				rhs--;
				++(*this);
			}
			return *this;
		}

		bool operator!=(const EdgeIterator& other) const
		{
			return this->pointer != other.pointer || this->adjList != other.adjList || this->e != other.e;
		}

		bool operator==(const EdgeIterator& other) const
		{
			return !(*this!=other);
		}
	};

	typedef EdgeIterator iterator;
	typedef EdgeIterator const_iterator;
	typedef pair<iterator, iterator> range;
	typedef pair<const_iterator, const_iterator> const_range;

	/**
	 * Constructor
	 * @param _numVertices Number of vertices
	 */
	Graph(const size_t _numVertices = 0) : numVertices(_numVertices), outEdges(_numVertices), inEdges(_numVertices), outSortedSizes(_numVertices), inSortedSizes(_numVertices) {}
	

	void setNumVertices(size_t _numVertices)
	{
		numVertices = _numVertices;
		vLabels.resize(numVertices, 0);
		outEdges.resize(numVertices);
		inEdges.resize(numVertices);
		outSortedSizes.resize(numVertices);
		inSortedSizes.resize(numVertices);
	}

	/**
	 * Insert an edge whose label is eLabel from x to y.
	 * If the same edge already exists, duplicateCount of the edge is increased by 1.
	 * @param x Source node
	 * @param y Target node
	 * @param eLabel Edge label
	 */
	void insertEdge(const int x, const int y, const int eLabel)
	{
		insertEdgeOneDirection(x, y, eLabel, true);
		insertEdgeOneDirection(y, x, eLabel, false);
	}
	
	/**
	 * Delete the edge whose label is eLabel from x to y (decrease the duplicateCount of the edge by 1).
	 * @param x Source vertex
	 * @param y Target vertex
	 * @param eLabel Edge label
	 */
	void deleteEdge(const int x, const int y, const int eLabel)
	{
		deleteEdgeOneDirection(x, y, eLabel, true);
		deleteEdgeOneDirection(y, x, eLabel, false);
	}
	
	/**
	 * Returns the range of vertices adjacent to vertexID whose label is vLabel and edge label is eLabel.
	 */
	const_range getAdj(const int vertexID, const int vLabel, const int eLabel, const bool direction)
	{
		adjListType& adjList = direction ? outEdges[vertexID] : inEdges[vertexID];
		size_t& sortedSize = direction ? outSortedSizes[vertexID] : inSortedSizes[vertexID];
		
		// Sort if unsorted region is large
		sortIfDeltaOverThreshold(adjList, sortedSize);

		const_range range;
		Edge e(-1, vLabel, eLabel);
		range.first = EdgeIterator(&sortedSize, &adjList, adjList.begin() - 1, e);
		range.second = EdgeIterator(&sortedSize, &adjList, adjList.end(), e);
		++range.first;
		return range;
	}

	/**
	 * If direction is true, returns the number of edges whose label is eLabel from vertexID1 to vertexID2.
	 * Otherwise, returns the number of edges whose label is eLabel from vertexID2 to vertexID1.
	 */
	int getEdgeCount(const int vertexID1, const int vertexID2, const int eLabel, const bool direction)
	{
		if(direction)
		{
			if(outEdges[vertexID1].size() > inEdges[vertexID2].size()) return getEdgeCount(vertexID2, vertexID1, eLabel, false);
		}
		else
		{
			if(inEdges[vertexID1].size() > outEdges[vertexID2].size()) return getEdgeCount(vertexID2, vertexID1, eLabel, true);
		}

		size_t& sortedSize = (direction ? outSortedSizes[vertexID1] : inSortedSizes[vertexID1]);
		adjListType& adjList = (direction ? outEdges[vertexID1] : inEdges[vertexID1]);

		sortIfDeltaOverThreshold(adjList, sortedSize);

		Edge e(vertexID2, vLabels[vertexID2], eLabel);
		if(unlikely(sortedSize > 2 * CACHE_LINE_SIZE / sizeof(Edge)))
		{
			int i = lower_bound(adjList.begin(), adjList.begin() + sortedSize, e) - adjList.begin();
			if(i == adjList.size()) return 0;
			if(adjList[i] == e) return adjList[i].getDuplicateCount();
			for(i = sortedSize; i < adjList.size(); i++)
			{
				if(adjList[i] == e) return adjList[i].getDuplicateCount();
			}
		}
		else
		{
			for(int i = 0; i < adjList.size(); i++)
			{
				if(adjList[i] == e) return adjList[i].getDuplicateCount();
			}
		}
		return 0;
	}
	

private:
	void insertEdgeOneDirection(const int x, const int y, const int eLabel, const bool direction)
	{
		size_t& sortedSize = (direction ? outSortedSizes[x] : inSortedSizes[x]);
		adjListType& adjList = (direction ? outEdges[x] : inEdges[x]);
		size_t unsortedSize = adjList.size() - sortedSize;
		Edge e(y, vLabels[y], eLabel, 1);

		if(unlikely(sortedSize > 2 * CACHE_LINE_SIZE / sizeof(Edge)))
		{
			int i = lower_bound(adjList.begin(), adjList.begin() + sortedSize, e) - adjList.begin();
			if(i != adjList.size() && adjList[i] == e)
			{
				if(adjList[i].duplicateCount < 0) adjList[i].duplicateCount = 0;
				adjList[i].duplicateCount++;
				return;
			}
		}
		else
		{
			for(int i = 0; i < sortedSize; i++)
			{
				if(adjList[i] == e)
				{
					if(adjList[i].duplicateCount < 0) adjList[i].duplicateCount = 0;
					adjList[i].duplicateCount++;
					return;
				}
			}
		}

		if(unsortedSize <= SORT_THRESHOLD)
		{
			for(int i = sortedSize; i < adjList.size(); i++)
			{
				if(adjList[i] == e)
				{
					if(adjList[i].duplicateCount < 0) adjList[i].duplicateCount = 0;
					adjList[i].duplicateCount++;
					return;
				}
			}
		}
		if(direction) outEdges[x].push_back(e);
		else inEdges[x].push_back(e);
	}

	void deleteEdgeOneDirection(const int x, const int y, const int eLabel, const bool direction)
	{
		size_t& sortedSize = (direction? outSortedSizes[x] : inSortedSizes[x]);
		adjListType& adjList = (direction ? outEdges[x] : inEdges[x]);
		size_t unsortedSize = adjList.size() - sortedSize;
		Edge e(y, vLabels[y], eLabel, -1);

		if(unlikely(sortedSize > 2 * CACHE_LINE_SIZE / sizeof(Edge)))
		{
			int i = lower_bound(adjList.begin(), adjList.begin() + sortedSize, e) - adjList.begin();
			if(i != adjList.size() && adjList[i] == e)
			{
				adjList[i].duplicateCount--;
				return;
			}
		}
		else
		{
			for(int i = 0; i < sortedSize; i++)
			{
				if(adjList[i] == e)
				{
					adjList[i].duplicateCount--;
					return;
				}
			}
		}

		if(unsortedSize <= SORT_THRESHOLD)
		{
			for(int i = sortedSize; i < adjList.size(); i++)
			{
				if(adjList[i] == e)
				{
					adjList[i].duplicateCount--;
					return;
				}
			}
		}
		
		if(direction) outEdges[x].push_back(e);
		else inEdges[x].push_back(e);
	}

	void sortIfDeltaOverThreshold(adjListType& adjList, size_t& sortedSize)
	{
		size_t unsortedSize = adjList.size() - sortedSize;
		if(unlikely(unsortedSize > SORT_THRESHOLD))
		{
			sort(adjList.begin(), adjList.end());
			int now = 0;
			for(int i = 1; i < adjList.size(); i++)
			{
				if(adjList[i] == adjList[now])
				{
					if(adjList[now].duplicateCount < 0) 
					{
						adjList[now].duplicateCount = 0;
					}
					adjList[now].duplicateCount += adjList[i].duplicateCount;
				}
				else
				{
					adjList[++now] = adjList[i];
				}
			}
			adjList.resize(now + 1);
			sortedSize = adjList.size();
		}
	}

	static void getNextEdge(typename adjListType::iterator& begin, const int& sortedSize, adjListType& adjList, Edge& p)
	{
		if(begin == adjList.end())
		{
			return;
		}

		int i;
		const int start = begin - adjList.begin();

		if(p.eLabel != -1)
		{
			if(p.eLabel != adjList[start].eLabel && unlikely((sortedSize - start) > 2 * CACHE_LINE_SIZE / (int)sizeof(Edge)))
			{
				i = lower_bound(begin, adjList.begin() + sortedSize, p) - adjList.begin();
				if(i != adjList.size() && adjList[i].eLabel != p.eLabel) i = sortedSize;
			}
			else
			{
				i = start;
			}
			for(; i < adjList.size(); i++)
			{
				if(adjList[i].eLabel == p.eLabel && (p.vLabel == -1 || p.vLabel == adjList[i].vLabel) && adjList[i].duplicateCount != 0)
				{
					begin = adjList.begin() + i;
					return;
				}
			}
		}
		else
		{
			for(i = start; i < adjList.size(); i++)
			{
				if((p.vLabel == -1 || p.vLabel == adjList[i].vLabel) && adjList[i].duplicateCount != 0)
				{
					begin = adjList.begin() + i;
					return;
				}
			}
		}
		begin = adjList.end();
	}
};

class DataGraph : public Graph
{
public:
	vector<int> numVerticesByLabel; // numVerticesByLabel[l] stores the number of vertices with label is l
	vector<int> vertexIDByLabel; // vertexIDByLabel[v] stores the index of v in verticesByLabel[l] where l is the label of v
	vector<vector<int>> verticesByLabel; // verticesByLabel[l] stores vertices with label l

	int* dataNLF; // dataNLF stores neighbor label frequency
	uint64_t* dataNLFBit;
	void buildNLF();
	void updateNLF(int leftVertex, int rightVertex, int eLabel, bool insert);
};

void readDataGraph(DataGraph& G, char* dataGraphPath)
{
	int numVertices = 0, numEdges = 0;
	int leftVertex, rightVertex, vLabel, eLabel, graphID, temp;
	int maxVLabel = -1;
	char ch;

	char* inFileBufferPtr = initFile(dataGraphPath);
	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'v':
			leftVertex = parseInt(&inFileBufferPtr);
			vLabel = parseInt(&inFileBufferPtr);
			maxDataVertex = max(maxDataVertex, leftVertex);
			maxVLabel = max(maxVLabel, vLabel);

			numVertices++;

			break;
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);
			maxDataVertex = max(maxDataVertex, max(leftVertex, rightVertex));

			break;
		case 't':
			ch = parseChar(&inFileBufferPtr); // #
			graphID = parseInt(&inFileBufferPtr);

			break;
		}
	}
	clearFile();
	
	G.setNumVertices(maxDataVertex + 1);
	numVertices = numEdges = 0;
	inFileBufferPtr = initFile(dataGraphPath);
	
	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'v':
			leftVertex = parseInt(&inFileBufferPtr);
			G.vLabels[leftVertex] = parseInt(&inFileBufferPtr);
			maxVLabel = max(maxVLabel, G.vLabels[leftVertex]);

			numVertices++;

			break;
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);

			if(queryEdgeHash.find(make_tuple(G.vLabels[leftVertex], G.vLabels[rightVertex], eLabel)) != queryEdgeHash.end())
			{
				// Insert edge, but sort later
				// Directly push into internal vector
				Edge e1 = Edge(rightVertex, G.vLabels[rightVertex], eLabel, 1);
				Edge e2 = Edge(leftVertex, G.vLabels[leftVertex], eLabel, 1);
				G.outEdges[leftVertex].push_back(e1);
				G.inEdges[rightVertex].push_back(e2);
			}

			break;
		case 't':
			ch = parseChar(&inFileBufferPtr); // #
			graphID = parseInt(&inFileBufferPtr);

			break;
		}
	}

	G.numVerticesByLabel.resize(maxVLabel + 1, 0);
	G.vertexIDByLabel.resize(G.numVertices);
	G.verticesByLabel.resize(maxVLabel + 1);

	for(int i = 0; i < G.numVertices; i++)
	{
		G.vertexIDByLabel[i] = G.numVerticesByLabel[G.vLabels[i]]++;
		G.verticesByLabel[G.vLabels[i]].push_back(i);
	}

	for(int i = 0; i < G.outEdges.size(); i++)
	{
		sort(G.outEdges[i].begin(), G.outEdges[i].end());
		int p = -1;
		for(int j = 0; j < G.outEdges[i].size(); j++)
		{
			if(p == -1 || G.outEdges[i][p] != G.outEdges[i][j])
			{
				G.outEdges[i][++p] = G.outEdges[i][j];
			}
			else
			{
				G.outEdges[i][p].duplicateCount++;
			}
		}
		G.outEdges[i].resize(p+1);
		G.outSortedSizes[i] = G.outEdges[i].size();
	}

	for(int i = 0; i < G.inEdges.size(); i++)
	{	
		sort(G.inEdges[i].begin(), G.inEdges[i].end());
		int p = -1;
		for(int j = 0; j < G.inEdges[i].size(); j++)
		{
			if(p == -1 || G.inEdges[i][p] != G.inEdges[i][j])
			{
				G.inEdges[i][++p] = G.inEdges[i][j];
			}
			else
			{
				G.inEdges[i][p].duplicateCount++;
			}
		}
		G.inEdges[i].resize(p+1);
		G.inSortedSizes[i] = G.inEdges[i].size();
	}
	clearFile();
}

void DataGraph::buildNLF()
{
	for(int i = 0; i < numVertices; i++)
	{
		for(auto e : outEdges[i])
		{
			int pairLabel = pairLabelHash[make_tuple(e.vLabel, e.eLabel, 1)];
			if(dataNLF[i * pairLabelHash.size() + pairLabel] < 4)
			{
				int b = pairLabel * 4 + dataNLF[i * pairLabelHash.size() + pairLabel];
				dataNLFBit[i * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
			}
			dataNLF[i * pairLabelHash.size() + pairLabel]++;
		}
		for(auto e : inEdges[i])
		{
			int pairLabel = pairLabelHash[make_tuple(e.vLabel, e.eLabel, 0)];
			if(dataNLF[i * pairLabelHash.size() + pairLabel] < 4)
			{
				int b = pairLabel * 4 + dataNLF[i * pairLabelHash.size() + pairLabel];
				dataNLFBit[i * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
			}
			dataNLF[i * pairLabelHash.size() + pairLabel]++;
		}
	}
}

void DataGraph::updateNLF(int leftVertex, int rightVertex, int eLabel, bool insert)
{
	if(insert)
	{
		int pairLabel = pairLabelHash[make_tuple(vLabels[rightVertex], eLabel, 1)];
		if(dataNLF[leftVertex * pairLabelHash.size() + pairLabel] < 4)
		{
			int b = pairLabel * 4 + dataNLF[leftVertex * pairLabelHash.size() + pairLabel];
			dataNLFBit[leftVertex * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
		}
		dataNLF[leftVertex * pairLabelHash.size() + pairLabel]++;

		pairLabel = pairLabelHash[make_tuple(vLabels[leftVertex], eLabel, 0)];
		if(dataNLF[rightVertex * pairLabelHash.size() + pairLabel] < 4)
		{
			int b = pairLabel * 4 + dataNLF[rightVertex * pairLabelHash.size() + pairLabel];
			dataNLFBit[rightVertex * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
		}
		dataNLF[rightVertex * pairLabelHash.size() + pairLabel]++;
	}
	else
	{
		int pairLabel = pairLabelHash[make_tuple(vLabels[rightVertex], eLabel, 1)];
		dataNLF[leftVertex * pairLabelHash.size() + pairLabel]--;
		if(dataNLF[leftVertex * pairLabelHash.size() + pairLabel] < 4)
		{
			int b = pairLabel * 4 + dataNLF[leftVertex * pairLabelHash.size() + pairLabel];
			dataNLFBit[leftVertex * NLFArraySize + (b >> 6)] ^= (1ULL << (b & 0x3F));
		}

		pairLabel = pairLabelHash[make_tuple(vLabels[leftVertex], eLabel, 0)];
		dataNLF[rightVertex * pairLabelHash.size() + pairLabel]--;
		if(dataNLF[rightVertex * pairLabelHash.size() + pairLabel] < 4)
		{
			int b = pairLabel * 4 + dataNLF[rightVertex * pairLabelHash.size() + pairLabel];
			dataNLFBit[rightVertex * NLFArraySize + (b >> 6)] ^= (1ULL << (b & 0x3F));
		}
	}
}

class QueryGraph : public Graph
{
public:
	vector<int> degree;
	vector<int> preMapping; // If u should be mapped to v, then preMapping[u]=v

	void setNumVertices(size_t _numVertices)
	{
		Graph::setNumVertices(_numVertices);
		preMapping.resize(_numVertices);
		degree.resize(_numVertices, 0);
	}

	vector<int> dagOrder;
	vector<int> dagOrderInv;
	vector<int> dagDepth;
	vector<vector<int>> dagParent;
	vector<vector<int>> dagChild;
	vector<vector<int>> dagIndex;
	vector<adjListType> dagParentOutEdges;
	vector<adjListType> dagParentInEdges;
	vector<adjListType> dagChildOutEdges;
	vector<adjListType> dagChildInEdges;

	vector<vector<int>> neighbors;
	vector<vector<int>> neighborIndex;

	int buildDAGWithGivenRoot(int rootVertex);
	int buildDAG();

	int* queryNLF;
	uint64_t* queryNLFBit;
	void buildNLF();

	bool* isLeafNode;
	int* NECMapping;
	int* NECInverse;
	int* NECSize;

	void constructLeafNEC();
};

void readQueryGraph(QueryGraph& Q, char* queryGraphPath)
{
	int numVertices = 0, numEdges = 0;
	int leftVertex, rightVertex, vLabel, eLabel, graphID, temp;
	char ch;

	char* inFileBufferPtr = initFile(queryGraphPath);
	
	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'v':
			leftVertex = parseInt(&inFileBufferPtr);
			vLabel = parseInt(&inFileBufferPtr);
			temp = parseInt(&inFileBufferPtr);

			numVertices++;
			maxQueryVertex = max(maxQueryVertex, leftVertex);
			break;
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);

			numEdges++;
			break;
		case 't':
			ch = parseChar(&inFileBufferPtr); // #
			ch = parseChar(&inFileBufferPtr); // s
			graphID = parseInt(&inFileBufferPtr);

			break;
		}
	}
	clearFile();

	inFileBufferPtr = initFile(queryGraphPath);

	Q.setNumVertices(maxQueryVertex + 1);
	numVertices = numEdges = 0;

	while(*inFileBufferPtr)
	{
		ch = parseChar(&inFileBufferPtr);
		switch(ch)
		{
		case 'v':
			leftVertex = parseInt(&inFileBufferPtr);
			Q.vLabels[numVertices] = parseInt(&inFileBufferPtr);
			Q.preMapping[numVertices] = parseInt(&inFileBufferPtr);

			// Temp code to handle vertex label issue
			// Assume that there is only one type of vertex label in the data graph
			if(Q.vLabels[numVertices] == -1) Q.vLabels[numVertices] = 0; 
			numVertices++;

			break;
		case 'e':
			leftVertex = parseInt(&inFileBufferPtr);
			rightVertex = parseInt(&inFileBufferPtr);
			eLabel = parseInt(&inFileBufferPtr);

			
			if(queryEdgeHash.find(make_tuple(Q.vLabels[leftVertex], Q.vLabels[rightVertex], eLabel)) == queryEdgeHash.end())
			{
				queryEdgeHash[make_tuple(Q.vLabels[leftVertex], Q.vLabels[rightVertex], eLabel)] = numEdges;
				numEdges++;
			}

			if(pairLabelHash.find(make_tuple(Q.vLabels[rightVertex], eLabel, 1)) == pairLabelHash.end())
			{
				int size = pairLabelHash.size();
				pairLabelHash[make_tuple(Q.vLabels[rightVertex], eLabel, 1)] = size;
			}
			if(pairLabelHash.find(make_tuple(Q.vLabels[leftVertex], eLabel, 0)) == pairLabelHash.end())
			{
				int size = pairLabelHash.size();
				pairLabelHash[make_tuple(Q.vLabels[leftVertex], eLabel, 0)] = size;
			}
			Q.insertEdge(leftVertex, rightVertex, eLabel);
			
			Q.degree[leftVertex]++;
			Q.degree[rightVertex]++;

			break;
		case 't':
			ch = parseChar(&inFileBufferPtr); // #
			ch = parseChar(&inFileBufferPtr); // s
			graphID = parseInt(&inFileBufferPtr);

			break;
		}
	}

	for(int i = 0; i < Q.outEdges.size(); i++)
	{
		sort(Q.outEdges[i].begin(), Q.outEdges[i].end());
		Q.outSortedSizes[i] = Q.outEdges[i].size();
	}

	for(int i = 0; i < Q.inEdges.size(); i++)
	{
		sort(Q.inEdges[i].begin(), Q.inEdges[i].end());
		Q.inSortedSizes[i] = Q.inEdges[i].size();
	}

	clearFile();
}

// Returns depth of the rooted DAG whose root is rootVertex
int QueryGraph::buildDAGWithGivenRoot(int rootVertex)
{
	dagOrder.clear();
	dagOrderInv.clear();
	dagDepth.clear();
	dagParent.clear();
	dagChild.clear();
	dagIndex.clear();
	dagParentOutEdges.clear();
	dagParentInEdges.clear();
	dagChildOutEdges.clear();
	dagChildInEdges.clear();
	neighbors.clear();
	neighborIndex.clear();

	dagOrder.resize(numVertices);
	dagOrderInv.resize(numVertices, -1);
	dagDepth.resize(numVertices, -1);
	dagParent.resize(numVertices);
	dagChild.resize(numVertices);
	dagIndex.resize(numVertices);
	dagParentOutEdges.resize(numVertices);
	dagParentInEdges.resize(numVertices);
	dagChildOutEdges.resize(numVertices);
	dagChildInEdges.resize(numVertices);
	for(int i = 0; i < numVertices; i++)
		dagIndex[i].resize(numVertices, -1);

	neighbors.resize(numVertices);
	neighborIndex.resize(numVertices);
	for(int i = 0; i < numVertices; i++)
		neighborIndex[i].resize(numVertices, -1);

	int cnt = 0;
	dagOrderInv[rootVertex] = cnt;
	dagDepth[rootVertex] = 0;
	dagOrder[cnt++] = rootVertex;

	int maxDepth = 0;
	for(int i = 0; i < cnt; i++)
	{
		int currVertex = dagOrder[i];
		maxDepth = dagDepth[currVertex];
		for(auto it : outEdges[currVertex])
		{
			int childVertex = it.v;
			if(dagOrderInv[childVertex] == -1) 
			{
				dagOrderInv[childVertex] = cnt;
				dagDepth[childVertex] = dagDepth[currVertex] + 1;
				dagOrder[cnt++] = childVertex;
			}
			if(dagOrderInv[currVertex] < dagOrderInv[childVertex])
			{
				if(dagIndex[currVertex][childVertex] == -1)
				{
					dagIndex[childVertex][currVertex] = dagParent[childVertex].size();
					dagParent[childVertex].push_back(currVertex);

					dagIndex[currVertex][childVertex] = dagChild[currVertex].size();
					dagChild[currVertex].push_back(childVertex);
				}
				dagChildOutEdges[currVertex].push_back(it);
				Edge e(currVertex, vLabels[currVertex], it.eLabel);
				dagParentInEdges[childVertex].push_back(e);
			}
		}
		for(auto it : inEdges[currVertex])
		{
			int childVertex = it.v;
			if(dagOrderInv[childVertex] == -1) 
			{
				dagOrderInv[childVertex] = cnt;
				dagDepth[childVertex] = dagDepth[currVertex] + 1;
				dagOrder[cnt++] = childVertex;
			}
			if(dagOrderInv[currVertex] < dagOrderInv[childVertex])
			{
				if(dagIndex[currVertex][childVertex] == -1)
				{
					dagIndex[childVertex][currVertex] = dagParent[childVertex].size();
					dagParent[childVertex].push_back(currVertex);

					dagIndex[currVertex][childVertex] = dagChild[currVertex].size();
					dagChild[currVertex].push_back(childVertex);
				}
				dagChildInEdges[currVertex].push_back(it);
				Edge e(currVertex, vLabels[currVertex], it.eLabel);
				dagParentOutEdges[childVertex].push_back(e);
			}
		}
	}
	

	for(int i = 0; i < numVertices; i++)
	{
		sort(dagChildOutEdges[i].begin(), dagChildOutEdges[i].end());
		sort(dagChildInEdges[i].begin(), dagChildInEdges[i].end());
		sort(dagParentOutEdges[i].begin(), dagParentOutEdges[i].end());
		sort(dagParentInEdges[i].begin(), dagParentInEdges[i].end());
		for(auto it : outEdges[i])
		{
			int neighborVertex = it.v;
			
			neighborIndex[i][it.v] = neighbors[i].size();
			neighbors[i].push_back(it.v);

			neighborIndex[it.v][i] = neighbors[it.v].size();
			neighbors[it.v].push_back(i);			
		}
	}
	return maxDepth;
}

// Build the rooted DAG and returns root vertex
int QueryGraph::buildDAG()
{
	int rootVertex = -1;
	int maxDepth = 0;
	
	for(int i = 0; i < numVertices; i++) {
		int depth = buildDAGWithGivenRoot(i);
		if(rootVertex == -1 || maxDepth < depth) 
		{
			rootVertex = i;
			maxDepth = depth;
		}
	}

	buildDAGWithGivenRoot(rootVertex);
	return rootVertex;
}

void QueryGraph::buildNLF()
{
	for(int i = 0; i < numVertices; i++)
	{
		for(auto e : outEdges[i])
		{
			int pairLabel = pairLabelHash[make_tuple(e.vLabel, e.eLabel, 1)];
			if(queryNLF[i * pairLabelHash.size() + pairLabel] < 4)
			{
				int b = pairLabel * 4 + queryNLF[i * pairLabelHash.size() + pairLabel];
				queryNLFBit[i * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
			}
			queryNLF[i * pairLabelHash.size() + pairLabel]++;
		}
		for(auto e : inEdges[i])
		{
			int pairLabel = pairLabelHash[make_tuple(e.vLabel, e.eLabel, 0)];
			if(queryNLF[i * pairLabelHash.size() + pairLabel] < 4)
			{
				int b = pairLabel * 4 + queryNLF[i * pairLabelHash.size() + pairLabel];
				queryNLFBit[i * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
			}
			queryNLF[i * pairLabelHash.size() + pairLabel]++;
		}
	}
}

void QueryGraph::constructLeafNEC()
{
	for(int i = 0; i < numVertices; i++)
	{
		NECMapping[i] = i;
		NECInverse[i] = i;
		NECSize[i] = 1;
	}

	for(int i = 0; i < numVertices; i++)
	{
		if(degree[i] == 1 && !isLeafNode[i])
		{
			isLeafNode[i] = true;
			if(preMapping[i] != -1) continue;
			
			Edge e = (inEdges[i].size() == 1) ? inEdges[i][0] : outEdges[i][0];
			bool direction = (inEdges[i].size() == 1) ? true : false;

			Graph::const_range range = getAdj(e.v, vLabels[i], e.eLabel, direction);
			for(auto it = range.first; it != range.second; it++)
			{
				int v = *it;
				if(degree[v] == 1 && !isLeafNode[v] && preMapping[v] == -1)
				{
					isLeafNode[v] = true;
					NECMapping[v] = i;
					NECInverse[i] = v;
					NECSize[i]++;
				}
			}
		}
	}
}
#endif
