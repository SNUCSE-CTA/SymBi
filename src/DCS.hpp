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

#ifndef _DCS_H
#define _DCS_H

#include <vector>
#include <algorithm>
#include <queue>
#include "graph.hpp"
#include "checktime.hpp"

using namespace std;

long long numMatches = 0;
long long numLocalMatched;

/**
 * This structure stores information about one query vertex in DCS structure.
 */
struct DCSNode
{
	int numCandidates; // C[u]
	// int* candidates; // use G.verticesByLabel
	int* mark; // mark[v]=0 if D_1[u,v]=0, mark[v]=1 if D_1[u,v]=1 and D_2[u,v]=0, mark[v]=2 if D_2[u,v]=1
	int* parentCount; // N^1_P[u,v]
	int* childCount; // N^2_C[u,v]
	int* validParentCand; // N^1_{u,v}[u_p]
	int* validChildCand; // N^2_{u,v}[u_c]
	int* validNeighborCand; // N^2_{u,v}[u']
	bool* NLFCheck; // NLFCheck[v]=true if <u,v> passes the NLF filter 
};

struct DCSEdge
{
	int u1, v1, u2, v2;
	
	DCSEdge(const int _u1 = -1, const int _v1 = -1, const int _u2 = -1, const int _v2 = -1) : u1(_u1), v1(_v1), u2(_u2), v2(_v2) {}

	DCSEdge(const DCSEdge& other) : u1(other.u1), v1(other.v1), u2(other.u2), v2(other.v2) {}
};

class DCSStructure
{
public:
	DataGraph G;
	QueryGraph Q;

	int numNodes;
	DCSNode* DCSNodes;

	int rootVertex;

	// priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> Qmin;
	// priority_queue<pair<int, int>> Qmax;
	queue<pair<int, int>> Qmin, Qmax; // Qmin for top-down, Qmax for bottom-up

	void allocate();
	void init();

	// TODO: inline function?
	vector<DCSEdge> findDCSChangedEdges(int v1, int v2, int eLabel);
	void insertEdge(vector<DCSEdge>& DCSEdges);
	void TopDownCheckPlus(int u, int v, int uc, int vc);
	void BottomUpCheckPlus(int up, int vp, int u, int v);
	void deleteEdge(vector<DCSEdge>& DCSEdges);
	void TopDownCheckMinus(int u, int v, int uc, int vc);
	void BottomUpCheckMinus(int up, int vp, int u, int v);

	//vector<vector<pair<int, int>>> Cm; // (Vertex, Time)
	long long matchedDuplicateCount;
	long long localMatchedDuplicateCount;
	int** CmVertices; // CmVertices[u] stores C_M(u) 
	int** CmVerticesDuplicateCount; // consider duplicate edge 
	int** CmSizes; // for adaptive matching order, stores E(u)
	int** CmMinNeighbor;
	int* CmCount; // CmCount[u] stores |C_M(u)|

	bool* visitedQuery;
	bool* visitedData;
	int* extendibleQuery; // number of matched neighbor vertices
	int* matchInfo; // matchInfo[u]=v if u is mapped to v


	bool* isIsolated;
	int* isolatedVertices; // [?] pair<pair<int, int>, int>* or compare function
	long long* matchedCountDP;

	// TODO : consider preMapping in backtracking
	long long validMatchCheck(int u, int v, int ub, long long duplicateCount);
	void computeCm(int u, int lvl);
	bool updateMatch(int u, int v, int lvl, int nextu, int prev);
	void restoreMatch(int u, int v, int lvl);
	void backtrackIsolated(int s, int e, int idx, int cnt, int start);
	void backtrack(int lvl);
	void findMatches(vector<DCSEdge>& DCSEdges);

	bool computeNLF(int u, int v);
};


/**
 * Returns whether <u,v> passes the NLF filter.
 */
bool DCSStructure::computeNLF(int u, int v)
{
	for(int i = 0; i < NLFArraySize; i++)
	{
		if((Q.queryNLFBit[u * NLFArraySize + i] & G.dataNLFBit[v * NLFArraySize + i]) != Q.queryNLFBit[u * NLFArraySize + i])
		{
			return false;
		}
	}
	return true;
}

void DCSStructure::allocate()
{
	// NLF
	NLFBitSize = pairLabelHash.size() * 4;
	NLFArraySize = (NLFBitSize + 63) >> 6;

	Q.queryNLF = new int[Q.numVertices * pairLabelHash.size()]();
	Q.queryNLFBit = new uint64_t[Q.numVertices * NLFArraySize]();

	G.dataNLF = new int[G.numVertices * pairLabelHash.size()]();
	G.dataNLFBit = new uint64_t[G.numVertices * NLFArraySize]();

	// Query NEC
	Q.isLeafNode = new bool[Q.numVertices]();
	Q.NECMapping = new int[Q.numVertices];
	Q.NECInverse = new int[Q.numVertices];
	Q.NECSize = new int[Q.numVertices];

	// DCS
	numNodes = Q.numVertices;
	DCSNodes = new DCSNode[numNodes];
	for(int i = 0; i < numNodes; i++)
	{
		if(Q.preMapping[i] == -1) DCSNodes[i].numCandidates = G.numVerticesByLabel[Q.vLabels[i]];
		else DCSNodes[i].numCandidates = 1;

		DCSNodes[i].mark = new int[DCSNodes[i].numCandidates]();
		DCSNodes[i].parentCount = new int[DCSNodes[i].numCandidates]();
		DCSNodes[i].childCount = new int[DCSNodes[i].numCandidates]();
		DCSNodes[i].validParentCand = new int[DCSNodes[i].numCandidates * Q.dagParent[i].size()]();
		DCSNodes[i].validChildCand = new int[DCSNodes[i].numCandidates * Q.dagChild[i].size()]();
		DCSNodes[i].validNeighborCand = new int[DCSNodes[i].numCandidates * Q.neighbors[i].size()]();
		DCSNodes[i].NLFCheck = new bool[DCSNodes[i].numCandidates]();
	}

	// Backtrack
	visitedQuery = new bool[Q.numVertices]();
	visitedData = new bool[G.numVertices]();
	extendibleQuery = new int[Q.numVertices]();
	matchInfo = new int[Q.numVertices]();
	memset(matchInfo, -1, sizeof(int) * Q.numVertices);

	CmVertices = new int*[Q.numVertices];
	CmVerticesDuplicateCount = new int*[Q.numVertices];
	for(int i = 0; i < Q.numVertices; i++) {
		CmVertices[i] = new int[G.numVertices]();
		CmVerticesDuplicateCount[i] = new int[G.numVertices]();
	}
	CmSizes = new int*[Q.numVertices];
	CmMinNeighbor = new int*[Q.numVertices];
	for(int i = 0; i < Q.numVertices; i++)
	{
		CmSizes[i] = new int[Q.numVertices]();
		CmMinNeighbor[i] = new int[Q.numVertices]();
	}
	CmCount = new int[Q.numVertices]();
	
	isIsolated = new bool[Q.numVertices]();
	isolatedVertices = new int[Q.numVertices]();
	matchedCountDP = new long long[Q.numVertices+1]();
}

void DCSStructure::init()
{
	for(int i = 0; i < Q.numVertices; i++)
	{
		for(int j = 0; j < DCSNodes[i].numCandidates; j++)
		{
			int v = (Q.preMapping[i] == -1) ? G.verticesByLabel[Q.vLabels[i]][j] : Q.preMapping[i];
			DCSNodes[i].NLFCheck[j] = computeNLF(i, v);
		}
	}
	for(int i = 0; i < DCSNodes[rootVertex].numCandidates; i++)
	{
		int v = (Q.preMapping[rootVertex] == -1) ? G.verticesByLabel[Q.vLabels[rootVertex]][i] : Q.preMapping[rootVertex];
		if(DCSNodes[rootVertex].NLFCheck[i])
			Qmin.emplace(Q.dagOrderInv[rootVertex], v);
	}
	vector<DCSEdge> DCSEdges;
	insertEdge(DCSEdges);
}

/**
 * Find E_{DCS} when the edge (v1,v2) is inserted or deleted.
 */
vector<DCSEdge> DCSStructure::findDCSChangedEdges(int v1, int v2, int eLabel)
{
	vector<DCSEdge> DCSEdges;
	for(int u = 0; u < Q.numVertices; u++)
	{
		if(G.vLabels[v1] != Q.vLabels[u]) continue;
		if(Q.preMapping[u] != -1 && Q.preMapping[u] != v1) continue;
		for(auto e : Q.dagChildOutEdges[u])
		{
			if(e.eLabel != eLabel) continue;
			int uc = e.v;
			if(G.vLabels[v2] != Q.vLabels[uc]) continue;
			if(Q.preMapping[uc] != -1 && Q.preMapping[uc] != v2) continue;

			DCSEdges.emplace_back(u, v1, uc, v2);
		}
	}

	for(int u = 0; u < Q.numVertices; u++)
	{
		if(G.vLabels[v2] != Q.vLabels[u]) continue;
		if(Q.preMapping[u] != -1 && Q.preMapping[u] != v2) continue;
		for(auto e : Q.dagChildInEdges[u])
		{
			if(e.eLabel != eLabel) continue;
			int uc = e.v;
			if(G.vLabels[v1] != Q.vLabels[uc]) continue;
			if(Q.preMapping[uc] != -1 && Q.preMapping[uc] != v1) continue;

			DCSEdges.emplace_back(u, v2, uc, v1);
		}
	}

	return DCSEdges;
}

/**
 * Update DCS structure when DCSEdges are inserted
 */
void DCSStructure::insertEdge(vector<DCSEdge>& DCSEdges)
{
	for(auto DCSEdge : DCSEdges) 
	{
		int v1Index = (Q.preMapping[DCSEdge.u1] == -1) ? G.vertexIDByLabel[DCSEdge.v1] : 0;
		int v2Index = (Q.preMapping[DCSEdge.u2] == -1) ? G.vertexIDByLabel[DCSEdge.v2] : 0;

		// NLF
		if(DCSNodes[DCSEdge.u1].mark[v1Index] == 0 && DCSNodes[DCSEdge.u1].parentCount[v1Index] == Q.dagParent[DCSEdge.u1].size() && !DCSNodes[DCSEdge.u1].NLFCheck[v1Index])
		{
			bool newNLFCheck = computeNLF(DCSEdge.u1, DCSEdge.v1);
			if(newNLFCheck)
			{
				DCSNodes[DCSEdge.u1].NLFCheck[v1Index] = true;
				Qmin.emplace(Q.dagOrderInv[DCSEdge.u1], DCSEdge.v1);
				if(DCSNodes[DCSEdge.u1].childCount[v1Index] == Q.dagChild[DCSEdge.u1].size())
				{
					Qmax.emplace(Q.dagOrderInv[DCSEdge.u1], DCSEdge.v1);
				}
			}
		}
		DCSNodes[DCSEdge.u1].NLFCheck[v1Index] = computeNLF(DCSEdge.u1, DCSEdge.v1);
		
		if(DCSNodes[DCSEdge.u2].mark[v2Index] == 0 && DCSNodes[DCSEdge.u2].parentCount[v2Index] == Q.dagParent[DCSEdge.u2].size() && !DCSNodes[DCSEdge.u2].NLFCheck[v2Index])
		{
			bool newNLFCheck = computeNLF(DCSEdge.u2, DCSEdge.v2);
			if(newNLFCheck)
			{
				DCSNodes[DCSEdge.u2].NLFCheck[v2Index] = true;
				Qmin.emplace(Q.dagOrderInv[DCSEdge.u2], DCSEdge.v2);
				if(DCSNodes[DCSEdge.u2].childCount[v2Index] == Q.dagChild[DCSEdge.u2].size())
				{
					Qmax.emplace(Q.dagOrderInv[DCSEdge.u2], DCSEdge.v2);
				}
			}
		}
		DCSNodes[DCSEdge.u2].NLFCheck[v2Index] = computeNLF(DCSEdge.u2, DCSEdge.v2);
		
		// Updated parent case 1
		if(DCSNodes[DCSEdge.u1].mark[v1Index] >= 1) TopDownCheckPlus(DCSEdge.u1, DCSEdge.v1, DCSEdge.u2,DCSEdge.v2);
		if(DCSNodes[DCSEdge.u1].mark[v1Index] == 2)
		{
			int neighborIndex = Q.neighborIndex[DCSEdge.u2][DCSEdge.u1];
			int neighborCandIndex = v2Index * Q.neighbors[DCSEdge.u2].size() + neighborIndex;
			DCSNodes[DCSEdge.u2].validNeighborCand[neighborCandIndex]++;
		}
		// Updated child case 1
		if(DCSNodes[DCSEdge.u2].mark[v2Index] == 2)
		{
			BottomUpCheckPlus(DCSEdge.u1, DCSEdge.v1, DCSEdge.u2, DCSEdge.v2);
			int neighborIndex = Q.neighborIndex[DCSEdge.u1][DCSEdge.u2];
			int neighborCandIndex = v1Index * Q.neighbors[DCSEdge.u1].size() + neighborIndex;
			DCSNodes[DCSEdge.u1].validNeighborCand[neighborCandIndex]++;
		}
	}

	// Updated parent case 2
	while(!Qmin.empty())
	{
		int u = Q.dagOrder[Qmin.front().first];
		int v = Qmin.front().second;
		int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
		Qmin.pop();

		DCSNodes[u].mark[vIndex] = 1;

		Edge pe(-1, -1, -1);
		Graph::const_range range;
		for(auto e : Q.dagChildOutEdges[u])
		{
			int uc = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, true);
			for(auto it = range.first; it != range.second; it++)
			{
				TopDownCheckPlus(u, v, uc, *it);
			}
			pe = e;
		}

		pe.vLabel = pe.eLabel = -1;
		for(auto e : Q.dagChildInEdges[u])
		{
			int uc = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, false);
			for(auto it = range.first; it != range.second; it++)
			{
				TopDownCheckPlus(u, v, uc, *it);
			}
			pe = e;
		}

		if(DCSNodes[u].childCount[vIndex] == Q.dagChild[u].size())
			Qmax.emplace(Q.dagOrderInv[u], v);
	}
	
	// Updated child case 2
	while(!Qmax.empty())
	{
		int u = Q.dagOrder[Qmax.front().first];
		int v = Qmax.front().second;
		int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
		Qmax.pop();

		DCSNodes[u].mark[vIndex] = 2;	

		Edge pe(-1, -1, -1);
		Graph::const_range range;
		for(auto e : Q.dagParentOutEdges[u])
		{
			int up = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, true);
			for(auto it = range.first; it != range.second; it++)
			{
				BottomUpCheckPlus(up, *it, u, v);
			}
			pe = e;
		}

		pe.vLabel = pe.eLabel = -1;
		for(auto e : Q.dagParentInEdges[u])
		{
			int up = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, false);
			for(auto it = range.first; it != range.second; it++)
			{
				BottomUpCheckPlus(up, *it, u, v);
			}
			pe = e;
		}

		for(auto e : Q.outEdges[u])
		{
			int un = e.v;
			range = G.getAdj(v, e.vLabel, e.eLabel, true);

			for(auto it = range.first; it != range.second; it++)
			{
				int vn = *it;
				if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
				int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
				int neighborIndex = Q.neighborIndex[un][u];
				int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
				DCSNodes[un].validNeighborCand[neighborCandIndex]++;
			}
		}
		
		for(auto e : Q.inEdges[u])
		{
			int un = e.v;
			range = G.getAdj(v, e.vLabel, e.eLabel, false);

			for(auto it = range.first; it != range.second; it++)
			{
				int vn = *it;
				if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
				int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
				int neighborIndex = Q.neighborIndex[un][u];
				int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
				DCSNodes[un].validNeighborCand[neighborCandIndex]++;
			}
		}
	}
}

/**
 * Update when <u,v> is updated parent of <uc,vc>
 * InsertionTopDown(<u,v>,<uc,vc>)
 */

void DCSStructure::TopDownCheckPlus(int u, int v, int uc, int vc)
{
	if(Q.preMapping[uc] != -1 && Q.preMapping[uc] != vc) return;
	int vcIndex = (Q.preMapping[uc] == -1) ? G.vertexIDByLabel[vc] : 0;
	int parentIndex = Q.dagIndex[uc][u];
	int parentCandIndex = vcIndex * Q.dagParent[uc].size() + parentIndex;

	if(DCSNodes[uc].validParentCand[parentCandIndex] == 0)
	{
		DCSNodes[uc].parentCount[vcIndex]++;
		if(DCSNodes[uc].parentCount[vcIndex] == Q.dagParent[uc].size())
		{
			if(DCSNodes[uc].NLFCheck[vcIndex])
				Qmin.emplace(Q.dagOrderInv[uc], vc);
		}
	}
	DCSNodes[uc].validParentCand[parentCandIndex]++;
}

/**
 * Update when <u,v> is updated child of <up,vp>
 * InsertionBottomUp(<u,v>,<up,vp>)
 */
void DCSStructure::BottomUpCheckPlus(int up, int vp, int u, int v)
{
	if(Q.preMapping[up] != -1 && Q.preMapping[up] != vp) return;
	int vpIndex = (Q.preMapping[up] == -1) ? G.vertexIDByLabel[vp] : 0;
	int childIndex = Q.dagIndex[up][u];
	int childCandIndex = vpIndex * Q.dagChild[up].size() + childIndex;

	if(DCSNodes[up].validChildCand[childCandIndex] == 0)
	{
		DCSNodes[up].childCount[vpIndex]++;
		if(DCSNodes[up].mark[vpIndex] == 1 && DCSNodes[up].childCount[vpIndex] == Q.dagChild[up].size())
		{
			Qmax.emplace(Q.dagOrderInv[up], vp);
		}
	}
	DCSNodes[up].validChildCand[childCandIndex]++;
}

/**
 * Update DCS structure when DCSEdges are deleted
 */
void DCSStructure::deleteEdge(vector<DCSEdge>& DCSEdges)
{
	for(auto DCSEdge : DCSEdges) 
	{
		int v1Index = (Q.preMapping[DCSEdge.u1] == -1) ? G.vertexIDByLabel[DCSEdge.v1] : 0;
		int v2Index = (Q.preMapping[DCSEdge.u2] == -1) ? G.vertexIDByLabel[DCSEdge.v2] : 0;

		// NLF
		if(DCSNodes[DCSEdge.u1].mark[v1Index] >= 1 && DCSNodes[DCSEdge.u1].NLFCheck[v1Index])
		{
			bool newNLFCheck = computeNLF(DCSEdge.u1, DCSEdge.v1);
			if(!newNLFCheck)
			{
				DCSNodes[DCSEdge.u1].NLFCheck[v1Index] = false;
				Qmin.emplace(Q.dagOrderInv[DCSEdge.u1], DCSEdge.v1);
			}
		}
		DCSNodes[DCSEdge.u1].NLFCheck[v1Index] = computeNLF(DCSEdge.u1, DCSEdge.v1);

		if(DCSNodes[DCSEdge.u2].mark[v2Index] >= 1 && DCSNodes[DCSEdge.u2].NLFCheck[v2Index])
		{
			bool newNLFCheck = computeNLF(DCSEdge.u2, DCSEdge.v2);
			if(!newNLFCheck)
			{
				DCSNodes[DCSEdge.u2].NLFCheck[v2Index] = false;
				Qmin.emplace(Q.dagOrderInv[DCSEdge.u2], DCSEdge.v2);
			}
		}
		DCSNodes[DCSEdge.u2].NLFCheck[v2Index] = computeNLF(DCSEdge.u2, DCSEdge.v2);
		
		// Updated parent case 1
		if(DCSNodes[DCSEdge.u1].mark[v1Index] >= 1) TopDownCheckMinus(DCSEdge.u1, DCSEdge.v1, DCSEdge.u2, DCSEdge.v2);
		if(DCSNodes[DCSEdge.u1].mark[v1Index] == 2)
		{
			int neighborIndex = Q.neighborIndex[DCSEdge.u2][DCSEdge.u1];
			int neighborCandIndex = v2Index * Q.neighbors[DCSEdge.u2].size() + neighborIndex;
			DCSNodes[DCSEdge.u2].validNeighborCand[neighborCandIndex]--;
		}
		// Updated child case 1
		if(DCSNodes[DCSEdge.u2].mark[v2Index] == 2) 
		{
			BottomUpCheckMinus(DCSEdge.u1, DCSEdge.v1, DCSEdge.u2, DCSEdge.v2);
			int neighborIndex = Q.neighborIndex[DCSEdge.u1][DCSEdge.u2];
			int neighborCandIndex = v1Index * Q.neighbors[DCSEdge.u1].size() + neighborIndex;
			DCSNodes[DCSEdge.u1].validNeighborCand[neighborCandIndex]--;
		}
	}

	// Updated parent case 2
	while(!Qmin.empty())
	{
		int u = Q.dagOrder[Qmin.front().first];
		int v = Qmin.front().second;
		int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
		Qmin.pop();

		Edge pe(-1, -1, -1);
		Graph::const_range range;
		if(DCSNodes[u].mark[vIndex] == 2)
		{
			for(auto e : Q.dagParentOutEdges[u])
			{
				int up = e.v;
				if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
					range = G.getAdj(v, e.vLabel, e.eLabel, true);
				for(auto it = range.first; it != range.second; it++)
				{
					BottomUpCheckMinus(up, *it, u, v);
				}
				pe = e;
			}

			pe.vLabel = pe.eLabel = -1;
			for(auto e : Q.dagParentInEdges[u])
			{
				int up = e.v;
				if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
					range = G.getAdj(v, e.vLabel, e.eLabel, false);
				for(auto it = range.first; it != range.second; it++)
				{
					BottomUpCheckMinus(up, *it, u, v);
				}
				pe = e;
			}

			for(auto e : Q.outEdges[u])
			{
				int un = e.v;
				range = G.getAdj(v, e.vLabel, e.eLabel, true);

				for(auto it = range.first; it != range.second; it++)
				{
					int vn = *it;
					if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
					int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
					int neighborIndex = Q.neighborIndex[un][u];
					int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
					DCSNodes[un].validNeighborCand[neighborCandIndex]--;
				}
			}
		
			for(auto e : Q.inEdges[u])
			{
				int un = e.v;
				range = G.getAdj(v, e.vLabel, e.eLabel, false);
	
				for(auto it = range.first; it != range.second; it++)
				{
					int vn = *it;
					if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
					int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
					int neighborIndex = Q.neighborIndex[un][u];
					int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
					DCSNodes[un].validNeighborCand[neighborCandIndex]--;
				}
			}
		}

		for(auto e : Q.dagChildOutEdges[u])
		{
			int uc = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, true);
			for(auto it = range.first; it != range.second; it++)
			{
				TopDownCheckMinus(u, v, uc, *it);
			}
			pe = e;
		}

		pe.vLabel = pe.eLabel = -1;
		for(auto e : Q.dagChildInEdges[u])
		{
			int uc = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, false);
			for(auto it = range.first; it != range.second; it++)
			{
				TopDownCheckMinus(u, v, uc, *it);
			}
			pe = e;
		}
		
		DCSNodes[u].mark[vIndex] = 0;
	}

	// Updated child case 2
	while(!Qmax.empty())
	{
		int u = Q.dagOrder[Qmax.front().first];
		int v = Qmax.front().second;
		int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
		Qmax.pop();
		if(DCSNodes[u].mark[vIndex] != 2) continue;

		Edge pe(-1, -1, -1);
		Graph::const_range range;
		for(auto e : Q.dagParentOutEdges[u])
		{
			int up = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, true);
			for(auto it = range.first; it != range.second; it++)
			{
				BottomUpCheckMinus(up, *it, u, v);
			}
			pe = e;
		}

		pe.vLabel = pe.eLabel = -1;
		for(auto e : Q.dagParentInEdges[u])
		{
			int up = e.v;
			if(!(e.vLabel == pe.vLabel && e.eLabel == pe.eLabel))
				range = G.getAdj(v, e.vLabel, e.eLabel, false);
			for(auto it = range.first; it != range.second; it++)
			{
				BottomUpCheckMinus(up, *it, u, v);
			}
			pe = e;
		}

		for(auto e : Q.outEdges[u])
		{
			int un = e.v;
			range = G.getAdj(v, e.vLabel, e.eLabel, true);

			for(auto it = range.first; it != range.second; it++)
			{
				int vn = *it;
				if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
				int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
				int neighborIndex = Q.neighborIndex[un][u];
				int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
				DCSNodes[un].validNeighborCand[neighborCandIndex]--;
			}
		}
		
		for(auto e : Q.inEdges[u])
		{
			int un = e.v;
			range = G.getAdj(v, e.vLabel, e.eLabel, false);

			for(auto it = range.first; it != range.second; it++)
			{
				int vn = *it;
				if(Q.preMapping[un] != -1 && Q.preMapping[un] != vn) continue;
				int vnIndex = (Q.preMapping[un] == -1) ? G.vertexIDByLabel[vn] : 0;
				int neighborIndex = Q.neighborIndex[un][u];
				int neighborCandIndex = vnIndex * Q.neighbors[un].size() + neighborIndex;
				DCSNodes[un].validNeighborCand[neighborCandIndex]--;
			}
		}

		DCSNodes[u].mark[vIndex] = 1;
	}
}


/**
 * Update when <u,v> is updated parent of <uc,vc>
 */
void DCSStructure::TopDownCheckMinus(int u, int v, int uc, int vc)
{
	if(Q.preMapping[uc] != -1 && Q.preMapping[uc] != vc) return;
	int vcIndex = (Q.preMapping[uc] == -1) ? G.vertexIDByLabel[vc] : 0;
	int parentIndex = Q.dagIndex[uc][u];
	int parentCandIndex = vcIndex * Q.dagParent[uc].size() + parentIndex;

	if(DCSNodes[uc].validParentCand[parentCandIndex] == 1)
	{
		if(DCSNodes[uc].mark[vcIndex] >= 1 && DCSNodes[uc].parentCount[vcIndex] == Q.dagParent[uc].size())
		{
			if(DCSNodes[uc].NLFCheck[vcIndex])
				Qmin.emplace(Q.dagOrderInv[uc], vc);
		}
		DCSNodes[uc].parentCount[vcIndex]--;
	}
	DCSNodes[uc].validParentCand[parentCandIndex]--;
}

/**
 * Update when <u,v> is updated child of <up,vp>
 */
void DCSStructure::BottomUpCheckMinus(int up, int vp, int u, int v)
{
	if(Q.preMapping[up] != -1 && Q.preMapping[up] != vp) return;
	int vpIndex = (Q.preMapping[up] == -1) ? G.vertexIDByLabel[vp] : 0;
	int childIndex = Q.dagIndex[up][u];
	int childCandIndex = vpIndex * Q.dagChild[up].size() + childIndex;

	if(DCSNodes[up].validChildCand[childCandIndex] == 1)
	{
		if(DCSNodes[up].mark[vpIndex] == 2 && DCSNodes[up].childCount[vpIndex] == Q.dagChild[up].size())
		{
			if(DCSNodes[up].NLFCheck[vpIndex])
			{
				Qmax.emplace(Q.dagOrderInv[up], vp);
			}
		}
		DCSNodes[up].childCount[vpIndex]--;
	}
	DCSNodes[up].validChildCand[childCandIndex]--;
}

/**
 * Check if v can belong to C_M(u)
 */
long long DCSStructure::validMatchCheck(int u, int v, int ub, long long duplicateCount)
{
	for(int i = 0; i < Q.outEdges[u].size(); i++)
	{
		auto e = Q.outEdges[u][i];
		int u2 = e.v;
		if(!visitedQuery[u2] || u2 == ub) continue;
		int v2 = (Q.preMapping[u2] == -1)? matchInfo[u2] : Q.preMapping[u2];

		int count = G.getEdgeCount(v, v2, e.eLabel, true);
		if(!count) return 0;
		duplicateCount *= count;
	}
	for(int i = 0; i < Q.inEdges[u].size(); i++)
	{
		auto e = Q.inEdges[u][i];
		int u2 = e.v;
		if(!visitedQuery[u2] || u2 == ub) continue;
		int v2 = (Q.preMapping[u2] == -1)? matchInfo[u2] : Q.preMapping[u2];

		int count = G.getEdgeCount(v, v2, e.eLabel, false);
		if(!count) return 0;
		duplicateCount *= count;
	}
	return duplicateCount;
}

/**
 * Compute C_M(u) when the backtracking level is lvl.
 */
void DCSStructure::computeCm(int u, int lvl)
{
	int ub = CmMinNeighbor[lvl][u]; // u_min
	int vb = matchInfo[ub]; // M(u_min)

	for(int i = 0; i < Q.outEdges[ub].size(); i++)
	{
		auto e = Q.outEdges[ub][i];
		if(e.v != u) continue;

		int size = 0;
		auto range = G.getAdj(vb, e.vLabel, e.eLabel, true);
		for(auto it = range.first; it != range.second; it++)
		{
			int v = *it;
			if(visitedData[v]) continue;
			if(Q.preMapping[u] != -1 && Q.preMapping[u] != v) continue;
			int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
			if(DCSNodes[u].mark[vIndex] != 2) continue;
			
			long long duplicateCount = validMatchCheck(u, v, ub, it.getDuplicateCount());
			if(!duplicateCount) continue;

			CmVertices[u][size] = v;
			CmVerticesDuplicateCount[u][size] = duplicateCount;
			size++;
		}
		CmCount[u] = size;
		return;
	}

	for(int i = 0; i < Q.inEdges[ub].size(); i++)
	{
		auto e = Q.inEdges[ub][i];
		if(e.v != u) continue;

		int size = 0;
		auto range = G.getAdj(vb, e.vLabel, e.eLabel, false);
		for(auto it = range.first; it != range.second; it++)
		{
			int v = *it;
			if(visitedData[v]) continue;
			if(Q.preMapping[u] != -1 && Q.preMapping[u] != v) continue;
			int vIndex = (Q.preMapping[u] == -1) ? G.vertexIDByLabel[v] : 0;
			if(DCSNodes[u].mark[vIndex] != 2) continue;

			long long duplicateCount = validMatchCheck(u, v, ub, it.getDuplicateCount());
			if(!duplicateCount) continue;

			CmVertices[u][size] = v;
			CmVerticesDuplicateCount[u][size] = duplicateCount;
			size++;
		}
		CmCount[u] = size;
		return;
	}
}

/**
 * Map u to v and update backtracking information
 */
bool DCSStructure::updateMatch(int u, int v, int lvl, int nextu = -1, int prev = -1)
{
	if(prev != -1) computeCm(prev, lvl - 1);
	matchInfo[u] = v;
	visitedQuery[u] = true;
	visitedData[v] = true;

	if(lvl > 0)
	{
		for(int i = 0; i < Q.numVertices; i++)
		{
			CmSizes[lvl][i] = CmSizes[lvl-1][i];
			CmMinNeighbor[lvl][i] = CmMinNeighbor[lvl-1][i];
		}
	}

	bool ret = true;
	int vIndex = (Q.preMapping[u] == -1)? G.vertexIDByLabel[v] : 0;
	
	// Update CmSizes, CmMinNeighbor, extendibleQuery, isIsolated
	int iterCount;
	for(iterCount = 0; iterCount < Q.dagParent[u].size(); iterCount++)
	{
		int u2 = Q.dagParent[u][iterCount];
		if(visitedQuery[u2] || Q.NECMapping[u2] != u2) continue;

		int candIndex = vIndex * Q.neighbors[u].size() + Q.neighborIndex[u][u2];
		if(extendibleQuery[u2] == 0 || CmSizes[lvl][u2] > DCSNodes[u].validNeighborCand[candIndex])
		{
			CmSizes[lvl][u2] = DCSNodes[u].validNeighborCand[candIndex];
			CmMinNeighbor[lvl][u2] = u;
			if(CmSizes[lvl][u2] < Q.NECSize[u2]) ret = false;
		}

		extendibleQuery[u2]++;
		if(extendibleQuery[u2] == Q.degree[u2])
		{
			isIsolated[u2] = true;
			if(ret == false) continue;
			if(u2 == nextu)
			{
				if(Q.NECSize[u2] != 1)
				{
					computeCm(u2, lvl);
					if(CmCount[u2] < Q.NECSize[u2]) ret = false;
				}
				continue;
			}
			if(CmSizes[lvl][u2] < Q.NECSize[u2]) CmCount[u2] = CmSizes[lvl][u2];
			else computeCm(u2, lvl);
			if(CmCount[u2] < Q.NECSize[u2]) ret = false;
		}
	}

	for(iterCount = 0; iterCount < Q.dagChild[u].size(); iterCount++)
	{
		int u2 = Q.dagChild[u][iterCount];
		if(visitedQuery[u2] || Q.NECMapping[u2] != u2) continue;

		int candIndex = vIndex * Q.neighbors[u].size() + Q.neighborIndex[u][u2];
		if(extendibleQuery[u2] == 0 || CmSizes[lvl][u2] > DCSNodes[u].validNeighborCand[candIndex])
		{
			CmSizes[lvl][u2] = DCSNodes[u].validNeighborCand[candIndex];
			CmMinNeighbor[lvl][u2] = u;
			if(CmSizes[lvl][u2] < Q.NECSize[u2]) ret = false;
		}

		extendibleQuery[u2]++;
		if(extendibleQuery[u2] == Q.degree[u2])
		{
			isIsolated[u2] = true;
			if(ret == false) continue;
			if(u2 == nextu)
			{
				if(Q.NECSize[u2] != 1)
				{
					computeCm(u2, lvl);
					if(CmCount[u2] < Q.NECSize[u2]) ret = false;
				}
				continue;
			}
			if(CmSizes[lvl][u2] < Q.NECSize[u2]) CmCount[u2] = CmSizes[lvl][u2];
			else computeCm(u2, lvl);
			if(CmCount[u2] < Q.NECSize[u2]) ret = false;
		}
	}
	return ret;
}

/**
 * Unmap u and update backtracking information
 */
void DCSStructure::restoreMatch(int u, int v, int lvl)
{
	visitedQuery[u] = false;
	visitedData[v] = false;
	matchInfo[u] = -1;

	// Update extendibleQuery, isIsolated
	int iterCount;
	for(iterCount = 0; iterCount < Q.dagParent[u].size(); iterCount++)
	{
		int u2 = Q.dagParent[u][iterCount];
		if(visitedQuery[u2] || Q.NECMapping[u2] != u2) continue;

		if(extendibleQuery[u2] == Q.degree[u2])
		{
			isIsolated[u2] = false;
		}
		extendibleQuery[u2]--;
	}

	for(iterCount = 0; iterCount < Q.dagChild[u].size(); iterCount++)
	{
		int u2 = Q.dagChild[u][iterCount];
		if(visitedQuery[u2] || Q.NECMapping[u2] != u2) continue;

		if(extendibleQuery[u2] == Q.degree[u2])
		{
			isIsolated[u2] = false;
		}
		extendibleQuery[u2]--;
	}
}

/**
 * Backtracking method for isolated vertices
 */
void DCSStructure::backtrackIsolated(int s, int e, int idx, int cnt, int start)
{
	if(idx == e)
	{
		numLocalMatched += localMatchedDuplicateCount;
		return;
	}
	if(idx == e-1)
	{
		int u = isolatedVertices[idx];
		for(int i = 0; i <= Q.NECSize[u]; i++)
			matchedCountDP[i] = 0;
		matchedCountDP[0] = localMatchedDuplicateCount;
		
		for(int i = 0; i < CmCount[u]; i++)
		{
			int v = CmVertices[u][i];
			if(visitedData[v]) continue;
			
			for(int j = Q.NECSize[u] - 1; j >= 0; j--)
				matchedCountDP[j+1] += matchedCountDP[j] * CmVerticesDuplicateCount[u][i];
		}
		numLocalMatched += matchedCountDP[Q.NECSize[u]];
	
		return;
	}

	int u = isolatedVertices[idx];
	for(int i = start; i < CmCount[u]; i++)
	{
		int v = CmVertices[u][i];
		if(visitedData[v]) continue;

		visitedData[v] = true;
		long long oldLocalMatchedDuplicateCount = localMatchedDuplicateCount;
		localMatchedDuplicateCount *= CmVerticesDuplicateCount[u][i];
		if(cnt < Q.NECSize[u]) backtrackIsolated(s, e, idx, cnt + 1, i + 1);
		else backtrackIsolated(s, e, idx + 1, 1, 0);
		visitedData[v] = false;
		localMatchedDuplicateCount = oldLocalMatchedDuplicateCount;
	}
}

void DCSStructure::backtrack(int lvl)
{
	if(lvl == Q.numVertices)
	{
		numMatches += matchedDuplicateCount;
		return;
	}

	// Choose next vertex according to matching order
	int u = -1;
	for(int i = 0; i < Q.numVertices; i++)
	{
		// 1. If u is already matched / not extendible / compressed NEC vertices, then pass
		if(visitedQuery[i] || !extendibleQuery[i] || Q.NECMapping[i] != i) continue;

		// 2. If u is isolated and NECSize[u] < |C_M(u)| then pass 
		if(isIsolated[i])
		{
			if(Q.NECSize[i] > CmCount[i]) return;
			else if(Q.NECSize[i] == CmCount[i])
			{
				if(u == -1) u = i;
				else
				{
					if(isIsolated[u] && CmCount[i] < CmCount[u]) u = i;
					else if(!isIsolated[u] && CmCount[i] < CmSizes[lvl-1][u]) u = i;
				}
				continue;
			}
			else continue;
		}

		// 3. Choose u whose |E(u)| is minimum
		if(u == -1 || CmSizes[lvl-1][i] < CmSizes[lvl-1][u]) u = i;

		// 4. Choose u which has the maximum number of non-matched neighbor when weight is same
		if(CmSizes[lvl-1][i] == CmSizes[lvl-1][u])
			if(Q.degree[i] - extendibleQuery[i] > Q.degree[u] - extendibleQuery[u])
				u = i;
	}

	// Only isolated vertices remain
	if(u == -1)
	{
		int size = 0;
		for(int i = 0; i < Q.numVertices; i++)
		{
			if(visitedQuery[i] || Q.NECMapping[i] != i) continue;
			isolatedVertices[size++] = i;
			// computeCm(i, lvl-1);
			int alreadyMatched = 0;
			for(int j = 0; j < CmCount[i]; j++)
			{
				int v = CmVertices[i][j];
				if(visitedData[v]) alreadyMatched++;
			}
			if(CmCount[i] - alreadyMatched < Q.NECSize[i]) return;
		}

		sort(isolatedVertices, isolatedVertices + size, [=](int a, int b)
		{
			if(Q.vLabels[a] != Q.vLabels[b]) return Q.vLabels[a] < Q.vLabels[b];
			return CmCount[a] < CmCount[b];
		});

		long long oldMatchedDuplicateCount = matchedDuplicateCount;
		int s = 0, e = 0;
		while(e < size && !visitedQuery[isolatedVertices[e]])
		{
			while(e < size
				&& !visitedQuery[isolatedVertices[e]]
				&& Q.vLabels[isolatedVertices[s]] == Q.vLabels[isolatedVertices[e]])
			{
				e++;
			}

			numLocalMatched = 0;
			localMatchedDuplicateCount = 1;
			backtrackIsolated(s, e, s, 1, 0);
			if(numLocalMatched == 0)
			{
				matchedDuplicateCount = oldMatchedDuplicateCount;
				return;
			}

			for(int i = s; i < e; i++) numLocalMatched *= factorization(Q.NECSize[isolatedVertices[i]]);
			matchedDuplicateCount *= numLocalMatched;

			s = e;
		}

		numMatches += matchedDuplicateCount;
		matchedDuplicateCount = oldMatchedDuplicateCount;
		return;
	}

	if(isIsolated[u])
	{
		if(CmCount[u] < Q.NECSize[u]) return;
		long long oldMatchedDuplicateCount = matchedDuplicateCount;
		for(int i = 0; i < CmCount[u]; i++)
		{
			int v = CmVertices[u][i];
			if(visitedData[v])
			{
				matchInfo[u] = -1;
				visitedQuery[u] = false;
				for(int j = 0; j < i; j++)
				{
					int v2 = CmVertices[u][j];
					visitedData[v2] = false;
				}
				matchedDuplicateCount = oldMatchedDuplicateCount;
				return;
			}

			matchedDuplicateCount *= CmVerticesDuplicateCount[u][i];
			matchInfo[u] = v;
			visitedQuery[u] = true;
			visitedData[v] = true;
		}

		matchedDuplicateCount *= factorization(Q.NECSize[u]);
		backtrack(lvl);

		for(int i = 0; i < CmCount[u]; i++)
		{
			int v = CmVertices[u][i];
			matchInfo[u] = -1;
			visitedQuery[u] = false;
			visitedData[v] = false;
		}
		matchedDuplicateCount = oldMatchedDuplicateCount;
	}
	else
	{
		computeCm(u, lvl-1);
		if(CmCount[u] < Q.NECSize[u]) return;
		for(int i = 0; i < CmCount[u]; i++) 
		{
			int v = CmVertices[u][i];
			
			long long oldMatchedDuplicateCount = matchedDuplicateCount;
			matchedDuplicateCount *= CmVerticesDuplicateCount[u][i];
			
			bool ret = updateMatch(u, v, lvl);
			if(ret) backtrack(lvl + 1);
			restoreMatch(u, v, lvl);
			
			matchedDuplicateCount = oldMatchedDuplicateCount;
		}
	}
}

void DCSStructure::findMatches(vector<DCSEdge>& DCSEdges)
{
	vector<int> idx;
	for(int i = 0; i < DCSEdges.size(); i++)
	{
		auto DCSEdge = DCSEdges[i];

		int v1Index = (Q.preMapping[DCSEdge.u1] == -1) ? G.vertexIDByLabel[DCSEdge.v1] : 0;
		int v2Index = (Q.preMapping[DCSEdge.u2] == -1) ? G.vertexIDByLabel[DCSEdge.v2] : 0;
		if(DCSNodes[DCSEdge.u1].mark[v1Index] != 2 || DCSNodes[DCSEdge.u2].mark[v2Index] != 2) continue;
		if(Q.NECMapping[DCSEdge.u1] != DCSEdge.u1 || Q.NECMapping[DCSEdge.u2] != DCSEdge.u2) continue;

		idx.push_back(i);
	}
	if(idx.empty()) return;
	
	matchedDuplicateCount = 1;
	int prev;
	bool ret1, ret2;

	for(int i = 0; i < idx.size(); i++) {
		auto DCSEdge = DCSEdges[idx[i]];

		// Handle if u1 or u2 is NEC 
		bool skipped = true;
		if(i == 0 || DCSEdge.u1 != DCSEdges[idx[i-1]].u1 || DCSEdge.v1 != DCSEdges[idx[i-1]].v1) 
		{
			ret1 = updateMatch(Q.NECInverse[DCSEdge.u1], DCSEdge.v1, 0, DCSEdge.u2);
			if(Q.NECInverse[DCSEdge.u1] != DCSEdge.u1)
			{
				matchedDuplicateCount *= Q.NECSize[DCSEdge.u1];
				Q.NECSize[DCSEdge.u1]--;
			}
			skipped = false;
			prev = DCSEdge.u2;
		}

		if(ret1)
		{
			if(skipped) ret2 = updateMatch(Q.NECInverse[DCSEdge.u2], DCSEdge.v2, 1, DCSEdge.u1, prev);
			else ret2 = updateMatch(Q.NECInverse[DCSEdge.u2], DCSEdge.v2, 1);
			if(Q.NECInverse[DCSEdge.u2] != DCSEdge.u2) 
			{
				matchedDuplicateCount *= Q.NECSize[DCSEdge.u2];
				Q.NECSize[DCSEdge.u2]--;
			}
			if(ret2) backtrack(2);
			if(Q.NECInverse[DCSEdge.u2] != DCSEdge.u2)
			{
				matchedDuplicateCount = Q.NECSize[DCSEdge.u1];
				Q.NECSize[DCSEdge.u2]++;
			}
			restoreMatch(Q.NECInverse[DCSEdge.u2], DCSEdge.v2, 1);
		}
			
		if(i + 1 == idx.size() || DCSEdge.u1 != DCSEdges[idx[i+1]].u1 || DCSEdge.v1 != DCSEdges[idx[i+1]].v1) 
		{
			if(Q.NECInverse[DCSEdge.u1] != DCSEdge.u1)
			{
				Q.NECSize[DCSEdge.u1]++;
				matchedDuplicateCount = 1;
			}
			restoreMatch(Q.NECInverse[DCSEdge.u1], DCSEdge.v1, 0);
		}
	}
}

#endif
