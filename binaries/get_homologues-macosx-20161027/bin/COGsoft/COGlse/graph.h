#ifndef _GRAPH_H
#define _GRAPH_H

#include<vector>
#include<map>
#include<set>
#include <sstream>
#include "reader.h"

using namespace std;

struct Hit;
struct Node;
struct Edge;
class Graph;

typedef vector<Node*> ConnectedComponent;

typedef map<long,Node*> NodesHash;
typedef map<string,Edge*> EdgesHash;	// hash by qId + "_" + tId


typedef vector<Edge*> EdgesVector;
typedef vector<Hit*> HitsVector;


string convert2id(long id1, long id2)
{
	stringstream ss;
	ss << id1 << "_" << id2;
	return ss.str();
}

//----------------------------------------------
struct Hit
{
	float score;
	bool isValid;
	// Other information about blast hit can be here

	//..........................................
	Hit(float score)
	{
		this->score = score;
		this->isValid = false;
	};
};

//----------------------------------------------
struct Node
{
	long id;
	EdgesVector edges;
	bool isMarked;
	//..........................................
	Node(long id)
	{
		this->id = id;
	}
	//..........................................
	void addEdge(Edge* edge)
	{
		edges.push_back(edge);
	}
};


//----------------------------------------------
struct Edge
{
	Node* qNode;
	Node* tNode;
	Hit* hit;
	//..........................................
	Edge(Node* qNode, Node* tNode, Hit* hit)
	{
		this->qNode = qNode;
		this->tNode = tNode;
		this->hit = hit;
	}
	//..........................................
	string getLabel()
	{
		return convert2id(qNode->id, tNode->id);
	}
	//..........................................
	bool isValid()
	{
		return hit->isValid;
	}
};


//----------------------------------------------
class Graph
{
public:

	NodesHash nodes;
	EdgesHash edges;	
	HitsVector hits;

public:
	//..........................................
	Graph()	
	{
	}
	//..........................................
	~Graph()	
	{
		for(NodesHash::iterator it = nodes.begin(); it != nodes.end(); it++) delete it->second;
		for(EdgesHash::iterator it = edges.begin(); it != edges.end(); it++) delete it->second;
		for(HitsVector::iterator it = hits.begin(); it != hits.end(); it++) delete (*it);
	}
	//..........................................
	void filterByEdgePresence(const char* fileName)
	{
		set<string> fltEdges;

		// load filtered edges
		LineEntryReader reader(fileName);
		FHitEnt fhe;
		while(reader.next(fhe))
		{
			string val = convert2id(fhe.query, fhe.target);
			fltEdges.insert(val);
		}

		// filter out graph edges 
		EdgesHash::iterator it;
		for(it = edges.begin(); it != edges.end(); it++)
		{
			Edge* edge = it->second;
			set<string>::iterator jt = fltEdges.find(edge->getLabel());
			if(jt == fltEdges.end())
			{
				edge->hit->isValid = false;
			}
		}
	}
	//..........................................
	void load(const char* fileName)
	{
		LineEntryReader reader(fileName);
		HitEnt he;
		while( reader.next(he) )
		{
			Node* qNode = getNode(he.query);
			Node* tNode = getNode(he.target);
			
			Hit* hit;
			// Check whether the opposite hit exists
			string key = convert2id(tNode->id,qNode->id);
			EdgesHash::iterator it = edges.find(key);
			if( it == edges.end() )
			{
				hit = new Hit((float)he.score);
				hits.push_back(hit);					
			}
			else
			{
				hit = it->second->hit;
				hit->isValid = true;
			}
		
			Edge* edge = new Edge(qNode, tNode, hit);
			edges[edge->getLabel()] = edge;

			qNode->addEdge(edge);
		}
		reader.close();
	}
	//..........................................
	bool checkEdgesValidity()	
	{
		bool res = true;
		EdgesHash::iterator it;
		for(it = edges.begin(); it != edges.end(); it++)
		{
			Edge* edge = it->second;
			if(!edge->isValid())
			{
				cerr << "Invalid edge! " 
					<< edge->qNode->id 
					<< " - "
					<< edge->tNode->id
					<< endl;
				res = false;
				break;
			}
		}
		return res;
	}
	//..........................................
	Node* getNode(long nodeId)
	{
		Node* node;
		NodesHash::iterator it = nodes.find(nodeId);
		if(it != nodes.end()) node = it->second;
		else
		{
			node = new Node(nodeId);
			nodes[nodeId] = node;
		}	
		return node;
	}
	//..........................................
	void collectConnectedComponents(vector<ConnectedComponent>& conComponents)
	{
		unmarkNodes();
		NodesHash::iterator it;
		for(it = nodes.begin(); it != nodes.end(); it++)
		{
			Node* node = it->second;
			if(node->isMarked) continue;
			ConnectedComponent conComponent;
			collectConnectedComponent(node, conComponent);
			conComponents.push_back(conComponent);			
		}
	};
protected:
	//..........................................
	void unmarkNodes()
	{
		NodesHash::iterator it;
		for(it = nodes.begin(); it != nodes.end(); it++)
		{
			it->second->isMarked = false;
		}
	};
	//..........................................
	void collectConnectedComponent(Node* startNode, ConnectedComponent& conComponent)
	{
		conComponent.push_back(startNode);
		startNode->isMarked = true;

		for (unsigned int i = 0; i < conComponent.size(); i++)
		{
			Node* node = conComponent[i];
			vector<Edge*>::iterator it;
			for(it = node->edges.begin(); it != node->edges.end(); it++)
			{				
				Edge* edge = *it;
				if( !edge->isValid() ) continue;
				if( edge->tNode->isMarked ) continue;
				conComponent.push_back(edge->tNode);
				edge->tNode->isMarked = true;
			}
		}
    }
};


#endif
