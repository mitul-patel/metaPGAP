#include "cogmaker.h"
#include "iomanip"
//#include <time.h>

#include <cassert>
#include <deque>
#include <queue>

/*#ifdef DEBUG
	#include <stdio.h>
	#include <time.h>
	clock_t start, end;
	void startClock(){ start = clock(); }
	void finishClock()
	{
		end = clock();
		printf("clocks: %f\n",((double)(end-start))/CLOCKS_PER_SEC);
		cout << (double)(end - start) << endl;
	}
#endif // DEBUG*/

void Filter::setFilter(const int query, map<int, long>* offsets, const string inFileName) 
{
	long offset;
	string prot;
	string str;
	string subject;

	static ifstream inFileF((inFileName + "query2subject.csv").c_str());
	if (!inFileF)
	{
		cerr << "File query2subject.csv not found; filter was not set" << endl;
	}	

	mask.clear();
	
	if ((*offsets).find(query) != (*offsets).end())
	{
		offset = (*(*offsets).find(query)).second;

		inFileF.clear();
		inFileF.seekg(offset);
		while (getline(inFileF, str))
		{
			prot = str.substr(0, str.find_first_of(','));
			if (atoi(prot.c_str()) == query)
			{
				subject = str.substr(str.find_first_of(',') + 1);
				mask.insert(atoi(subject.c_str()));
			}
			else
			{
				break;
			}
		}
	}
	inFileF.clear();
}

bool Filter::isRelevant(const int subj)
{
	if (mask.find(subj) != mask.end() /*|| mask.size() == 0*/) // if there is no mask still process
	{
		return true;
	}
	return false;
}

HitSet::HitSet(	const double	expect, 
				const double	overlap,
				map<string, int>* domHash,
				map<int, string>*	prot2org,
				map<int, int>* prot2len,
				map<int, long>*	offsetsF,
				map<int, long>*	offsetsU,
				map<int, Lse>*	allLse,
				map<int, int>*	prot2lse,
				const string cName,
				const int		startNum,
				const string inFile)
:organism("0")
{
	thresholds.overlap		= overlap;
	thresholds.expect		= expect;
	db.domHash				= domHash;
	db.prot2org				= prot2org;
	db.prot2len				= prot2len;
	db.offsetsF				= offsetsF;
	db.offsetsU				= offsetsU;
	db.allLse				= allLse;
	db.prot2lse				= prot2lse;
	cogName					= cName;
	cogStartNum				= startNum;
	inFileName				= inFile;

	assert(cogStartNum > 0);	// printCOG assumes this!
}

void HitSet::clear()
{
	organism = "0";
	lseHits.clear();
	conflict.clear();
}

void HitSet::insert(set<int>* curLse) 
{
	int curQuery;
	long offset;
	string str;
	vector<string> v;
	int query;
	int subject;
	string org;
	double score;
	double evalue;
	int lse;
	double curLseSize;
	double qStart;
	double qEnd;
	double sStart;
	double sEnd;
	map<string, int> bet4q; // org, subject: BeTs for current query
	map<string, double> bestScore; // org, best current score
	map<string, int>::iterator p_m;
	
	static ifstream inFile((inFileName + "hits.csv").c_str());
	if (!inFile)
	{
		cerr << "File hits.csv not found: can't create cogs" << endl;
		exit(1);
	}	
	
	set<int>::iterator p_s = (*curLse).begin();
	if ((*db.prot2org).find(*p_s) == (*db.prot2org).end())
	{
		cerr << "organism for " << (*p_s) << " not found" << endl;
		return;
	}

	organism = (*(*db.prot2org).find(*p_s)).second; // set organism for a given Lse

	while (p_s != (*curLse).end()) // loop for all proteins in given Lse
	{
		curQuery = *p_s;
		filter.setFilter(curQuery, db.offsetsF, inFileName);

		if ((*db.offsetsU).find(curQuery) != (*db.offsetsU).end())
		{
			offset = (*(*db.offsetsU).find(curQuery)).second; // find curQuery in hits.csv

			inFile.clear();
			inFile.seekg(offset);
			
			while (getline(inFile, str)) // loop for all hits for a given query protein
			{
#ifdef DEBUG1
				cout << str << endl;
#endif
				split(str, ",", &v);
				query = atoi(v[0].c_str());
				subject = atoi(v[1].c_str());
				if (query == curQuery && (*db.prot2org).find(subject) != (*db.prot2org).end() && (*db.prot2lse).find(subject) != (*db.prot2lse).end())
				{
					if (filter.isRelevant(subject) == true)
					{
						org = (*(*db.prot2org).find(subject)).second;
						score = atof(v[7].c_str());
						evalue = atof(v[6].c_str());
						if (org != organism && evalue <= thresholds.expect) // skip hits to self genome
						{
							qStart = atoi(v[2].c_str());
							qEnd = atoi(v[3].c_str());
							sStart = atoi(v[4].c_str());
							sEnd = atoi(v[5].c_str());

							if ((qEnd - qStart)/getLen(query) >= thresholds.overlap &&
								(sEnd - sStart)/getLen(subject) >= thresholds.overlap)
							{
								if (bet4q.find(org) == bet4q.end()) // first hit in this genome
								{
									bestScore.insert(pair<string, double>(org, score));
									bet4q.insert(pair<string, int>(org, subject));
								}
								else
								{
									if (score > (*bestScore.find(org)).second) // replace hit for better one
									{
										bestScore.erase(bestScore.find(org));
										bestScore.insert(pair<string, double>(org, score));

										bet4q.erase(bet4q.find(org));
										bet4q.insert(pair<string, int>(org, subject));
									}
								}
							}						
						}
					}
				}
				else
				{
					break;
				}
			}// while loop for all hits for a given query protein

#ifdef DEBUG1
			cout << "LSE " << *(*curLse).begin() << " has " << bet4q.size() << " best hits for protein " << *p_s << endl;
#endif

			p_m = bet4q.begin();
			while (p_m != bet4q.end())
			{
				subject = (*p_m).second;
				org = (*p_m).first;
				score = (*bestScore.find(org)).second;
				lse = (*(*db.prot2lse).find(subject)).second; // substitute subject protein with lse for this protein

				Hsp* hsp = new Hsp;

				hsp->score = score;
				hsp->subjectLse = lse;
				hsp->subjectOrg = org;

				lseHits.insert(pair<int, Hsp>(lse, *hsp));

				delete hsp;
				p_m++;
			}
			bestScore.clear();
			bet4q.clear();

		}// if offset found for a given protein

		p_s++;

	} // while loop for all proteins in given Lse

	curLseSize = (*curLse).size();

#ifdef DEBUG1
			cout << "LSE " << *(*curLse).begin() << " has " << lseHits.size() << " best hits for " << curLseSize << " proteins" << endl;
#endif

	multimap<int, Hsp>::iterator p_mm = lseHits.begin();
	multimap<string, Hsp>::iterator p_last;
	int lseInserted = -1;

	while (p_mm != lseHits.end()) // insert sides for future triangles
	{
		lse = (*p_mm).first; // lse for subject
		
		if (lseHits.count(lse) > curLseSize/2 && lse != lseInserted)
		{
			sides.insert(pair<int, int>(*(*curLse).begin(), lse));
			lseInserted = lse;
		}

		if (lseHits.count(lse) == curLseSize/2)
		{
			org = (*p_mm).second.subjectOrg;
			if (lse != lseInserted)
			{
				p_last = conflict.insert(pair<string, Hsp>(org, (*p_mm).second));
				lseInserted = lse;
			}
			else // not first hit into questionable LSE
			{
				score = (*p_mm).second.score;
				(*p_last).second.score += score; // compute total score for a questionable LSE
			}
		}
		p_mm++;

	} // while loop for inserting sides 

	multimap<string, Hsp>::iterator p_c = conflict.begin();
	while (p_c != conflict.end())
	{
		lse = (*p_c).second.subjectLse;
		
		if (conflict.count((*p_c).first) == 1) // only one questionable LSE for the organism
		{
			sides.insert(pair<int, int>(*(*curLse).begin(), lse));
		}
		else // should be exactly 2 btw
		{

			score = (*p_c).second.score;
			
			p_c++; // presuming the next is the second example of organism
			if (score < (*p_c).second.score)
			{
				lse = (*p_c).second.subjectLse; // change lse for better one
			}
			sides.insert(pair<int, int>(*(*curLse).begin(), lse));
		}
		p_c++;
	}
	inFile.clear();
	clear();

	map<string, int>::iterator mp = (*db.domHash).begin();
	while (mp != (*db.domHash).end())
	{
		keys.insert(pair<int, string>((*mp).second, (*mp).first));
		mp++;
	}
	(*db.domHash).clear();
}

string HitSet::getOrg(int prot)
{
	if ((*db.prot2org).find(prot) != (*db.prot2org).end())
	{
		return (*(*db.prot2org).find(prot)).second;
	}
	return "-1";
}

int HitSet::getLen(int prot)
{
	if ((*db.prot2len).find(prot) != (*db.prot2len).end())
	{
		return (*(*db.prot2len).find(prot)).second;
	}
	return -1;
}


#if 0
void HitSet::printRest(ofstream* fOut, set<int>* curLse)
{
	p_s = (*curLse).begin();
	while (p_s != (*curLse).end())
	{	
		if (processedLse.find(*p_s) == processedLse.end())
		{
			curQuery = (*p_s);
			
			*fOut << (*p_lse) << ',' << organism << ',';
				*fOut << (*p_lse) << ',' << len << ',' << '1' << ',' << len << ',';
				*fOut << cogName << setw(5) << setfill('0') << cogN << ',' << endl;
				(*fOut).clear(); 
		}
	}
}
#endif


typedef set< pair<int, int> > edges_type;
typedef vector< vector<int> > vertex_map_type;
typedef vector<int>::const_iterator v_int_cit;

// Edge-based search: As each new edge is added, search for the new
// triangles that extend the COG that might result.

static inline void try_insert_edge(int v1, int v2, edges_type &cog_edges,
								   queue< pair<int, int> > &edge_queue)
{
	const pair<int, int> edge(v1, v2), edge_rev(v2, v1);

	if (cog_edges.find(edge) == cog_edges.end())
	{
		cog_edges.insert(edge);
		cog_edges.insert(edge_rev);
		edge_queue.push(edge);
	}
}

static inline void check_edge(const pair<int, int> &edge,
							  const vertex_map_type &vertex_map,
							  const edges_type &symmetric_edges,
							  edges_type &cog_edges,
							  queue< pair<int, int> > &edge_queue)
{
	assert(cog_edges.find(edge) != cog_edges.end());

	// seed is edge endpoint having smaller degree (minor speedup)
	int seed = edge.first, nonseed = edge.second;
	if (vertex_map[edge.first].size() > vertex_map[edge.second].size())
		swap(seed, nonseed);

	for (v_int_cit p1=vertex_map[seed].begin();
		 p1 != vertex_map[seed].end(); p1++)
	{
		const int third = *p1;
		if (third == nonseed)
			continue;
		if (symmetric_edges.find(make_pair(third, nonseed))
			== symmetric_edges.end())
			continue;
		try_insert_edge(third, seed, cog_edges, edge_queue);
		try_insert_edge(third, nonseed, cog_edges, edge_queue);
	}
}

// An alternate (and perhaps better) strategy would be to interleave the
// second phase adds within the first phase, instead of explicitly switching
// between the two at coarse granularity, as is done here.

static inline void makeCOG(const int root_vertex, const int root_adj_vertex,
						   const vertex_map_type &vertex_map,
						   const edges_type &symmetric_edges,
						   edges_type &cog_edges, edges_type &processed)
{
	queue< pair<int, int> > edge_queue;
	edge_queue.push(make_pair(root_vertex, root_adj_vertex));

	// In the first phase, try to form triangles from each edge in
	// edge_queue, adding new triangle edges to the queue.
	while (not edge_queue.empty())
	  {
		const pair<int, int> edge = edge_queue.front();
		edge_queue.pop();

		processed.insert(edge);
		processed.insert(make_pair(edge.second, edge.first));

		check_edge(edge, vertex_map, symmetric_edges, cog_edges, edge_queue);
	  }
}

void HitSet::makeCOGs(ofstream* fOut)
{
	edges_type edges(sides.begin(), sides.end());

	// symmetric, non-self edges only
	edges_type symmetric_edges;

	const int max_vertex = db.prot2org->rbegin()->first;

	ofstream edgefile("all-edges.txt");

	for (edges_type::iterator ep=edges.begin(); ep != edges.end(); ep++)
	{
		int from = ep->first, to = ep->second;
		if (edges.find(make_pair(to, from)) != edges.end()
			and from < to)
		{
			// check that protein orgs are okay (too paranoid?)
			map<int, string>::const_iterator from_p = db.prot2org->find(from);
			map<int, string>::const_iterator to_p = db.prot2org->find(to);
			assert(from_p != db.prot2org->end());
			assert(to_p != db.prot2org->end());
			assert(from_p->second != to_p->second);

			symmetric_edges.insert(make_pair(from, to));
			symmetric_edges.insert(make_pair(to, from));

			string from_prot, to_prot;
			map<int, string>::iterator mp;
			mp = keys.find(from);
			if (mp != keys.end())
				from_prot = mp->second;
			mp = keys.find(to);
			if (mp != keys.end())
				to_prot = mp->second;

			edgefile << from << "," << from_prot << ","
					 << to << "," << to_prot << endl;
			edgefile << to << "," << to_prot << ","
					 << from << "," << from_prot << endl;

			assert(0 <= from and from <= max_vertex);
			assert(0 <= to and to <= max_vertex);
		}
	}
	edgefile.close();
	edges.clear();				// reclaim memory

	// vertex_map: map a vertex to a vector of adjacent vertices
	//             (same graph as symmetric_edges)
	vertex_map_type vertex_map(max_vertex+1);

	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
		vertex_map[ep->first].push_back(ep->second);

	// true iff an edge is already part of a COG or TWOG
	// (only contains ascending half-edge; v1 -> v2 where v1 < v2)
	edges_type processed;
	edges_type triprocessed;

	int cogN = cogStartNum;		// number assigned to cog in output

	ofstream cogfile("cog-edges.txt");

	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
	{
		const int root_vertex = ep->first, root_adj_vertex = ep->second;
		if (not (root_vertex < root_adj_vertex)
			or (processed.find(*ep) != processed.end()))
			continue;

		edges_type cog_edges;
		// root_vertex--root_adj_vertex are endpoints of initial COG edge
		cog_edges.insert(make_pair(root_vertex, root_adj_vertex));
		cog_edges.insert(make_pair(root_adj_vertex, root_vertex));

		// now expand from that initial edge
		makeCOG(root_vertex, root_adj_vertex, vertex_map, symmetric_edges,
				cog_edges, processed);

		// if we found a COG, print it
		set<int> cog_vertices;	// vertices in cog_edges
		for (edges_type::const_iterator ep=cog_edges.begin();
			 ep != cog_edges.end(); ep++)
		  cog_vertices.insert(ep->first);

		if (cog_edges.size() >= 3)
		{
			for (edges_type::const_iterator ep=cog_edges.begin();
				 ep != cog_edges.end(); ep++)
			{
				triprocessed.insert(make_pair(ep->first, ep->second));
				triprocessed.insert(make_pair(ep->second, ep->first));
				int from = ep->first;
				int to = ep->second;
				map<int, string>::const_iterator from_p = keys.find(from);
				map<int, string>::const_iterator to_p = keys.find(to);

				cogfile << cogName << setw(5) << setfill('0') << cogN << ","
						<< from_p->second << "," << getOrg(from) << ","
						<< to_p->second << "," << getOrg(to) << endl;
			}
			printCOG(cogN++, &cog_vertices, fOut);
		}
	}

	cogfile.close();

#ifndef NO_TWOG_SINGLETON
	// print TWOGs (unfortunately a hack)
	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
	{
		const int root_vertex = ep->first, root_adj_vertex = ep->second;
		if (not (root_vertex < root_adj_vertex)
			or (triprocessed.find(*ep) != triprocessed.end()))
			continue;

		edges_type cog_edges;
		// root_vertex--root_adj_vertex are endpoints of initial COG edge
		cog_edges.insert(make_pair(root_vertex, root_adj_vertex));
		cog_edges.insert(make_pair(root_adj_vertex, root_vertex));

		// if we found a COG, print it
		set<int> cog_vertices;	// vertices in cog_edges
		for (edges_type::const_iterator ep=cog_edges.begin();
			 ep != cog_edges.end(); ep++)
		  cog_vertices.insert(ep->first);

		printCOG(cogN++, &cog_vertices, fOut);
	}

	set<int> processed_vertices;
	for (edges_type::const_iterator p=processed.begin();
		 p != processed.end(); p++)
	{
		processed_vertices.insert(p->first);
		processed_vertices.insert(p->second);
	}

	for (map<int, Lse>::const_iterator p_lse = db.allLse->begin();
		 p_lse != db.allLse->end(); p_lse++)
	{
		const int vertex = p_lse->first;
		if (processed_vertices.find(vertex) == processed_vertices.end())
		{
			set<int> cog_vertices;
			cog_vertices.insert(vertex);
			printCOG(0, &cog_vertices, fOut);
		}
	}
#endif
}

void HitSet::printCOG(int cogN, set<int>* cog, ofstream* fOut)
{
	string prot;
	set<int>::iterator p_lse;
	map<int, string>::iterator mp;

	string cogSuffix = "";
//	if (cog->size() == 2)
//		cogSuffix = "T";
//	else if (cog->size() == 1)
//		cogSuffix = "S";

	for (set<int>::const_iterator p_cog = cog->begin();
		 p_cog != cog->end(); p_cog++)
	{
		const int lse = *p_cog;

		const set<int> lses = db.allLse->find(lse)->second.lse;
		for (set<int>::const_iterator p_lse = lses.begin();
			 p_lse != lses.end(); p_lse++)
		{
			mp = keys.find(*p_lse);
			if (mp != keys.end())
			{
				prot = mp->second;
				*fOut << prot << ',' << getOrg(*p_lse) << ',';
				*fOut << prot << ',' << getLen(*p_lse) << ',' << '1' << ',' << getLen(*p_lse) << ',';
				if (cogN == 0)
				{
					*fOut << ',' << endl;
				}
				else
				{
					*fOut << cogName << cogSuffix << setw(5) << setfill('0') << cogN << ',' << endl;
				}
			}
		}
	}
}

void makeOffsets(map<int, long>* offsets, ifstream& inFile)
{
	long offset;
	string str;
	string curProt;
	string prot;
	size_t indx;
	
	offset = inFile.tellg();
	while (getline(inFile, str))
	{
		indx = str.find_first_of(",");
		prot = str.substr(0, indx);
		if (prot != "")
		{
			if (curProt != prot)
			{
				(*offsets).insert(pair<int, long>(atoi(prot.c_str()), offset));
				curProt = prot;
			}
		}
		offset = inFile.tellg();
	}
}

void populate (map<string, int>* domHash, map<int, string>* p2ent, ifstream& inFile)
{
	vector<string> v;
	string str;
	string protStr;
	map<string, int>::iterator mp;
	int prot;	
	string ent;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
        protStr = v[0];
		mp = (*domHash).find(protStr);
		if (mp != (*domHash).end())
		{
			prot = (*mp).second;
			ent = v[1];
			(*p2ent).insert(pair<int, string>(prot, ent));
		}
	}
}

void populate(map<int, string>* p2ent, ifstream& inFile)
{
	vector<string> v;
	string str;
	int prot;	
	string ent;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		prot = atoi(v[0].c_str());
		ent = v[1];
		
		(*p2ent).insert(pair<int, string>(prot, ent));
	}
}

void populate(map<string, string>* p2ent, ifstream& inFile)
{
	vector<string> v;
	string str;
	string prot;	
	string ent;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		prot = v[0];
		ent = v[1];
		
		(*p2ent).insert(pair<string, string>(prot, ent));
	}
}

void populate(map<int, int>* p2len, ifstream& inFile)
{
	vector<string> v;
	string str;
	int prot;	
	int len;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		prot = atoi(v[0].c_str());
		len = atoi(v[1].c_str());
		
		(*p2len).insert(pair<int, int>(prot, len));
	}
}

void populate (map<string, int>* domHash, map<int, int>* prot2lse, map<int, Lse>* allLse, ifstream& inFile)
{
	vector<string> v;
	string str;
	string protStr;
	map<string, int>::iterator mp;
	int prot;
	int firstProt; // first protein in LSE uses as an ID for its LSE
	unsigned int i;

	while(getline(inFile, str))
	{
		if (str != "")
		{
			Lse* curLse = new Lse;
			split(str, ",", &v);
			
			for (i = 0; i < v.size(); i++)
			{
				protStr = v[i];
				mp = (*domHash).find(protStr);
				
				if(mp != (*domHash).end())
				{
					prot = (*mp).second;
					(*curLse).lse.insert(prot);
				}
				else
					continue;
			}

			firstProt = (*(*curLse).lse.begin());

			set<int>::iterator p_s = (*curLse).lse.begin();
			while (p_s != (*curLse).lse.end())
			{
				(*prot2lse).insert(pair<int, int>((*p_s), firstProt));
				p_s++;
			}

			(*allLse).insert(pair<int, Lse>(firstProt, (*curLse)));

			delete curLse;
		}
	}
}

int split(const string str, const string sepStr, vector<string>* v)
{
    string	token;
    int		nToken = 0;
    size_t	index;
    size_t	posSep;
    size_t	strLen;
    size_t	tokenLen;

	(*v).clear();

    strLen = str.length();

    index = str.find_first_not_of(sepStr);// skip leading separators

    do
    {
	posSep = str.find_first_of(sepStr, index);
	if (posSep == string::npos)
	{// separator not found
	    tokenLen = strLen - index;
	}
	else
	{
	    tokenLen = posSep - index;
	}
	token = str.substr(index, tokenLen);
	(*v).push_back(token);
	nToken++;
//	index = str.find_first_not_of(sepStr, posSep + 1);//skip all contin. separators
	index = posSep + 1;
    }
	while(posSep != string::npos && index != string::npos);
    
    return nToken;
}
