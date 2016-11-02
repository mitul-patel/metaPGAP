// COGMAKER.H
#ifndef COGMAKER_H
#define COGMAKER_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cstdlib>

using namespace std;

//#define NO_TWOG_SINGLETON
//#define DEBUG

#ifdef DEBUG
	#include <stdio.h>
	#include <time.h>
struct Stopwatch
{
	clock_t start, end;
	void startClock(){ start = clock(); }
	void finishClock()
	{
		end = clock();
		printf("clocks: %f\n",((double)(end-start))/CLOCKS_PER_SEC);
		cout << (double)(end - start) << endl;
	}
};
#endif // DEBUG

struct Thresholds
{
	double	overlap; // minimum overlap between 2 sectors of the same COG to be joined
	double	expect; // e-value cut-off
};

struct Lse
{
	set<int> lse;
};

struct DB
{
	map<string, int>*	domHash;
	map<int, string>*	prot2org;
	map<int, int>*		prot2len;
	map<int, long>*		offsetsF; // for filtered hits
	map<int, long>*		offsetsU; // for unfiltered hits
	map<int, Lse>*		allLse; // every protein has lse, at least for itself (lse.size() = 1)
	map<int, int>*		prot2lse;
};

struct Side
{
	int subjct;
	float score;
};


struct Hsp
{
	int		subjectLse;
	string	subjectOrg;
	double	score;
};

class Filter
{
public:
	void	setFilter	(const int query, map<int, long>* offset, const string inFileName);
	bool	isRelevant	(const int subject);
private:
	set<int> mask; // filtered hits for a given lse
};

class HitSet // includes all or some HSPs for a given lse
{
public:
			HitSet	(const double	expect,
					const double	overlap,
					map<string, int>* domHash,
					map<int, string>*	prot2org,
					map<int, int>* prot2len,
					map<int, long>*	offsetsF,
					map<int, long>*	offsetsU,
					map<int, Lse>*	allLse,
					map<int, int>*	prot2lse,
					const string	cName,
					const int		startNum,
					const string	inFile);

			void insert	(set<int>* curLse);
			void makeCOGs(ofstream* fOut);
			void printRest(ofstream* fOut, set<int>* curLse);
			void clear	();	
	
private:
	string				cogName;
	int					cogStartNum;
	string				inFileName;
	Thresholds			thresholds;
	DB					db;
	string					organism;	// query (LSE) org
	Filter				filter;
	map<int, string> keys;
	multimap<int, Hsp>	lseHits; // for current LSE only: subjectLse, Hsp
	multimap<string, Hsp>	conflict; // subjectOrg, Hsp
	multimap<int, int>	sides; // apexQ, apexT
	set<int> processedLse;
	string getOrg(int prot);
	int getLen(int prot);
	void printCOG(int cogN, set<int>* cog, ofstream* fOut);
};

int		split		(const string str, const string sepStr, vector<string>* v);
void	makeOffsets	(map<int, long>* offsets, ifstream& inFile);
void	populate	(map<int, string>* p2ent, ifstream& inFile);
void	populate	(map<string, int>* domHash, map<int, string>* p2ent, ifstream& inFile);
void	populate	(map<int, int>* p2len, ifstream& inFile);
void	populate	(map<string, string>* p2ent, ifstream& inFile);
void	populate	(map<string, int>* domHash, map<int, int>* prot2lse, map<int, Lse>* allLse, ifstream& inFile);

#endif
