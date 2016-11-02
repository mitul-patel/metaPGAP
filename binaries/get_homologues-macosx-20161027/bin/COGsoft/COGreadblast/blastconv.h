#ifndef BLASTCONV_H
#define BLASTCONV_H

#include <list>
#include "blastconvglob.h"

struct Hit
{
	long q64;
	long t64;
	long qStart;
	long qEnd;
	long tStart;
	long tEnd;
	double evalue;
	double score;

	Hit(long, long, long, long, long, long, double, double);
	Hit();

	bool operator== (const Hit& h) const
	{
		return (q64 == h.q64 && t64 == h.t64 && qStart == h.qStart && qEnd == h.qEnd && tStart == h.tStart && tEnd == h.tEnd);
	}
};

struct QnS
{
	long q64;
	long t64;

	QnS(long, long);
	QnS();

	bool operator== (const QnS& qns) const
	{
		return (q64 == qns.q64 && t64 == qns.t64);
	}
};

struct cmpScore
{
	bool operator() (const Hit& h1, const Hit& h2) const
	{
		return (h1.score > h2.score);
	}
};

struct cmpT
{
	bool operator() (const Hit& h1, const Hit& h2) const
	{
		if (h1.t64 == h2.t64)
		{
			if (h1.qStart == h2.qStart)
			{
				if (h1.qEnd == h2.qEnd)
				{
					if (h1.tStart == h2.tStart)
					{
						if (h1.tEnd == h2.tEnd)
						{
							return (h1.score > h2.score);
						}
						return (h1.tEnd < h2.tEnd);
					}
					return (h1.tStart < h2.tStart);
				}
				return (h1.qEnd < h2.qEnd);
			}
			return (h1.qStart < h2.qStart);
		}
		return (h1.t64 < h2.t64);
	}
};

struct cmpT1
{
	bool operator() (const QnS& qns1, const QnS& qns2) const
	{
		return (qns1.t64 < qns2.t64);
	}
};

class HitSet
{
private:
	list<Hit> hits4q;
public:
	void insert(Hit);
	void correct();
	void printAll(fstream*);
	void print(fstream*);
	void sortOnScore();
	void sortOnT();
	void clear();
};

class QSet
{
private:
	list<QnS> s4q;
public:
	void insert(QnS);
	void correct();
	void printAll(fstream*);
	void print(fstream*);
	void sortOnT();
	void clear();
};

class Offfile
{
private:
	multimap<int, long> offsets;
	multimap<int, long>::iterator p;
	bool reverse;
	string curProt;
	string prot;
	int proti;
	int curProti;
	bool endof;
	ifstream inFile;
	string fn;
	void makeOffsets();
public:
	Offfile(const bool rev, const char* fileName);
	string nextLine();
	void add();
};
#endif
