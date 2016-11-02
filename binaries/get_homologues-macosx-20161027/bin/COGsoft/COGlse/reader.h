#ifndef _READER_H
#define _READER_H


#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

struct Entry
{
	virtual void init(const string& str, int size) = 0;
	virtual ostream& write(ostream& os) = 0;
	virtual ~Entry(){};
};

ostream& operator << (ostream& os, Entry& e);


struct KeyEnt : public Entry
{	
	long id;
	string name;
	string org;
	string cog;

	KeyEnt();
	KeyEnt(const long, const string, const string, const string);
	KeyEnt(const KeyEnt&);

	
	void init(const string& str, int size);
	ostream& write(ostream& os)
	{
		return os << id
			<< "," << name 
			<< "," << org;
	};
};

// filtered hit
struct FHitEnt  : public Entry
{
	long query;
	long target;

	FHitEnt();
	FHitEnt(const long, const long);
	FHitEnt(const FHitEnt&);
	void init(const string& str, int size);

	bool operator<(const FHitEnt& f) const
    {
       if (query < f.query)
	   {
		   return true;
	   }
	   return false;			   
    }

	ostream& write(ostream& os)
	{
		return os << query
			<< "," << target;
	}

};

struct SelfHitEnt : public Entry
{
	long id;
	int len;
	double score;

	SelfHitEnt();
	SelfHitEnt(const long, const int, double score);
	SelfHitEnt(const SelfHitEnt&);

	void init(const string& str, int size);
	ostream& write(ostream& os)
	{
		return os << id
			<< "," << len
			<< "," << score;
	}

};


struct HitEnt : public Entry
{
	long query;
	long target;
	int qStart;
	int qEnd;
	int tStart;
	int tEnd;
	double evalue;
	double score;

	HitEnt();
	HitEnt(const long, const long, const int, const int, const int, const int, double, double);
	HitEnt(const HitEnt&);

	void init(const string& str, int size);

	bool operator<(const HitEnt& h) const
    {
       if (query < h.query)
	   {
		   return true;
	   }
	   else if (query == h.query)
	   {
		   if (score > h.score)
		   {
			   return true;
		   }
	   }
	   return false;			   
    }

	ostream& write(ostream& os)
	{
		return os << query
			<< "," << target
			<< "," << qStart
			<< "," << qEnd
			<< "," << tStart
			<< "," << tEnd
			<< "," << evalue
			<< "," << score;
	}

};

struct PairEnt : public Entry
{	
	string code1;
	string code2;

	PairEnt() : Entry(), code1(""), code2(""){};
	PairEnt(const string& code1, const string& code2) : Entry()
	{
		this->code1 = code1;
		this->code2 = code2;
	};

	PairEnt(const PairEnt& e)
	{
		code1 = e.code1;
		code2 = e.code2;
	}

	void init(const string& str, int size)
	{
		size_t indx1 = str.find_first_of(',');
		code1 = str.substr(0, indx1);
		
		indx1++;
		size_t indx2 = str.find_first_of(',', indx1);
		code2 = str.substr(indx1, indx2 - indx1);
	}
	ostream& write(ostream& os)
	{
		return os << code1
			<< "," << code2;
	}
};

struct OrgEnt : public Entry
{	
	string protCode;
	string orgCode;

	OrgEnt() : Entry(), protCode(""), orgCode(""){};
	OrgEnt(const string& protCode, const string& orgCode) : Entry()
	{
		this->protCode = protCode;
		this->orgCode = orgCode;
	};

	OrgEnt(const OrgEnt& e)
	{
		protCode = e.protCode;
		orgCode = e.orgCode;
	}

	void init(const string& str, int size)
	{
		size_t indx1 = str.find_first_of(',');
		protCode = str.substr(0, indx1);
		
		indx1++;
		size_t indx2 = str.find_first_of(',', indx1);
		orgCode = str.substr(indx1, indx2 - indx1);
	}
	ostream& write(ostream& os)
	{
		return os << protCode
			<< "," << orgCode;
	}
};


class LineEntryReader
{
protected:
	ifstream inFile;
	string str;
	size_t indx1;
	size_t indx2;
	int size; // number of tokens in the (first) string
	ifstream::pos_type curFilePos;
	ifstream::pos_type endFilePos;

public:
	LineEntryReader(const char* fileName);
	~LineEntryReader();
	bool next(Entry& entry);
	void close();
	int readPercent();
};

#endif
