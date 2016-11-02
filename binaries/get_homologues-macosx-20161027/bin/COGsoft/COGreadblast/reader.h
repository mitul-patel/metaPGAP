#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

struct KeyEnt
{	
	long id;
	string name;

	KeyEnt();
	KeyEnt(const long, const string);
	KeyEnt(const KeyEnt&);
};

struct FHitEnt // filtered hit
{
	long query;
	long target;

	FHitEnt();
	FHitEnt(const long, const long);
	FHitEnt(const FHitEnt&);
	bool operator<(const FHitEnt& f) const
    {
       if (query < f.query)
	   {
		   return true;
	   }
	   return false;			   
    }
};

struct SelfHitEnt
{
	long id;
	int len;
	double score;

	SelfHitEnt();
	SelfHitEnt(const long, const int, double score);
	SelfHitEnt(const SelfHitEnt&);
};

struct HitEnt
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
};

class Reader
{
protected:
	ifstream inFile;
	string str;
	size_t indx1;
	size_t indx2;
	int size; // number of tokens in the (first) string
public:
	Reader(const char* fileName);
	~Reader();
	void close();
};

class KeyReader : public Reader
{
public:
	KeyReader(const char* fileName);
	bool next(KeyEnt&);
};

class FHitReader : public Reader
{
public:
	FHitReader(const char* fileName);
	bool next(FHitEnt&);
};

class SelfHitReader : public Reader
{
public:
	SelfHitReader(const char* fileName);
	bool next(SelfHitEnt&);
};

class HitReader : public Reader
{
public:
	HitReader(const char* fileName);
	bool next(HitEnt&);
};

#endif
