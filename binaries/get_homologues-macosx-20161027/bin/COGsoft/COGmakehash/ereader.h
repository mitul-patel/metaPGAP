#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#ifndef EREADER_H
#define EREADER_H

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

class Reader
{
protected:
	ifstream inFile;
	string str;
	size_t indx1;
	size_t indx2;
//	int size; // number of tokens in the (first) string
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


class EasyReader : public Reader
{
public:
	EasyReader(const char* fileName, const char sprtr, const int block1);
	~EasyReader();
	bool next(KeyEnt&);
private:
	int counter;
	char sprtr;
	unsigned int nameBlock;
    string	token;
    int		nToken;
    size_t	strLen;
    size_t	tokenLen;
	
	vector<string> v;

	int split(const string str);
};

#endif

