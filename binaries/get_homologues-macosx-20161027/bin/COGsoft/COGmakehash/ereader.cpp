#include "ereader.h"

KeyEnt::KeyEnt():id(0), name("")
{}

KeyEnt::KeyEnt(const long entId, const string entName)
{
	id = entId; 
	name = entName;
}

KeyEnt::KeyEnt(const KeyEnt& ent)
{
	id = ent.id;
	name = ent.name;
}

Reader::Reader(const char* fileName):indx1(0), indx2(0)//, size(0)
{
	inFile.open(fileName);
	if (!inFile)
	{
		cerr << "Can't open input file " << fileName << endl;
		exit(1);
	}
}

void Reader::close()
{
	inFile.close();
}

Reader::~Reader()
{
	close();
}

KeyReader::KeyReader(const char* fileName)
: Reader(fileName)
{}

bool KeyReader::next(KeyEnt& ke)
{
	if(getline(inFile, str))
	{
		indx1 = str.find_first_of(',');
		ke.id = atoi((str.substr(0, indx1)).c_str());
		ke.name = str.substr(indx1 + 1);
		
		return true;
	}
	else return false;
}


EasyReader::EasyReader(const char* fileName, const char sep, const int block1)
: Reader(fileName)
{
	sprtr = sep;
	nameBlock = block1 - 1;
	if (nameBlock < 0)
	{
		cerr << "Incorrect index of the name field: " << block1 << endl;
		exit(1);
	}
	nToken = 0;
	counter = 0;

}

EasyReader::~EasyReader()
{
	v.clear();
}

int EasyReader::split(const string str)
{
	v.clear();

    strLen = str.length();

    indx1 = str.find_first_not_of(sprtr);// skip leading separators

    do
    {
	indx2 = str.find_first_of(sprtr, indx1);
	if (indx2 == string::npos)
	{// separator not found
	    tokenLen = strLen - indx1;
	}
	else
	{
	    tokenLen = indx2 - indx1;
	}
	token = str.substr(indx1, tokenLen);
	v.push_back(token);
	nToken++;
//	indx1 = str.find_first_not_of(sprtr, indx2 + 1);//skip all contin. separators
	indx1 = indx2 + 1;
    }
	while(indx2 != string::npos && indx1 != string::npos);
    
    return nToken;
}

bool EasyReader::next(KeyEnt& ke)
{
	if(getline(inFile, str))
	{
		if (str.find('#') != str.npos) 
		{
			while (getline(inFile, str))
			{
				if (str.find('#') == str.npos)
					break;
			}
		}

		split(str);
		if (v.size() > nameBlock) 
		{
			if (v[nameBlock] != "")
				ke.name = v[nameBlock];
		}	
		else
		{
			cerr << "String " << str << " does not contain " << nameBlock << " blocks " << endl;
			exit(1);
		}
		return true;
	}
	return false;
}
