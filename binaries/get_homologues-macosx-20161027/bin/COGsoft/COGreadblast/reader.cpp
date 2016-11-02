#include "reader.h"

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

FHitEnt::FHitEnt():query(0), target(0)
{}

FHitEnt::FHitEnt(const long Id1 , const long Id2)
{
	query = Id1; 
	target = Id2;
}

FHitEnt::FHitEnt(const FHitEnt& ent)
{
	query = ent.query;
	target = ent.target;
}

SelfHitEnt::SelfHitEnt():id(0), len(0), score(0.0)
{}

SelfHitEnt::SelfHitEnt(const long i, const int l, const double s)
{
	id = i; 
	len = l;
	score = s;
}

SelfHitEnt::SelfHitEnt(const SelfHitEnt& ent)
{
	id = ent.id;
	len = ent.len;
	score = ent.score;
}

HitEnt::HitEnt():query(0), target(0), qStart(0), qEnd(0), tStart(0), tEnd(0), evalue(0.0), score(0.0)
{}

HitEnt::HitEnt(const long q, const long t, const int qS, const int qE, const int tS, const int tE, const double e, const double s)
{
	query = q;
	target = t; 
	qStart = qS;
	qEnd = qE;
	tStart = tS;
	tEnd = tE;
	evalue = e;
	score = s;
}

HitEnt::HitEnt(const HitEnt& ent)
{
	query = ent.query;
	target = ent.target; 
	qStart = ent.qStart;
	qEnd = ent.qEnd;
	tStart = ent.tStart;
	tEnd = ent.tEnd;
	evalue = ent.evalue;
	score = ent.score;
}

Reader::Reader(const char* fileName):indx1(0), indx2(0), size(0)
{
	char c = 0;
	char beforec = 0;

	inFile.open(fileName);
	if (!inFile)
	{
		cerr << "Can't open input file " << fileName << endl;
		exit(1);
	}
	
	while(c != '\n') // analize first line
	{	
		beforec = c;
		c = inFile.get();

		if (c == ',')
		{
			++size;
		}
	}

	if (size > 0 && beforec != ',')
	{
		++size; // last token
	}
	
	close();
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

FHitReader::FHitReader(const char* fileName)
: Reader(fileName)
{}

bool FHitReader::next(FHitEnt& fe)
{
	if(getline(inFile, str))
	{
		indx1 = str.find_first_of(',');
		fe.query = atoi((str.substr(0, indx1)).c_str());
		fe.target = atoi((str.substr(indx1 + 1)).c_str());
		return true;
	}
	else return false;
}


SelfHitReader::SelfHitReader(const char* fileName)
: Reader(fileName)
{}

bool SelfHitReader::next(SelfHitEnt& se)
{
	if(getline(inFile, str))
	{
		indx1 = str.find_first_of(',');
		se.id = atoi((str.substr(0, indx1)).c_str());
		indx2 = str.find_first_of(',', indx1 + 1);
		se.len = atoi((str.substr(indx1 + 1, indx2)).c_str());
		se.score = atof((str.substr(indx2 + 1)).c_str());
		return true;
	}
	else return false;
}

HitReader::HitReader(const char* fileName)
: Reader(fileName)
{}

bool HitReader::next(HitEnt& h)
{
	if(getline(inFile, str))
	{
		indx1 = str.find_first_of(',');
		h.query = atoi((str.substr(0, indx1)).c_str());
		
		indx2 = str.find_first_of(',', indx1 + 1);
		h.target = atoi((str.substr(indx1 + 1, indx2)).c_str());
		
		indx1 = indx2;
		indx2 = str.find_first_of(',', indx1 + 1);
		h.qStart = atoi((str.substr(indx1 + 1, indx2)).c_str());
		
		indx1 = indx2;
		indx2 = str.find_first_of(',', indx1 + 1);
		h.qEnd = atoi((str.substr(indx1 + 1, indx2)).c_str());
		
		indx1 = indx2;
		indx2 = str.find_first_of(',', indx1 + 1);
		h.tStart = atoi((str.substr(indx1 + 1, indx2)).c_str());

		indx1 = indx2;
		indx2 = str.find_first_of(',', indx1 + 1);
		h.tEnd = atoi((str.substr(indx1 + 1, indx2)).c_str());

		indx1 = indx2;
		indx2 = str.find_first_of(',', indx1 + 1);
		h.evalue = atof((str.substr(indx1 + 1, indx2)).c_str());

		h.score = atof((str.substr(indx2 + 1)).c_str());

		return true;
	}
	else return false;
}
