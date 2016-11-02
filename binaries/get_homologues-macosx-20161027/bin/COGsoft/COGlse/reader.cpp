#include "reader.h"


ostream& operator << (ostream& os, Entry& e)
{
	return e.write(os);
}

KeyEnt::KeyEnt(): Entry(), id(0), name(""), org(""), cog("")
{}

KeyEnt::KeyEnt(const long entId, const string entName, const string entOrg, const string entCog) : Entry()
{
	id = entId; 
	name = entName;
	org = entOrg;
	cog = entCog;
}

KeyEnt::KeyEnt(const KeyEnt& ent) : Entry()
{
	id = ent.id;
	name = ent.name;
	org = ent.org;
	cog = ent.cog;
}


void KeyEnt::init(const string& str, int size)
{
	size_t indx1 = str.find_first_of(',');
	id = atoi((str.substr(0, indx1)).c_str());

	indx1++;
	if (size == 3)
	{
		size_t indx2 = str.find_first_of(',', indx1);
		name = str.substr(indx1, indx2 - indx1);

		indx2++;
		org = str.substr(indx2);
	}

	if (size == 2) // no data about organisms
	{
		name = str.substr(indx1);
	}
}


FHitEnt::FHitEnt(): Entry(), query(0), target(0)
{}

FHitEnt::FHitEnt(const long Id1 , const long Id2) : Entry()
{
	query = Id1; 
	target = Id2;
}

FHitEnt::FHitEnt(const FHitEnt& ent) : Entry()
{
	query = ent.query;
	target = ent.target;
}

void FHitEnt::init(const string& str, int size)
{
	size_t indx1 = str.find_first_of(',');
	query = atoi((str.substr(0, indx1)).c_str());

	indx1++;
	target = atoi((str.substr(indx1)).c_str());
}



SelfHitEnt::SelfHitEnt(): Entry(),id(0), len(0), score(0.0)
{}

SelfHitEnt::SelfHitEnt(const long i, const int l, const double s): Entry()
{
	id = i; 
	len = l;
	score = s;
}

SelfHitEnt::SelfHitEnt(const SelfHitEnt& ent): Entry()
{
	id = ent.id;
	len = ent.len;
	score = ent.score;
}

void SelfHitEnt::init(const string& str, int size)
{
	size_t indx1 = str.find_first_of(',');
	id = atoi((str.substr(0, indx1)).c_str());

	indx1++;
	size_t indx2 = str.find_first_of(',', indx1);
	len = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx2++;
	score = atof((str.substr(indx2)).c_str());
}



HitEnt::HitEnt(): Entry(),query(0), target(0), qStart(0), qEnd(0), tStart(0), tEnd(0), evalue(0.0), score(0.0)
{}

HitEnt::HitEnt(const long q, const long t, const int qS, const int qE, const int tS, const int tE, const double e, const double s): Entry()
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

HitEnt::HitEnt(const HitEnt& ent): Entry()
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

void HitEnt::init(const string& str, int size)
{
	size_t indx1 = str.find_first_of(',');
	query = atoi((str.substr(0, indx1)).c_str());
	
	indx1++;
	size_t indx2 = str.find_first_of(',', indx1);
	target = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx1 = indx2 + 1;
	indx2 = str.find_first_of(',', indx1);
	qStart = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx1 = indx2 + 1;
	indx2 = str.find_first_of(',', indx1);
	qEnd = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx1 = indx2 + 1;
	indx2 = str.find_first_of(',', indx1);
	tStart = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx1 = indx2 + 1;
	indx2 = str.find_first_of(',', indx1);
	tEnd = atoi((str.substr(indx1, indx2 - indx1)).c_str());

	indx1 = indx2 + 1;
	indx2 = str.find_first_of(',', indx1);
	evalue = atof((str.substr(indx1, indx2 - indx1)).c_str());
	
	indx1 = indx2 + 1;
	score = atof((str.substr(indx1)).c_str());
}


LineEntryReader::LineEntryReader(const char* fileName):indx1(0), indx2(0), size(0)
{
	char c = 0;

	inFile.open(fileName);
	if (!inFile)
	{
		cerr << "Can't open input file " << fileName << endl;
		exit(1);
	}
	curFilePos = inFile.tellg();
	inFile.seekg( 0, ios_base::end );
	endFilePos = inFile.tellg();
	inFile.seekg( 0, ios_base::beg);

	while(c != '\n' ) // analize first line
	{
		c = inFile.get();
		if( c == char_traits<char>::eof()) break;
		if (c == ',')
		{
			++size;
		}
	}

	if (size > 0)
	{
		++size; // last token
	}
	
	inFile.seekg( 0, ios_base::beg);
}

void LineEntryReader::close()
{
	inFile.close();
}

LineEntryReader::~LineEntryReader()
{
	close();
}

bool LineEntryReader::next(Entry& entry)
{
	if(getline(inFile, str))
	{
		curFilePos = inFile.tellg();
		entry.init(str, size);
		return true;
	}
	else return false;
}

int LineEntryReader::readPercent()
{
	streamoff curPos = curFilePos;
	streamoff endPos = endFilePos;

	return (endPos == 0) ? 0 : (int) (curPos *100.0/endPos);
}

