#include "blastconv.h"

Hit::Hit(long q, long t, long qs, long qe, long ts, long te, double e, double s)
{
	q64 = q;
	t64 = t;
	qStart = qs;
	qEnd = qe;
	tStart = ts;
	tEnd = te;	
	evalue = e;
	score = s;
}

Hit::Hit()
{
	q64 = 0;
	t64 = 0;
	qStart = 0;
	qEnd = 0;
	tStart = 0;
	tEnd = 0;	
	evalue = 0.0;
	score = 0.0;
}

void HitSet::insert(Hit h)
{
	hits4q.push_back(h);
}

void HitSet::correct() 
{
	hits4q.unique();
}

void HitSet::printAll(fstream* f)
{
	list<Hit>::iterator pl = hits4q.begin();
	while (pl != hits4q.end())
	{
		*f << (*pl).q64 << ',' << (*pl).t64 << ',' << (*pl).qStart << ',' << (*pl).qEnd << ',' << (*pl).tStart << ',' << (*pl).tEnd << ',' << (*pl).evalue << ',' << (*pl).score << endl;
		if ((*pl).q64 != (*pl).t64)
		{
			*f << (*pl).t64 << ',' << (*pl).q64 << ',' << (*pl).tStart << ',' << (*pl).tEnd << ',' << (*pl).qStart << ',' << (*pl).qEnd << ',' << (*pl).evalue << ',' << (*pl).score << endl;
		}
		pl++;
	}
}

void HitSet::print(fstream* f)
{
	list<Hit>::iterator pl = hits4q.begin();
	while (pl != hits4q.end())
	{
		*f << (*pl).q64 << ',' << (*pl).t64 << ',' << (*pl).qStart << ',' << (*pl).qEnd << ',' << (*pl).tStart << ',' << (*pl).tEnd << ',' << (*pl).evalue << ',' << (*pl).score << endl;
		pl++;
	}
}

void HitSet::sortOnScore()
{
	hits4q.sort(cmpScore());
}

void HitSet::sortOnT()
{
	hits4q.sort(cmpT());
}

void HitSet::clear()
{
	hits4q.clear();
}

QnS::QnS()
{
	q64 = 0;
	t64 = 0;
}

QnS::QnS(long q, long t)
{
	q64 = q;
	t64 = t;
}

void QSet::insert(QnS qns)
{
	s4q.push_back(qns);
}

void QSet::correct() 
{
	s4q.unique();
}

void QSet::printAll(fstream* f)
{
	list<QnS>::iterator pl = s4q.begin();
	while (pl != s4q.end())
	{
		*f << (*pl).q64 << ',' << (*pl).t64  << endl;
		if ((*pl).q64 != (*pl).t64)
		{
			*f << (*pl).t64 << ',' << (*pl).q64  << endl;
		}
		pl++;
	}
}

void QSet::print(fstream* f)
{
	list<QnS>::iterator pl = s4q.begin();
	while (pl != s4q.end())
	{
		*f << (*pl).q64 << ',' << (*pl).t64  << endl;
		pl++;
	}
}

void QSet::sortOnT()
{
	s4q.sort(cmpT1());
}

void QSet::clear()
{
	s4q.clear();
}

Offfile::Offfile(const bool rev, const char* fileName)
{
	fn = fileName;
	reverse = rev;
	endof = false;
	proti = 0;
	curProti = 0;
	inFile.open(fileName);
	if (!inFile)
	{
		cerr << "can't find input file " << fileName << endl;
		exit(1);
	}
	makeOffsets();
	inFile.clear();
}
		

void Offfile::makeOffsets()
{
	long offset;
	string str;
	size_t indx1;

	offsets.clear();
	
	offset = inFile.tellg();
	while (getline(inFile, str))
	{
		indx1 = str.find_first_of(',');		
		prot = str.substr(0, indx1);

		if (curProt != prot) // keep offset for the first entry in group only
		{
			offsets.insert(pair<int, long>(atoi(prot.c_str()), offset));
			curProt = prot;
		} 
		offset = inFile.tellg();
	}
	p = offsets.begin();
	inFile.clear();
	inFile.seekg((*p).second);
	getline(inFile, str);
	indx1 = str.find_first_of(',');
	curProt = str.substr(0, indx1);
	inFile.seekg((*p).second);
}

string Offfile::nextLine()
{
	long offset;
	string str;
	size_t indx1;

	if (p != offsets.end())
	{
		if (endof == true || curProti == 0)
		{
			inFile.clear();

			offset = (*p).second;

			if (p != offsets.end())
			{
			offset = (*p).second;			
			}
			else
			{
				cout << "offset not found" << endl;
			}
			inFile.seekg(offset);		
			
			curProti = (*p).first;
			endof = false;
		}

			if(getline(inFile, str))
			{ 
				indx1 = str.find_first_of(',');
				prot = str.substr(0, indx1);
				proti = atoi(prot.c_str());

				if (curProti != proti)
				{
					++p;
					if (p != offsets.end())
					{
						offset = (*p).second;
						inFile.seekg(offset); // position to the first line in a same-prot-group
						getline(inFile, str);
						indx1 = str.find_first_of(',');
						prot = str.substr(0, indx1);
						proti = atoi(prot.c_str());
						curProti = proti;
					}
				}
				if (str != "")
				{
					return str;
				}
			}

		endof = true;
		p++;
		return(nextLine()); // don't stop on eof
	}
	offsets.clear();
	return "";
}
