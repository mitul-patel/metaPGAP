//The program uses OS-specific class BeadCounter (see bc.h)

//#define DEBUG

#ifdef DEBUG
	#include <stdio.h>
	#include <time.h>
	clock_t start, end;
	void startClock(){ start = clock(); }
	void finishClock()
	{
		end = clock();
		printf("clocks: %f\n",((double)(end-start))/CLOCKS_PER_SEC);
	}
#endif // DEBUG

#include <iostream>
#include <fstream>
#include <vector> 
#include <set>
#include <map>
#include <algorithm>
#include <list>

#include "blastconv.h"
#include "reader.h"

using namespace std;
const char* fMaskBlast = FMASK;//mandatory extension of input BLAST files is .tab; FMASK is OS-dependent and defines in bc.h
vector<string> v;

int main(int argc, char *argv[])
{
#ifdef DEBUG
	startClock();
#endif
    string dirout = ".";
	dirout = dirout + SPRTR + "conv";
    string unfiltered = "blan";
    string filtered = "blaf";
    string self = "blan";
    double e = 10.0;
	bool quiet = true;
	bool hash = true;
	bool overwrite = true;
	bool reverse = false;

	char c;
    string fileName;
    string str;
    string query;
    string subject;
    string qStart;
    string qEnd;
    string tStart;
    string tEnd;
    string evalue;
    double score;
    int Q;
    int T;
    int qID = 2;
    int tID = 2;
	int len;
	int nProt = 0;

	map<string, int> domHash;
	map<string, int>::iterator p_m;

	while(--argc > 0 && (*++argv)[0] == '-')
    {
		c = *++argv[0];

		switch (c)
		{
		case 'd':
			str = argv[0];
			dirout = str.substr(str.find('=') + 1);
			break;
		case 'u':
			str = argv[0];
			unfiltered = str.substr(str.find('=') + 1);
			break;
		case 'f':
			str = argv[0];
			filtered = str.substr(str.find('=') + 1);
			break;
		case 's':
			str = argv[0];
			self = str.substr(str.find('=') + 1);
			break;
		case 'e':
			str = argv[0];
			e = atof(str.substr(str.find('=') + 1).c_str());
			break;
		case 'q':
			str = argv[0];
			Q = atoi(str.substr(str.find('=') + 1).c_str());
			qID = Q > 0 ? Q : qID;
			break;
		case 't':
			str = argv[0];
			T = atoi(str.substr(str.find('=') + 1).c_str());
			tID = T > 0 ? T : tID;
			break;
		case 'v':
			quiet = false;
			break;
		case 'a':
			overwrite = false;
			break;
		case 'n':
			hash = false;
			break;
		case 'r':
			reverse = true;
			break;
		case 'h':
			commandLineHelp();
			break;
		default:
			break;
		}//switch (c)
    }

    const char* dUnfiltered = unfiltered.c_str();
    const char* dFiltered = filtered.c_str();
    const char* dSelf = self.c_str();

    string curPathU = dirout + SPRTR + "hits.csv";
    string curPathF = dirout + SPRTR + "query2subject.csv";
    string curPathS = dirout + SPRTR + "self.csv";
    string keys = dirout + SPRTR + "hash.csv";

    fstream fileU;
	fstream fileF;		
    fstream fileS;

	if (overwrite == true)
	{
		fileU.open(curPathU.c_str(), ios::out);
		fileF.open(curPathF.c_str(), ios::out);
		fileS.open(curPathS.c_str(), ios::out);
	}
	else
	{
		fileU.open(curPathU.c_str(), ios::out | ios::app);
		fileF.open(curPathF.c_str(), ios::out | ios::app);
		fileS.open(curPathS.c_str(), ios::out | ios::app);
	}

    if (!fileU)
    {
    	cerr << "Error opening file: " << curPathU << endl;
    	exit(1);
    }
    if (!fileF)
    {
       	cerr << "Error opening file: query2subject.csv" << endl;
		exit(1);
    }
    if (!fileS)
    {
		cerr << "Error opening file: self.csv" << endl;
		exit(1);
    }

	KeyReader kr(keys.c_str());
	KeyEnt ke;
	while (kr.next(ke))
	{
		domHash.insert(pair<string, int>(ke.name, ke.id));
	}
	
	string curPathL = dirout + SPRTR + "COGreadblast.log";
	ofstream log(curPathL.c_str());
	if (!log)
	{
		cerr << "Error opening log file: COGreadblast.log" << endl;
		exit(1);
	}

	if (quiet == true)
	{
		streambuf *psbuf;
		psbuf = log.rdbuf();
		cout.rdbuf(psbuf);
	}

	cout << "Processing unfiltered hits" << endl;
	{
	set<string> domains;
	BeadCounter bcy(dUnfiltered, fMaskBlast);
    do
    {
		fileName = bcy.getFile();

		ifstream inFile(fileName.c_str());
		if (!inFile)
		{
			cerr << "Error opening input unfiltered blast file " << fileName << endl;
			exit(1);
		}
		else
		{
			while (getline(inFile, str))
			{
				if (str[0] != '#')
				{
					if (split(str, "\t", &v) != 12)
					{
						cout << "corrupted BLAST line in " << fileName << endl;
						cout << fileName << " processing abandoned, go to next file" << endl;
						break;
					}

					evalue = v[10];
					if (atof(evalue.c_str()) > e)
						continue;

					query = truncate(v[0], qID);
					if (domHash.find(query) == domHash.end())
						continue;

					subject = truncate(v[1], tID);
					if (domHash.find(subject) == domHash.end())
						continue;

					qStart = v[6];
					qEnd = v[7];
					tStart = v[8];
					tEnd = v[9];
					evalue = v[10];
					score = atof((v[11]).c_str());
		
					if (hash == true)
					{
						fileU << (*(domHash.find(query))).second << ',' << (*(domHash.find(subject))).second << ',' << qStart << ',' << qEnd << ',' << tStart << ',' << tEnd << ',' << evalue << ',' << score  << endl;
					}
					else
					{
						fileU << query << ',' << subject << ',' << qStart << ',' << qEnd << ',' << tStart << ',' << tEnd << ',' << evalue << ',' << score  << endl;
					}
				}
			}
			inFile.close();
			inFile.clear();
			inFile.open(bcy.getFile());
		}
    }
    while(bcy.findNextFile());


	}
	
	cout << "Processing filtered hits" << endl;
	{
	BeadCounter bcf(dFiltered, fMaskBlast);
    do
    {
		fileName = bcf.getFile();

		ifstream inFile(fileName.c_str());
		if (!inFile)
		{
			cerr << "Error opening input filtered blast file " << fileName << endl;
			exit(1);
		}
		else
		{
			while (getline(inFile, str))
			{
				if (str[0] != '#')
				{
					if (split(str, "\t", &v) != 12)
					{
						cout << "corrupted BLAST line in " << fileName << endl;
						cout << fileName << " processing abandoned, go to next file" << endl;
						break;
					}

					evalue = v[10];
					if (atof(evalue.c_str()) > e)
						continue;

					query = truncate(v[0], qID);
					if (domHash.find(query) == domHash.end())
						continue;

					subject = truncate(v[1], tID);
					if (domHash.find(subject) == domHash.end())
						continue;
					
					if (hash == true)
					{
						fileF << (*(domHash.find(query))).second << ',' << (*(domHash.find(subject))).second << endl;
					}
					else
					{	
						fileF << query << ',' << subject << endl;
					}
				}				

			}	
			inFile.close();
			inFile.clear();
			inFile.open(bcf.getFile());
		}

	}
    while(bcf.findNextFile());
	}


 	cout << "Processing self hits" << endl; 
	{
    string curQuery = "";
    int lenBest;
    double scoreBest;
    bool found = false;

	BeadCounter bcx(dSelf, fMaskBlast);
    do
    {
		fileName = bcx.getFile();

		ifstream inFile(fileName.c_str());
		if (!inFile)
		{
			cerr << "Error opening input filtered blast file " << fileName << endl;
			exit(1);
		}
		else
		{
			while (getline(inFile, str))
			{
				if (str[0] != '#')
				{
					if (split(str, "\t", &v) != 12)
					{
						cout << "corrupted BLAST line in " << fileName << endl;
						cout << fileName << " processing abandoned, go to next file" << endl;
						break;
					}

					query = truncate(v[0], qID);
					if (domHash.find(query) == domHash.end())
						continue;

					subject = truncate(v[1], qID);//use here qID instead of tID
					if (domHash.find(subject) == domHash.end())
						continue;

					if (curQuery != query)
					{//gets here only with first hit into query
						if (found == true)
						{//self hit for the previous query was found (normal situation)
							found = false;
						}
						else if (curQuery != "") 
						{//if self hit for the previous query missed for some reason then write best hit for the previous query
							fileS << (*domHash.find(curQuery)).second << ',' << lenBest << ',' << scoreBest << endl;
							nProt++;
						}

						scoreBest = atof((v[11]).c_str());
						lenBest = atoi((v[7]).c_str()) - atoi((v[6]).c_str()) + 1;
						curQuery = query;
					}// if (curQuery != query)
					
					if (query == subject)
					{	
						qStart = v[6];
						qEnd = v[7];
						score = atof((v[11]).c_str());
						len = atoi((v[7]).c_str()) - atoi((v[6]).c_str()) + 1;
						fileS << (*domHash.find(query)).second << ',' << len << ',' << score  << endl;
						nProt++;
						found = true;	
					}
				}// if (str[0] != '#')
			}// while(getline(inFile, str))
				
			if (found == false && nProt != 0)
			{//if self hit for the last query missed for some reason then write best hit for the previous query
				fileS << (*domHash.find(curQuery)).second << ',' << lenBest << ',' << scoreBest << endl;
				nProt++;
			}

			inFile.close();
			inFile.clear();
			inFile.open(bcx.getFile());
	    }	    

	}
    while(bcx.findNextFile());
	}
 
	cout << endl;
	cout << nProt << " self hits found" << endl;
#ifdef DEBUG
finishClock();
#endif
	//*************sort, add and delete
	if (overwrite == false)
	{
		string curPathUU = dirout + SPRTR + "hitsUnique.csv";
		string curPathFU = dirout + SPRTR + "q2sUnique.csv";
		string curPathSU = dirout + SPRTR + "selfUnique.csv";
		string curPathUS = dirout + SPRTR + "hitsSorted.csv";
		string curPathFS = dirout + SPRTR + "q2sSorted.csv";

		domHash.clear();
		multimap<int, long> offsets;
		multimap<int, long>::iterator pm;

		Hit* h;
		QnS* s;

		int q64 = 0;
		int curQ64 = 0;
//		long offset;
		string boof;

		fileU.close();
		fileF.close();
		fileS.close();

		fileU.open(curPathU.c_str(), ios::in);
		fileF.open(curPathF.c_str(), ios::in);
		fileS.open(curPathS.c_str(), ios::in);	
		
		fstream fileUUniq;
		fstream fileFUniq;
#ifdef DEBUG
startClock();
#endif
{
		Offfile of(reverse, curPathU.c_str());
		HitSet hits4q;

		fileUUniq.open(curPathUU.c_str(), ios::out);


		while("" != (boof = of.nextLine()))
		{
			split(boof, ",", &v);
			q64 = atoi(v[0].c_str());
				if (curQ64 != q64)
			{
				hits4q.correct();
				hits4q.sortOnScore();
				hits4q.print(&fileUUniq);
				curQ64 = q64;
				hits4q.clear();
			}
			h = new Hit(atoi(v[0].c_str()), atoi(v[1].c_str()), atol(v[2].c_str()), atol(v[3].c_str()), atol(v[4].c_str()), atol(v[5].c_str()), atof(v[6].c_str()), atof(v[7].c_str()));
			hits4q.insert(*h);
			delete h;
		}
		fileUUniq.close();
}// end of Offile of lifetime
#ifdef DEBUG
finishClock();
#endif
		system((COPY + (' ' + curPathUU + ' ' + curPathU)).c_str());
		system((DEL + (' ' + curPathUU)).c_str());
#ifdef DEBUG
startClock();
#endif

		//************* filtered
{
		Offfile of1(reverse, curPathF.c_str());
		QSet subs4q;

		fileFUniq.open(curPathFU.c_str(), ios::out);

		while("" != (boof = of1.nextLine()))
		{
			split(boof, ",", &v);
			q64 = atoi(v[0].c_str());
				if (curQ64 != q64)
			{
				subs4q.correct();
				subs4q.print(&fileFUniq);
				curQ64 = q64;
				subs4q.clear();
			}
			s = new QnS(atoi(v[0].c_str()), atoi(v[1].c_str()));
			subs4q.insert(*s);
			delete s;
		}
		fileFUniq.close();
}// end of Offfile of1 lifetime
#ifdef DEBUG
finishClock();
#endif
		system((COPY + (' ' + curPathFU + ' ' + curPathF)).c_str());
		system((DEL + (' ' + curPathFU)).c_str());

		//********************* Self

		fstream fileSUniq;
		fileSUniq.open(curPathSU.c_str(), ios::out);
		map<string, string> selfhits;

		map<string, string>::iterator msp;
		double selfscore1;
		double selfscore2;
		string str2;

		while(getline(fileS, str))
		{
			size_t indx = str.find_first_of(",");
			string prot = str.substr(0, indx);
			if (selfhits.find(prot) == selfhits.end())
			{
				selfhits.insert(pair<string, string>(prot, str));
			}
			else
			{
			
				selfscore1 = atof((str.substr(str.find_last_of(",") + 1)).c_str());
				msp = selfhits.find(prot);
				str2 = (*msp).second;
				
				selfscore2 = atof((str2.substr(str2.find_last_of(",") + 1)).c_str());

				if (selfscore1 > selfscore2)
				{
					(*msp).second = str;
				}
			}
				
		}

		msp = selfhits.begin();
		while (msp != selfhits.end())
		{
			fileSUniq << (*msp).second << endl;
			msp++;
		}

		fileS.close();
		fileSUniq.close();

		system((COPY + (' ' + curPathSU + ' ' + curPathS)).c_str());
		system((DEL + (' ' + curPathSU)).c_str());
	}
	return 0;
}
