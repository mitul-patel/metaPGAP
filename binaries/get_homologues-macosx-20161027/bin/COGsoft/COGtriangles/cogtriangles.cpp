// all proteins to be processed MUST be found in prot2org input file AND lse.csv file unless
// there is no LSE input file (the program still has to be runnable)

#include "cogmaker.h"
#include "os.h"

#ifdef DEBUG
	Stopwatch sw;	
#endif

int main(int argc, char *argv[])
{
	double	expectThreshold	= 10;
	double	overlapThreshold = 0.75;
	bool quiet = true;

    string dirpath = ".";
	dirpath = dirpath + SPRTR + "conv";
	string self = "self.csv";
	string fout ="out.csv";
	string q2o = "prot2org.csv";
	string lsef = "lse.csv";
	string cogName = "c.01.";
	int startNum = 1;

	map<string, int> domHash; 
	size_t indx;
	string subj;
	int numID;

	string inFile;

	char	c;
	string	str;
	double	oT;
	double	eT;

	while(--argc > 0 && (*++argv)[0] == '-') 
	{
		c = *++argv[0];		
		switch (c)
		{
			case 'i': // input directory (processed BLAST)
				str = argv[0];
				dirpath = str.substr(str.find('=') + 1);
				break;
			case 'q': // queries2orgs
				str = argv[0];
				q2o = str.substr(str.find('=') + 1);
				break;
			case 'l':
				str = argv[0];
				lsef = str.substr(str.find('=') + 1);
				break;
			case 'o': // output file
				str = argv[0];
				fout = str.substr(str.find('=') + 1);
				break;
			case 't':
				str = argv[0];
				oT = atof(str.substr(str.find('=') + 1).c_str());
				overlapThreshold = oT > 0 && oT <= 1 ? oT : overlapThreshold;
				break;
			case 'e':
				str = argv[0];
				eT = atof(str.substr(str.find('=') + 1).c_str());
				expectThreshold = eT > 0 ? eT : expectThreshold;
				break;
			case 'n':
				str = argv[0];
				cogName = str.substr(str.find('=') + 1);
				break;
			case 's':
				str = argv[0];
				startNum = atoi(str.substr(str.find('=') + 1).c_str());
				break;
			case 'v':
				quiet = false;
				break;
			default:
				cerr << "Usage: COGtriangles [Options]" << endl << endl;

				cerr << "Options:" << endl << endl;

				cerr << "-i=dconv	input directory (must contain hash.csv; default ./conv)" << endl << endl;

				cerr << "-q=fqp2o	query file (protein-to-organism data; default prot2org.csv)" << endl << endl;

				cerr << "-l=flse		file with lineage-specific expansions data (default lse.csv)" << endl << endl;

				cerr << "-o=fout		output file (default out.csv)" << endl << endl;

				cerr << "-t=othr		hit coverage threshold (default 0.75)" << endl << endl;

				cerr << "-e=evalue	threshold for a BLAST hits (default 10)" << endl << endl;

				cerr << "-n=name		constant part of the claster name (default c.01.)" << endl << endl;

				cerr << "-s=number	first number in the numerical part of cluster name (default 1)" << endl << endl;

				cerr << "-v		verbose mode (mostly debugging output to STDOUT)" << endl << endl;
				exit(1);
		}
	}

	string curPathL = dirpath + SPRTR + "COGtriangles.log";
	ofstream log(curPathL.c_str());
	if (!log)
	{
		cerr << "Error opening output file: COGtriangles.log" << endl;
		exit(1);
	}

	if (quiet == true)
	{
		streambuf *psbuf;
		psbuf = log.rdbuf();
		cout.rdbuf(psbuf);
	}

	#ifdef DEBUG
		cout << "Start the program" << endl;
		sw.startClock();
	#endif

	ofstream fOut(fout.c_str());
    if (!fOut)
    {
    	cerr << "Error opening output file: " << fout << endl;
    	exit(1);
    }

	inFile = dirpath + SPRTR;

	ifstream inDomHash((inFile + "hash.csv").c_str());
	if (!inDomHash)
	{
		cout << "File " << inFile + "hash.csv not found" << endl;
		cout << "Native target IDs must be numeric" << endl;
		cout << "Trying to proceed..." << endl << endl;
	}
	else
	{
		while(getline(inDomHash, str))
		{
			indx = str.find_last_of(',');
			subj = str.substr(indx + 1);
			numID = atoi((str.substr(0, indx)).c_str());
			domHash.insert(pair<string, int>(subj, numID));
		}
	}





    map<int, string>	prot2org;
	map<int, int>		prot2len;
	map<int, long>		offsetsU; // Presuming all hits for any query are successive
	map<int, long>		offsetsF;
	map<int, Lse>		allLse;
	map<int, int>		prot2lse;


	ifstream inP2O((q2o).c_str()); 
	if (!inP2O)
	{
		cerr << "Error opening input file: " << q2o << endl;
		exit(1);
	}

	populate(&domHash, &prot2org, inP2O);

	ifstream inSelf((inFile + self).c_str()); 
	if (!inSelf)
	{
		cerr << "Error opening input file: self.csv" << endl;
	}
	else
	{
		populate(&prot2len, inSelf);
	}

	ifstream inQ2S((inFile + "query2subject.csv").c_str());
	if (!inQ2S)
	{
		cout << "File " << inFile + "query2subect.csv" << " not found; all hits will be processed" << endl;
	}
	else
	{
		makeOffsets(&offsetsF, inQ2S);
		if (offsetsF.size() == 0)
		{
			cout << "File query2subject.csv is empty; all hits will be processed" << endl;
		}
	}

	ifstream inHitFile((inFile + "hits.csv").c_str());
	if (!inHitFile)
	{
		cerr << "Error opening input file: " << inFile + "hits.csv" << endl;
		exit(1);
	}
	else
	{
		makeOffsets(&offsetsU, inHitFile);
		if (offsetsU.size() == 0)
		{
			cerr << "No hits found" << endl;
			exit(1);
		}
	}

	ifstream inLseFile((lsef).c_str());
	if (!inLseFile)
	{
		cerr << "Error opening input file: " << lsef << endl;
		exit(1);
	}
	else
	{
		populate(&domHash, &prot2lse, &allLse, inLseFile);
	}
	cout << "LSE completed" << endl;

#ifdef DEBUG
	cout << "Number of LSE: " << allLse.size() << endl;
	cout << "Ready for hitset" << endl;
	sw.finishClock();
	sw.startClock();
#endif

	HitSet hs(
		expectThreshold,
		overlapThreshold,
		&domHash,
		&prot2org,
		&prot2len,
		&offsetsF,
		&offsetsU,
		&allLse,
		&prot2lse,
		cogName,
		startNum,
		inFile);


	map<int, Lse>::iterator p_m = allLse.begin();
	int counter = 0;
	while (p_m != allLse.end())
	{
		++counter;
		hs.insert(&(*p_m).second.lse);
		hs.clear();
		p_m++;
	}
	cout << "Hitset completed" << endl;

	#ifdef DEBUG
		sw.finishClock();
	#endif

	hs.makeCOGs(&fOut);

//	p_m = allLse.begin();
//	while (p_m != allLse.end())
//	{
//		hs.printRest(&fOut, &(*p_m).second.lse);
//	}
	exit(EXIT_SUCCESS);
}
