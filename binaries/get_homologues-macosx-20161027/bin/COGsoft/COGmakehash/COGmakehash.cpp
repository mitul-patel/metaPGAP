#include "ereader.h"
#include <set>

#if defined(__GNUC__)

	#define SPRTR '/'

	using namespace std;
//********************** Win32
	#elif defined(_MSC_VER)
		#define SPRTR '\\'
			#else
			#error "Platform not supported. Need to update source code"
			#endif 
		


int main(int argc, char *argv[])
{
	bool overwrite = true;
	ofstream outFile;
	string str;
	string infile;
	string ofile;
	string dout = "conv";
	dout = SPRTR + dout;
	dout = '.' + dout;
	char sprtr = ',';
	int name = 1;
	//int org;
	//int cluster;
	char c;

	set<string> keys;

	while(--argc > 0 && (*++argv)[0] == '-')
    {
		c = *++argv[0];

		switch (c)
		{
		case 'h':
			cerr << "Usage: COGmakehash [Options]" << endl << endl;
			cerr << "where:" << endl << endl;
			cerr << "-i=flst		input file (list of IDs, default data.csv)" << endl << endl;
			cerr << "-o=dout		output directory (creates hash.csv, default ./conv)" << endl << endl;
			cerr << "-s=sepchar	separator character (default ,)" << endl << endl;
			cerr << "-n=nblock	index of the name field in the input file (1-based, default 1)" << endl << endl;
		//	cerr << "G index of genome block (0 for none)" << endl << endl;
		//	cerr << "C index of cluster block (0 for none)" << endl << endl;
			exit(1);
			break;
		case 'i':
			str = argv[0];
			infile = str.substr(str.find('=') + 1);
			break;
		case 'o':
			str = argv[0];
			dout = str.substr(str.find('=') + 1);
			break;
		case 's':
			str = argv[0];
			sprtr = *((str.substr(str.find('=') + 1)).data());
			break;
		case 'n':
			str = argv[0];
			name = atoi(str.substr(str.find('=') + 1).c_str());
			break;
	/*	case 'g':
			str = argv[0];
			org = atoi(str.substr(str.find('=') + 1).c_str());
			break;
		case 'c':
			str = argv[0];
			cluster = atoi(str.substr(str.find('=') + 1).c_str());
			break;*/
		case 'a':
			overwrite = false;
			break;
		default:
			break;
		}
	}

	ofile = dout + SPRTR + "hash.csv";
	
	outFile.open(ofile.c_str());
	if (!outFile)
	{
		cerr << "Can't open output file " << ofile << endl;
		exit(1);
	}
	
	KeyEnt k;

	EasyReader er(infile.c_str(), sprtr, name);

	int counter = 0;
	while(er.next(k))
	{
		if (keys.insert(k.name).second == true && k.name != "")
			outFile << ++counter << ',' << k.name << endl; //',' << k.org << ',' << k.cog << endl;
	}	
	er.close();
}
