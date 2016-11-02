#include "blastconvglob.h"

void commandLineHelp()
{
	cerr << "Usage: COGreadblast [Options]" << endl << endl;

	cerr << "Otptions:" << endl << endl;
	
	cerr << "-d=dconv	directory for converted data (must contain hash.csv; default ./conv)" << endl << endl;
	
	cerr << "-u=dunfilt	directory with the unfiltered BLAST results (default ./blan)" << endl << endl;
	
	cerr << "-f=dfilt	directory with the filtered BLAST results (default ./blaf)" << endl << endl;
	
	cerr << "-s=dself	directory with the self-BLAST results (default ./blan)" << endl << endl;
	
	cerr << "-e=evalue	e-value threshold for BLAST hits (default 10)" << endl << endl;
	
	cerr << "-q=nblock	index of the sequence ID field for the BLAST query (default 2)" << endl << endl;
	
	cerr << "-t=nblock	index of the sequence ID field for the BLAST target (default 2)" << endl << endl;
	
	cerr << "-a		append/aggregate mode (use if BLAST hits from one query do not form " << endl; 
	cerr << "		a contiguous block in the BLAST output files)" << endl << endl;
	
	cerr << "-r		symmetrize reciprocal hits (use when BLAST search has not been run " << endl; 
	cerr << "		in a fully symmetrical all-against-all manner)" << endl << endl;
	
	cerr << "-v		verbose mode (mostly debugging output to STDOUT)" << endl << endl;
	exit(1);
}

int split(const string str, const string sepStr, vector<string>* v)
{
    string	token;
    int		nToken = 0;
    size_t	index;
    size_t	posSep;
    size_t	strLen;
    size_t	tokenLen;

	(*v).clear();

    strLen = str.length();
	if (strLen == 0)
		return 0;

    index = str.find_first_not_of(sepStr);// skip leading separators

    do
    {
	posSep = str.find_first_of(sepStr, index);
	if (posSep == string::npos)
	{// separator not found
	    tokenLen = strLen - index;
	}
	else
	{
	    tokenLen = posSep - index;
	}
	token = str.substr(index, tokenLen);
	(*v).push_back(token);
	nToken++;
	//index = str.find_first_not_of(sepStr, posSep + 1);//skip all contin. separators
	index = posSep + 1;
    }
	while(posSep != string::npos && index != string::npos);
    
    return nToken;
}

string truncate(const string name, const int ID)
{
	int count;
	string trName = name;

	if (trName.find('|') != string::npos && ID != 0) 
	{
		count = 1;
		
		while (count != ID)
		{
			trName = trName.substr(trName.find('|') + 1);
			count++;
		}
		if (trName.find('|') != string::npos)	
		{
			trName = trName.substr(0, trName.find('|'));
		}
	}
	return trName;
}

void makeOffsets(multimap<int, long>* offsets, bool reverse, fstream& inFile)
{
	long offset;
	string str;
	string prot;
	string subjct;
	size_t indx1;
	size_t indx2;

	(*offsets).clear();
	
	offset = inFile.tellg();
	while (getline(inFile, str))
	{
		indx1 = str.find_first_of(",");
		indx2 = str.find_first_of(",", indx1 + 1);
		prot = str.substr(0, indx1);
		if (indx2 != str.npos)
			subjct = str.substr(++indx1, indx2 - indx1);
		else
			subjct = str.substr(++indx1);

		if (prot <= subjct || reverse == false)
		{
			(*offsets).insert(pair<int, long>(atoi(prot.c_str()), offset));
		}
		else
		{
			(*offsets).insert(pair<int, long>(atoi(subjct.c_str()), offset));
		}
		offset = inFile.tellg();
	}
}

void makeOffsetsUniq(multimap<int, long>* offsets, fstream& inFile)
{
	long offset;
	string str;
	string prot;
	string subjct;
	size_t indx; 

	(*offsets).clear();
	
	offset = inFile.tellg();
	while (getline(inFile, str))
	{
		indx = str.find_first_of(",");
		prot = str.substr(0, indx);
		subjct = str.substr(indx + 1);

		(*offsets).insert(pair<int, long>(atoi(prot.c_str()), offset));

		offset = inFile.tellg();
	}
}
