#include <set>
#include <algorithm>
#include "graph.h"
#include "reader.h"
#include "logger.h"

typedef set<string> FltGenomeSet;
typedef map<long, KeyEnt> GeneIDOrgHash;

struct LSEGenomePrm
{
	string genomeName;
	FltGenomeSet fltGenomes;
	ofstream outStream;

	LSEGenomePrm(string& genomeName)
	{
		this->genomeName = genomeName;
	}
	LSEGenomePrm(const LSEGenomePrm& lgp)
	{
		this->genomeName = lgp.genomeName;
		this->fltGenomes = lgp.fltGenomes;
//		this->outStream = lgp.outStream;
	}
};

typedef map<string, LSEGenomePrm> LseGenomeHash;


struct MKLSE
{
	string version;
	string versionDate;
	string dataDir;
	string idsFileName;
	string orgFileName;
	string idsOrgFileName;
	string hitsFileName;
	string fhitsFileName;
	string jobFileName;
	string lseFileName;
	string tmpDir;
	bool useFilteredBlastHits;


	GeneIDOrgHash idsOrgHash;
	map<long, SelfHitEnt> selfHitsHash;

	LseGenomeHash genomes;
	ofstream lseStream;
	Logger logger;


	//----------------------------------------------
	MKLSE(): logger(cout, 5000)
	{
		version = "1.0";
		versionDate = "26 December 2006";
		dataDir = "./data";
		tmpDir = "./data/tmp";

		idsFileName = "/hash.csv";
		hitsFileName = "/hits.csv";
		fhitsFileName = "/query2subject.csv";
		idsOrgFileName = "/_hash.tmp";

		jobFileName = "";
		orgFileName = "";
		lseFileName = "";
		useFilteredBlastHits = true;
	}
	//----------------------------------------------
	void run()
	{
		logger.taskName = "Build genomes LSE";
		logger.subtaskName = "...";
		LseGenomeHash::iterator it;
		for(it = genomes.begin(); it != genomes.end(); it++)
		{
			buildGenomeLSE(it);
		}

		logger.taskName = "Write single proteins";
		logger.subtaskName = "...";
		logger.percentDone = 0;
		writeSingleProteins();
		logger.finalPrint();
	}
	//----------------------------------------------
	void initData()
	{
		lseStream.open(lseFileName.c_str());

		loadJobGenomesFile();
		mapOrgsToIds();
		loadIdsOrgHash();
		splitHitsByGenomes();
		if(useFilteredBlastHits) splitFHitsByGenomes();
	}
	//----------------------------------------------
	void deinitData()
	{
		lseStream.close();
	}

	//----------------------------------------------
	/**
	*  Generate hash.csv like file, but with information about genomes
	*/
	void mapOrgsToIds()
	{

		map<string,string> protOrgHash;
		loadProtOrgHash(orgFileName.c_str(), protOrgHash);

		logger.taskName = "Map genome names to ids";
		logger.subtaskName = "...";
		ofstream os((tmpDir + idsOrgFileName).c_str());
		LineEntryReader reader((dataDir + idsFileName).c_str());
		KeyEnt ke;
		while( reader.next(ke) )
		{
			map<string,string>::iterator it = protOrgHash.find(ke.name);
			if(it != protOrgHash.end())
			{
				ke.org = it->second;
				os << ke << "\n";
			}
			logger.percentDone = reader.readPercent();
			logger.print();
		}
		reader.close();

		os.flush();
		os.close();
		logger.finalPrint();

	}
	//----------------------------------------------
	void loadProtOrgHash(const char* fileName, map<string,string>& protOrgHash)
	{
		logger.taskName = "Load orgs hash";
		logger.subtaskName = "...";

		LineEntryReader reader(fileName);
		OrgEnt oe;
		while( reader.next(oe) )
		{
			protOrgHash[oe.protCode] = oe.orgCode;
			logger.percentDone = reader.readPercent();
			logger.print();
		}
		reader.close();
		logger.finalPrint();
	}


	//----------------------------------------------
	string getHitsTMPFileName(string& org)
	{
		return "/_hits_" + org + ".tmp";
	}
	//----------------------------------------------
	string getFHitsTMPFileName(string& org)
	{
		return "/_fhits_" + org + ".tmp";
	}
	//----------------------------------------------
	void loadJobGenomesFile()
	{
		logger.taskName = "Load job file";
		logger.subtaskName = "...";

		genomes.clear();
		LineEntryReader reader(jobFileName.c_str());
		PairEnt pe;
		while( reader.next(pe) )
		{
			string lseGenome = pe.code1;
			string fltGenome = pe.code2;
			LseGenomeHash::iterator it = genomes.find(lseGenome);
			if(it == genomes.end())
			{
				it = genomes.insert( LseGenomeHash::value_type(lseGenome, LSEGenomePrm(lseGenome)) ).first;
			}
			it->second.fltGenomes.insert(fltGenome);
			logger.percentDone = reader.readPercent();
			logger.print();
		}
		reader.close();
		logger.finalPrint();
	}

	//----------------------------------------------
	void loadIdsOrgHash()
	{
		logger.taskName = "Loading ids hash";
		logger.subtaskName = "...";
		LineEntryReader reader( (tmpDir + idsOrgFileName).c_str());
		KeyEnt ke;
		while( reader.next(ke) )
		{
			idsOrgHash[ke.id] = ke;
			logger.percentDone = reader.readPercent();
			logger.print();
		}
		reader.close();
		logger.finalPrint();
	}
	//----------------------------------------------
	void openLSEGenomeHitsTMPFiles()
	{
		LseGenomeHash::iterator it = genomes.begin();
		for(; it != genomes.end(); it++)
		{
			string outFileName = tmpDir + getHitsTMPFileName(it->second.genomeName);
			it->second.outStream.open(outFileName.c_str());
		}
	}
	//----------------------------------------------
	void openLSEGenomeFHitsTMPFiles()
	{
		LseGenomeHash::iterator it = genomes.begin();
		for(; it != genomes.end(); it++)
		{
			string outFileName = tmpDir + getFHitsTMPFileName(it->second.genomeName);
			it->second.outStream.open(outFileName.c_str());
		}
	}
	//----------------------------------------------
	void closeLSEGenomeTMPFiles()
	{
		LseGenomeHash::iterator it = genomes.begin();
		for(; it != genomes.end(); it++)
		{
			it->second.outStream.close();
		}
	}

	//----------------------------------------------
	void splitHitsByGenomes()
	{
		logger.taskName = "Splitting hits by LSE genomes";
		logger.subtaskName = "...";

		openLSEGenomeHitsTMPFiles();

		LineEntryReader reader((dataDir + hitsFileName).c_str());
		HitEnt he;

		string qOrg;
		string tOrg;
		LseGenomeHash::iterator lsegIt;
		GeneIDOrgHash::iterator qkeIt;

		long curQId = -1;
		long curTId = -1;
		bool queryShouldOmit = false;

		while( reader.next(he) )
		{
			logger.percentDone = reader.readPercent();
			logger.print();
			// Avoid hits to themselves
			if(he.query == he.target) continue;

			// To init data at the new block of query
			if(curQId != he.query)
			{
				queryShouldOmit = false;
				curQId = he.query;
				curTId = -1;

				//  define genome for given query of hit
				qkeIt = idsOrgHash.find(he.query);
				if(qkeIt == idsOrgHash.end())
				{
					queryShouldOmit = true;
					continue;
				}
				else qOrg = qkeIt->second.org;

				// if genome has not been found in LSE genomes hash => omit all hits with such query
				lsegIt = genomes.find(qOrg);
				if(lsegIt == genomes.end())
				{
					queryShouldOmit = true;
					continue;
				}
			}

			if(queryShouldOmit) continue;
			if(curTId == he.target) continue;
			curTId = he.target;

			// Check wheather the target belongs to Filter genomes
			GeneIDOrgHash::iterator tkeIt = idsOrgHash.find(he.target);
			if(tkeIt == idsOrgHash.end()) continue;
			tOrg = tkeIt->second.org;
			FltGenomeSet::iterator fltIt = lsegIt->second.fltGenomes.find(tOrg);
			if(fltIt != lsegIt->second.fltGenomes.end())
			{
				queryShouldOmit = true;
				continue;
			}

			// Means that if the target genome is the same as query, then it is our case!
			if(qOrg.compare(tOrg) == 0 )
			{
				lsegIt->second.outStream << he << endl;
			}
		}

		closeLSEGenomeTMPFiles();
		logger.finalPrint();
	}
	//----------------------------------------------
	void splitFHitsByGenomes()
	{
		logger.taskName = "Splitting fhits by LSE genomes";
		logger.subtaskName = "...";

		openLSEGenomeFHitsTMPFiles();

		LineEntryReader reader((dataDir + fhitsFileName).c_str());
		FHitEnt he;
		LseGenomeHash::iterator lsegIt;
		GeneIDOrgHash::iterator qkeIt;
		GeneIDOrgHash::iterator tkeIt;

		while( reader.next(he) )
		{
			logger.percentDone = reader.readPercent();
			logger.print();

			if(he.query == he.target) continue;

			qkeIt = idsOrgHash.find(he.query);
			if(qkeIt == idsOrgHash.end()) continue;
			lsegIt = genomes.find(qkeIt->second.org);
			if(lsegIt == genomes.end()) continue;

			tkeIt = idsOrgHash.find(he.target);
			if(tkeIt == idsOrgHash.end()) continue;


			if(qkeIt->second.org.compare(tkeIt->second.org) == 0 )
			{
				lsegIt->second.outStream << he << endl;
			}
		}
		closeLSEGenomeTMPFiles();
		logger.finalPrint();
	}


	//----------------------------------------------
	void buildGenomeLSE(LseGenomeHash::iterator& gmIt)
	{
		logger.subtaskName = gmIt->second.genomeName;
		logger.percentDone = 0;

		string fileName;
		Graph graph;

		fileName = tmpDir + getHitsTMPFileName(gmIt->second.genomeName);
		graph.load(fileName.c_str());

		if( useFilteredBlastHits )
		{
			fileName = tmpDir + getFHitsTMPFileName(gmIt->second.genomeName);
			graph.filterByEdgePresence(fileName.c_str());
		}

		vector<ConnectedComponent> ccomps;
		graph.collectConnectedComponents(ccomps);

		writeComponents(ccomps, true);
		logger.finalPrint();
	}
	//----------------------------------------------
	void writeComponents(vector<ConnectedComponent>& ccomps, bool removeProcessedIds)
	{
		vector<ConnectedComponent>::iterator it = ccomps.begin();
		for(; it != ccomps.end(); it++)
		{
			ConnectedComponent::iterator jt =  it->begin();
			for(; jt != it->end(); jt ++)
			{
				if(jt != it->begin()) lseStream << ",";

				GeneIDOrgHash::iterator gIt = idsOrgHash.find((*jt)->id);
				if(gIt != idsOrgHash.end())
				{
					lseStream << gIt->second.name;
					if(removeProcessedIds) idsOrgHash.erase(gIt);
				}
			}
			lseStream << endl;
		}
	}
	//----------------------------------------------
	void writeSingleProteins()
	{
		GeneIDOrgHash::iterator it = idsOrgHash.begin();
		for(; it != idsOrgHash.end(); it++)
		{
			lseStream << it->second.name << endl;
		}
	}
	//----------------------------------------------
	bool readParams(int argc, char* argv[])
	{
		for(int i = 1; i < argc; i++)
		{
			string str = argv[i];
			if(str.compare(0,3,"-j=") == 0)
			{
				jobFileName = str.substr(3);
			}
			else if( str.compare(0,3,"-p=") == 0 )
			{
				orgFileName = str.substr(3);
			}
			else if( str.compare(0,3,"-o=") == 0 )
			{
				lseFileName = str.substr(3);
			}
			else if ( str.compare(0,3,"-d=") == 0 )
			{
				dataDir = str.substr(3);
			}
			else if ( str.compare(0,3,"-t=") == 0 )
			{
				tmpDir = str.substr(3);
			}
			else if ( str.compare ("-e") == 0 )
			{
				useFilteredBlastHits = false;
			}
		}

		return (jobFileName.size() > 0)
			&& (orgFileName.size() > 0)
			&& (lseFileName.size() > 0);
	}
	//----------------------------------------------
	void printMAN()
	{
		/*
		cout
		 << "\n"
		 << "mklse " << version << "; " << versionDate << "\n"
		 << "\n"
		 << "Generates LSEs for each genome listed in <lsetask> file \n"
		 << "basing on blastp results.\n"
		 << "\n"
		 << "Options:\n"
		 << "\n"
		 << " -j=filename\n"
		 << "     File sets the task for LSE generator. It describes the set\n"
		 << "     of genomes for which LSEs should be generated, and list of\n"
		 << "     filter genomes (for each LSE genome) that should be used\n"
		 << "     to filter out too weak blastp hits. Each line of <lsetask>\n"
		 << "     file should contain two comma separated columns: first \n"
		 << "     column - the name of genome for which LSE should be \n"
		 << "     generated; the second column - the name of genome, which\n"
		 << "     should be used as a filter genome. More than one line can be\n"
		 << "     specified for one LSE genome.\n"
		 << "\n"
		 << " -p=filename\n"
		 << "     The standard <prot2org> file that maps external (real) protein \n"
		 << "     ids to genome names.\n"
		 << "\n"
		 << " -o=filename\n"
		 << "     Output file name \n"
		 << "\n"
		 << "  -d=dirname    (default ./data)\n"
		 << "     The full path of data directory. This directory should contain \n"
		 << "     predefined set of files, mainly the results of blastp analysis. \n"
		 << "     The list of files:\n"
		 << "     hits.csv          - file with blast hits (internal protein ids \n"
		 << "                         used)\n"
		 << "     query2subject.csv - file which shows the presence of query to\n"
		 << "                          target hits in filtered balstp(internal \n"
		 << "                          protein ids used)\n"
		 << "     hash.csv          - correspondance between internal and external(real\n"
		 << "                         gis, sw...) protein ids\n"
		 << "\n"
		 << "  -t=dirname    (default ./data/tmp)\n"
		 << "     The full path of temporary directory.\n"
		 << "\n"
		 << "  -e\n"
		 << "     Do not use filtered balstp data from query2subject.csv file\n"
		 << "" << endl;
		 */
		cout
		 << "\n"
		 << "Usage: COGlse [Options]\n"
		 << "\n"
		 << "Options:\n"
		 << "\n"
		 << "-d=dconv\tdirectory for converted data (must contain hash.csv; default ./conv)\n"
		 << "\n"
		 << "-j=jobfile\tLSE maker job description (default jobLSE.csv)\n"
		 << "\n"
		 << "-p=fqp2o\tquery file (protein-to-organism data; default prot2org.csv)\n"
		 << "\n"
		 << "-e\t\tignore the filtered search results\n"
		 << "\n"
		 << "-o=fout\t\toutput file (default lse.csv)\n"
		 << "\n"
		 << "-t=dtemp\ttemporary directory (default ./data/tmp)\n"
		 << "" << endl;

	}
};


//----------------------------------------------
int main(int argc, char* argv[])
{
	MKLSE mklse;
	if( mklse.readParams(argc, argv) )
	{
		mklse.initData();
		mklse.run();
		mklse.deinitData();
	}
	else
	{
		mklse.printMAN();
	}
	return 0;
}

