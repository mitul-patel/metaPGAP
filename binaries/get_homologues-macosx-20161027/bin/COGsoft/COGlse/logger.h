#ifndef _LOGGER_H
#define _LOGGER_H

#include <ostream>
#include <fstream>
#include <string>



class Logger
{
public:
	string taskName;
	string subtaskName;
	int percentDone;

private:
	ostream& os;

	bool enabled;
	int iterIndex;
	int printIndex;
	int period;

public:
	//----------------------------------------------
	Logger(ostream& _os, int period): os(_os)
	{
		this->period = period;
		init();
	}

	//----------------------------------------------
	void setEnabled(bool enabled)
	{
		this->enabled = enabled;
	}	

	//----------------------------------------------
	void print()
	{
		print(false);
	}


	//----------------------------------------------
	void finalPrint()
	{
		print(true);
		os << "\n" << flush;
		init();
	}
private:
	//----------------------------------------------
	void init()
	{
		iterIndex = 0;
		printIndex = 0;
	}

	//----------------------------------------------
	void print(bool isLast)
	{
		char procSymbols[] = "|/-\\";
		if(enabled)
		{
			if(isLast)
			{
				os  << "\r "
				    << " "
					<< "  " << taskName 
					<< "\t"   << subtaskName
					<< "\t" << 100 << "% " 
					<< flush;

			}
			else if( ((iterIndex++) % period) == 0)
			{
				os << "\r "
					<< procSymbols[(++printIndex) %4]
					<< "  " << taskName 
					<< "\t"   << subtaskName
					<< "\t" << percentDone << "% " 
					<< flush;
			}
		}
	}
	
};

#endif
