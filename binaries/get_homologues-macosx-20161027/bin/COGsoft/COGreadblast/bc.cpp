#include "bc.h"

#if defined(__GNUC__)

	void myOpenDir(const char* dirpath)
	{
		chdir(dirpath);
	}

	BeadCounter::BeadCounter(const char* dirPath, const char* dMask, const char* fMask)
	{
		counter = 0;
		bool success = false;

		fileMask = fMask;
		dir = opendir(dirPath);
		if (dir == NULL)
		{
		cerr << "Unable to locate the directory: " << dirPath << endl;
		exit(1);
		}
		else
		{
		subdirName = string(dirPath) + "/" + dMask + "/";
		
		dir = opendir(subdirName.c_str()); //open desired subdirectory
		if (dir == NULL)
		{
			cerr << "Unable to locate subdirectory: " << subdirName << endl;
			exit(1);
		}
		else
		{
			while(NULL != (dp = readdir(dir)))
			{
			curFilePath = subdirName + dp->d_name;

			if (curFilePath.find(fileMask) != string::npos)
			{
				success = true;
				counter++;
				break;
			}
			}
			if (success == false)
			{
			cerr << "Directory " << subdirName << " does not contain " << fileMask << " files" << endl;
			exit(1);
			}
		}
		}
	}

	BeadCounter::BeadCounter(const char* dirPath, const char* fMask) // if there're no subdirectoryies
	{
		counter = 0;
		bool success = false;

		fileMask = fMask;
		subdirName = string(dirPath) + "/";
		dir = opendir(dirPath);
		if (dir == NULL)
		{
		cerr << "Unable to locate the directory: " << dirPath << endl;
		exit(1);
		}
		else
		{
			while(NULL != (dp = readdir(dir)))
			{
			curFilePath = subdirName + dp->d_name;

			if (curFilePath.find(fileMask) != string::npos)
			{
				success = true;
				counter++;
				break;
			}
			}
			if (success == false)
			{
			cerr << "Directory " << subdirName << " does not contain " << fileMask << " files" << endl;
			exit(1);
			}
		}
	}

	BeadCounter::~BeadCounter()
	{
		closedir(dir);

		if (counter > 1)
			cout << "found " << counter << " " << fileMask << " files" << endl;
		else if (counter = 1)
			cout << "found " << counter << " " << fileMask << " file" << endl;
	}

	const char* BeadCounter::getFile() const
	{
		return curFilePath.c_str();
	}

	bool BeadCounter::findNextFile()
	{
		while(NULL != (dp = readdir(dir)))
		{
		curFilePath = subdirName + dp->d_name;

		if (curFilePath.find(fileMask) != string::npos)
		{
			counter++;
			return true;		    
		}

		}
		return false;
	}

	#elif defined(_MSC_VER)

		void myOpenDir(const char* dirpath)
		{
			_chdir(dirpath);
		}

BeadCounter::BeadCounter(const char* dirPath, const char* dMask, const char* fMask)
{
	counter = 0;

   if( _getcwd(startDir, _MAX_PATH ) == NULL )
   {
	   cerr << "Can't get path to the current directory" << endl;
	   exit(1);
   }

	fileMask = new char[strlen(fMask) + 1];
	strcpy(fileMask, fMask);

	if (_chdir(dirPath) == -1L) //change working directory
	{
		cerr << "Unable to locate the directory: " << dirPath << endl;
		exit(1);
	}

	if ((hSubdir = _findfirst(dMask, &subdir)) == -1L)	//Find first subdirecotory
	{
		cerr << "The directory " << dirPath << " is empty" << endl;
		exit(1);
	}
	else
	{
		string subdirName = subdir.name;
		string curPath = ".\\" + subdirName; //full path to directory
		
		_chdir(curPath.c_str());

		if ((hFile = _findfirst(fMask, &file)) == -1L) // Find first file in current directory 
		{
			cerr << "No " << fMask << " files in subdirectory " << curPath << endl;
			exit(1);
		}
		else
		{
			curFilePath = file.name;
			counter++;
		}
	}
}

BeadCounter::BeadCounter(const char* dirPath, const char* fMask) // if there're no subdirectoryies
{
	hSubdir = -1L;
	counter = 0;

	if( _getcwd(startDir, _MAX_PATH ) == NULL )
	{
		cerr << "Can't get path to the current directory" << endl;
		exit(1);
	}

	fileMask = new char[strlen(fMask) + 1];
	strcpy(fileMask, fMask);

	if (_chdir(dirPath) == -1L) //change working directory
	{
		cerr << "Unable to locate the directory: " << dirPath << endl;
		exit(1);
	}

	if ((hFile = _findfirst(fMask, &file)) == -1L) // Find first file in current directory 
	{
		cerr << "No " << fMask << " files in directory " << dirPath << endl;
		exit(1);
	}
	else
	{
		curFilePath = file.name;
		counter++;
	}
}

BeadCounter::~BeadCounter()
{
	if (_chdir(startDir) == -1L)
	{
		cerr << "Unable to return to the directory: " << startDir << endl;
		exit(1);
	}
	if (counter > 1)
		cout << "found " << counter << " " << fileMask << " files" << endl;
	else if (counter = 1)
		cout << "found " << counter << " " << fileMask << " file" << endl;
	delete fileMask;
}

char* BeadCounter::getFile() const
{
	return curFilePath;
}

bool BeadCounter::findNextFile()
{
	if (_findnext(hFile, &file) == 0)
	{
		curFilePath = file.name;
		counter++;
		return true;
	}
	else if (hSubdir != -1L)
	{// there are subdirectories
		_findclose(hFile);
		_chdir("..");

		if (_findnext(hSubdir, &subdir) == 0)
		{
			string subdirName = subdir.name;
			string curPath = ".\\" + subdirName;
			_chdir(curPath.c_str());
			
			if ((hFile = _findfirst(fileMask, &file)) == -1L)
			{
				cerr << "No " << fileMask << " files in current directory! " << curPath << endl;
				exit(1);
			}
			else
			{
				curFilePath = file.name;
				counter++;
				return true;
			}
		}
		else
		{
			_findclose(hSubdir);
			return false;
		}
	}
	else //_hSubdir == -1L
	{// there're no subdirectories
		_findclose(hFile);
	}
	return false;
}
#endif
