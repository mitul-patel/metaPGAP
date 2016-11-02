#ifndef BC_H
#define BC_H

#include <cstdlib>

void myOpenDir(const char* dirpath);

#if defined(__GNUC__)
	#include <dirent.h>	
	#include <iostream>
	#include <sys/types.h>
	#include <string>

	#define INT64 long long
	#define UINT64 unsigned long long
	#define ATOI64 atoll
	#define SPRTR '/'
	#define FMASK ".tab"
	#define DEL "rm"
	#define COPY "cp"

	using namespace std;

	class BeadCounter
	{
	public:
		BeadCounter	(const char* dirPath, const char* dMask, const char* fMask);
		BeadCounter	(const char* dirPath, const char* fMask);
		~BeadCounter();
		const char* getFile() const;
		bool findNextFile();
	      
	private:
		int counter;
		string subdirName;
		string curFilePath;
		string fileMask;
		char startDir[PATH_MAX];
		struct dirent *dp;
		DIR *dir;
	};
//********************** Win32
	#elif defined(_MSC_VER)
		#include <io.h>
		#include <direct.h>
		#include <iostream>
		#include <fstream>
		#include <string>

		#define INT64 __int64
		#define UINT64 unsigned __int64
		#define ATOI64 _atoi64
		#define SPRTR '\\'
		#define FMASK "*.tab" 
		#define DEL "del"
		#define COPY "copy"

		using namespace std;

		class BeadCounter
		{
		public:
					BeadCounter	(const char* dirPath, const char* dMask, const char* fMask);
					BeadCounter	(const char* dirPath, const char* fMask);
					~BeadCounter();
			char*	getFile		() const;
			bool	findNextFile();

		private:
			int		counter;
			char*	curFilePath;
			char*	fileMask;
			char	startDir[_MAX_PATH];

			struct _finddata_t subdir;	//structure for subdirs in the working directory
			struct _finddata_t file;	//structure for files in the processing subdirectory	
			intptr_t hSubdir;			//handle for subdirectories;
			intptr_t hFile;				//handle for files
		};

			#else
			#error "Platform not supported. Need to update source code"
			#endif 

#endif		
