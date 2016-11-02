// global functions
#ifndef BLASTCONVGLOB_H
#define BLASTCONVGLOB_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "bc.h"

using namespace std;

int split(const string str, const string sepStr, vector<string>* v);

void commandLineHelp();

string truncate(const string name, const int ID);

void makeOffsets(multimap<int, long>* offsets, bool reverse, fstream& inFile);

void makeOffsetsUniq(multimap<int, long>* offsets, fstream& inFile);

#endif
