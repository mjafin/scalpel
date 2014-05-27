#ifndef MICROASSEMBLER_HH
#define MICROASSEMBLER_HH 1

/******************************************************************
** Microassembler.hh
**
** Tool for localized assembly of genomic regions using
** the de Bruijn paradigm to detect genetic variants
** 
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <unistd.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include "api/BamReader.h"
#include "api/BamWriter.h"

#include "align.hh"
#include "util.hh"
#include "Mer.hh"
#include "Ref.hh"
#include "ReadInfo.hh"
#include "ReadStart.hh"
#include "Transcript.hh"
#include "Edge.hh"
#include "ContigLink.hh"
#include "Node.hh"
#include "Path.hh"
#include "Graph.hh"

using namespace std;
using namespace HASHMAP;
using namespace BamTools;

#define bvalue(value) ((value ? "true" : "false"))

#define FATHER 1
#define MOTHER 2
#define SELF 3
#define SIBLING 4

// constants
//////////////////////////////////////////////////////////////////////////

const char Graph_t::CODE_MAPPED = 'M';
const char Graph_t::CODE_BASTARD = 'B';

const string Graph_t::COLOR_ALL    = "darkorange3"; //"black";
const string Graph_t::COLOR_LOW    = "deepskyblue4"; //"grey85";
const string Graph_t::COLOR_NOVO   = "red";
const string Graph_t::COLOR_SOURCE = "green\" style=\"filled";
const string Graph_t::COLOR_SINK   = "yellow\" style=\"filled";
const string Graph_t::COLOR_TOUCH  = "pink";

class Microassembler {

public:
	
	int BUFFER_SIZE;
	
	// Configuration
	//////////////////////////////////////////////////////////////////////////

	bool VERBOSE;
	bool PRINT_DOT_READS;

	bool PRINT_RAW;
	bool PRINT_ALL;

	bool PRINT_DENOVO;
	bool PRINT_REFPATH;

	int MIN_QV;
	int QV_RANGE;
	int MIN_QUAL;
	int MIN_MAP_QUAL;

	string MAPFILE;
	string RG_FILE;
	string REFFILE;
	string READSET;

	string PREFIX;

	int minK;
	int maxK;
	int MAX_TIP_LEN;
	unsigned int MIN_THREAD_READS;
	int COV_THRESHOLD;
	double MIN_COV_RATIO;
	int TIP_COV_THRESHOLD;

	int INCLUDE_BASTARDS;

	bool SCAFFOLD_CONTIGS;
	int  INSERT_SIZE;
	int  INSERT_STDEV;

	int  QUAD_ASM;
    int  FASTQ_ASM;
	int  NODE_STRLEN;

	int  BAMFILE;

	int DFS_LIMIT;
	int PATH_LIMIT;
	int MAX_INDEL_LEN;
	int MAX_MISMATCH;
	
	// data structures
	//////////////////////////////////////////////////////////////////////////
	
	int graphCnt;
	set<string> readgroups;
	set<string> RG_father;
	set<string> RG_mother;
	set<string> RG_self;
	set<string> RG_sibling;
	
	map<string, Ref_t *> reftable; // table of references (exomes) to analyze
	
	Microassembler() { 
		graphCnt = 0;
		
		BUFFER_SIZE = 10*1024;

		VERBOSE         = false;
		PRINT_DOT_READS = true;

		PRINT_RAW       = true;
		PRINT_ALL       = false;

		PRINT_DENOVO    = false;
		PRINT_REFPATH   = false;

		MIN_QV         = 10;
		QV_RANGE       = '!';
		MIN_QUAL       = MIN_QV + QV_RANGE;
		MIN_MAP_QUAL   = 0;

		READSET = "qry";

		minK = 10;
		maxK = 100;
		MAX_TIP_LEN = minK;
		MIN_THREAD_READS = 3;
		COV_THRESHOLD = 5;
		MIN_COV_RATIO = 0.01;
		TIP_COV_THRESHOLD = 1;

		INCLUDE_BASTARDS = 0;

		SCAFFOLD_CONTIGS = 0;
		INSERT_SIZE = 150;
		INSERT_STDEV = 15;

		QUAD_ASM = 0;
        FASTQ_ASM = 0;
		NODE_STRLEN = 1000;

		BAMFILE = 0;
		RG_FILE = "";

		DFS_LIMIT = 1000;
		PATH_LIMIT = 0;
		MAX_INDEL_LEN = 250;
		MAX_MISMATCH = 1;
	}
		
	~Microassembler() { }
	
	void printConfiguration(ostream & out);	
	void loadRefs(const string & filename);
	void loadRG(const string & filename, int member);
	void processGraph(Graph_t & g, const string & refname, const string & prefix, int minK, int maxK);
	void fastqAsm(Graph_t & g, const string & prefix);
    void greedyAssemble(Graph_t & g, const string & prefix);
	int run(int argc, char** argv);
};

#endif
