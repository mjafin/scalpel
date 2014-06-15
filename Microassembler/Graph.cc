#include "Graph.hh"

/******************************************************************
** Graph.cc
**
** Class for representing and storing a bi-directed de Bruijn graph  
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

// clear
//////////////////////////////////////////////////////////////

void Graph_t::clear(bool flag)
{
	if(flag) {
		readid2info.clear();
	}
	totalreadbp_m = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		delete mi->second;
	}

	nodes_m.clear();

	source_m = NULL;
	sink_m = NULL;

	if (ref_m && flag == true)
	{
		delete ref_m;
		ref_m = NULL;
	}
}

// loadSequence
//////////////////////////////////////////////////////////////

void Graph_t::loadSequence(int readid, const string & seq, bool isRef, int trim5)
{	
	if (!isRef)
	{
		totalreadbp_m += seq.length();
	}

	CanonicalMer_t uc;
	CanonicalMer_t vc;

	set<Mer_t> readmers;

	int end = seq.length() - K;
	int offset = 0;
	for (; offset < end; offset++)
	{
		uc.set(seq.substr(offset,   K));
		vc.set(seq.substr(offset+1, K));

		//cerr << readid << "\t" << offset << "\t" << uc << "\t" << vc << endl;

		MerTable_t::iterator ui = nodes_m.find(uc.mer_m);
		MerTable_t::iterator vi = nodes_m.find(vc.mer_m);

		if (ui == nodes_m.end())
		{
			ui = nodes_m.insert(make_pair(uc.mer_m, new Node_t(uc.mer_m))).first; 
		}

		if (vi == nodes_m.end())
		{
			vi = nodes_m.insert(make_pair(vc.mer_m, new Node_t(vc.mer_m))).first; 
		}

		//ui->second->appendRefFlag(isRef);
		//vi->second->appendRefFlag(isRef);

		if (!isRef)
		{
			if (offset == 0) 
			{ 
				ui->second->cov_m += 1; 
				ui->second->updateCovDistr((int)(ui->second->cov_m));
				if (uc.ori_m == F)
				{
					ui->second->addReadStart(readid, 0, trim5, uc.ori_m);
				}
				else
				{
					ui->second->addReadStart(readid, K-1, trim5, uc.ori_m);
				}
			}
			
			vi->second->cov_m += 1;
			vi->second->updateCovDistr((int)(vi->second->cov_m));
		}

		Edgedir_t fdir, rdir;

		if      (uc.ori_m == F && vc.ori_m == F) { fdir = FF; rdir = RR; }
		else if (uc.ori_m == F && vc.ori_m == R) { fdir = FR; rdir = FR; }
		else if (uc.ori_m == R && vc.ori_m == F) { fdir = RF; rdir = RF; }
		else if (uc.ori_m == R && vc.ori_m == R) { fdir = RR; rdir = FF; }

		readmers.insert(uc.mer_m);

		if (readmers.find(vc.mer_m) != readmers.end())
		{
			if (VERBOSE) { cerr << "cycle detected in read " << readid << " offset: " << offset << " : " << seq << endl; }

			if (readid > -1)
			{
				if (readCycles == 0)
				{
					cerr << "WARNING: Cycles detected in the reads" << endl << endl;
				}

				readCycles++;
			}

			//readid = -1;
		}

		ui->second->addEdge(vc.mer_m, fdir, readid);
		vi->second->addEdge(uc.mer_m, rdir, readid);
	}
}


// addSeqCov
//////////////////////////////////////////////////////////////

void Graph_t::addSeqCov(int readid, const string & seq, int cov)
{
 // cerr << "AddSeqCov: " << seq << " " << cov << endl;

	CanonicalMer_t uc;
	CanonicalMer_t vc;

	int end = seq.length() - K;
	for (int offset = 0; offset < end; offset++)
	{
		uc.set(seq.substr(offset,   K));
		vc.set(seq.substr(offset+1, K));

		//cerr << readid << "\t" << offset << "\t" << uc << "\t" << vc << endl;

		MerTable_t::iterator ui = nodes_m.find(uc.mer_m);
		MerTable_t::iterator vi = nodes_m.find(vc.mer_m);

		if (ui == nodes_m.end())
		{
			ui = nodes_m.insert(make_pair(uc.mer_m, new Node_t(uc.mer_m))).first; 
		}

		if (vi == nodes_m.end())
		{
			vi = nodes_m.insert(make_pair(vc.mer_m, new Node_t(vc.mer_m))).first; 
		}

	    ui->second->cov_m += cov;
		vi->second->cov_m += cov;

		Edgedir_t fdir, rdir;

		if      (uc.ori_m == F && vc.ori_m == F) { fdir = FF; rdir = RR; }
		else if (uc.ori_m == F && vc.ori_m == R) { fdir = FR; rdir = FR; }
		else if (uc.ori_m == R && vc.ori_m == F) { fdir = RF; rdir = RF; }
		else if (uc.ori_m == R && vc.ori_m == R) { fdir = RR; rdir = FF; }

		ui->second->addEdge(vc.mer_m, fdir, readid);
		vi->second->addEdge(uc.mer_m, rdir, readid);
	}
}

// trim
//////////////////////////////////////////////////////////////

void Graph_t::trim(int readid, const string & seq, const string & qv, bool isRef)
{
	int trim3 = 0;
	int trim5 = 0;
	int len = seq.length();
	string cseq = seq;

	for (int i = 0; i < len; i++) { cseq[i] = toupper(cseq[i]); }

	while ((!isDNA(seq[trim5]) || (qv[trim5] < MIN_QUAL)) && (trim5 < len)) { trim5++; }

	if (trim5 < len) {
		while ((!isDNA(seq[len-1-trim3]) || (qv[len-1-trim3] < MIN_QUAL)) && (trim3 < len)) { trim3++; }

		readid2info[readid].isjunk = false;
		
		for (int i = trim5; i < len-trim3; i++)
		{
			if (!isDNA(seq[i])) {
				// skip the junk
				readid2info[readid].isjunk = true;
				break;
			}
		}
	}
	else { readid2info[readid].isjunk = true; }
	
	readid2info[readid].trm5 = trim5;
	readid2info[readid].trm3 = trim3;
}

// trimAndLoad
//////////////////////////////////////////////////////////////

void Graph_t::trimAndLoad(int readid, const string & seq, const string & qv, bool isRef)
{
	int len = seq.length();
	string cseq = seq;

	for (int i = 0; i < len; i++) { cseq[i] = toupper(cseq[i]); }

	int trim5 = 0;
	while ((!isDNA(seq[trim5]) || (qv[trim5] < MIN_QUAL)) && (trim5 < len)) { trim5++; }

	if (trim5 < len)
	{
		int trim3 = 0;
		while ((!isDNA(seq[len-1-trim3]) || (qv[len-1-trim3] < MIN_QUAL)) && (trim3 < len)) { trim3++; }

		bool cleanRead = true;

		for (int i = trim5; i < len-trim3; i++)
		{
			if (!isDNA(seq[i]))
			{
				// skip the junk
				cleanRead = false;
				break;
			}
		}

		if (cleanRead)
		{
			if (trim5 || trim3) { cseq = seq.substr(trim5, len-trim5-trim3); }
			loadSequence(readid, cseq, isRef, trim5);
		}
	}
}

int Graph_t::countBastardReads()
{
	int bastards = 0;

	for (unsigned int i = 0; i < readid2info.size(); i++)
	{
		if (readid2info[i].code_m == CODE_BASTARD)
		{
			bastards++;
		}
	}

	return bastards;
}

int Graph_t::countMappedReads()
{
	int mapped = 0;

	for (unsigned int i = 0; i < readid2info.size(); i++)
	{
		if (readid2info[i].code_m == CODE_MAPPED)
		{
			mapped++;
		}
	}

	return mapped;
}

// addRead
////////////////////////////////////////////////////////////////

ReadId_t Graph_t::addRead(const string & set, const string & readname, const string & seq, char code)
{
	ReadId_t retval = readid2info.size();
	readid2info.push_back(ReadInfo_t(set, readname, seq, code));
	return retval;
}

// addMates
////////////////////////////////////////////////////////////////

void Graph_t::addMates(ReadId_t r1, ReadId_t r2)
{
	readid2info[r1].mateid_m = r2;
	readid2info[r2].mateid_m = r1;
}

// printReads
////////////////////////////////////////////////////////////////

void Graph_t::printReads()
{
	for (unsigned int i = 0; i < readid2info.size(); i++)
	{
		cout << i << "\t" << readid2info[i].readname_m << "\t" << readid2info[i].set_m << endl;
	}
}

// addPair
//////////////////////////////////////////////////////////////

void Graph_t::addPair(const string & set,
	const string & readname,
	const string & seq1, const string & qv1,
	const string & seq2, const string & qv2,
	char code)
{
	if (INCLUDE_BASTARDS || (code == CODE_MAPPED))
	{
		int rd1 = addRead(set, readname+"_1", seq1, code);
		trimAndLoad(rd1, seq1, qv1, false);

		int rd2 = addRead(set, readname+"_2", seq2, code);
		trimAndLoad(rd2, seq2, qv2, false);

		addMates(rd1, rd2);
	}
}

// addUnpaired
//////////////////////////////////////////////////////////////

void Graph_t::addUnpaired(const string & set,
	const string & readname,
	const string & seq,
	const string & qv,
	char code)
{
	int rd = addRead(set, readname, seq, code);
	//trimAndLoad(rd, seq, qv, false);
	trim(rd, seq, qv, false);
	
}

// addpaired
//////////////////////////////////////////////////////////////

void Graph_t::addpaired(const string & set,
	const string & readname,
	const string & seq,
	const string & qv,
	const int mate_id,
	char code)
{
	int rd;
	if(mate_id == 1) {
		rd = addRead(set, readname+"_1", seq, code);
		trim(rd, seq, qv, false);
	}
	else if(mate_id == 2) {
		rd = addRead(set, readname+"_2", seq, code);
		trim(rd, seq, qv, false);
	}
	else { // no mate
		rd = addRead(set, readname+"_0", seq, code);
		trim(rd, seq, qv, false);
	}
}


// build graph by trimming and loading the reads
//////////////////////////////////////////////////////////////

void Graph_t::buildgraph()
{
	for (unsigned int i = 0; i < readid2info.size(); i++)
	{
		if ( !(readid2info[i].isjunk) ) { // skip junk (not A,C,G,T)
			string seq = readid2info[i].seq_m;
			int len = seq.length();
			int t5 = readid2info[i].trm5;
			int t3 = readid2info[i].trm3;
			string cseq = seq;
			if (t5 || t3) { cseq = seq.substr(t5, len-t5-t3); }
			loadSequence(i, cseq, false, t5);
		}
	}	
}

// loadReadsSFA
//////////////////////////////////////////////////////////////

void Graph_t::loadReadsSFA(const string & filename)
{
	cerr << "LoadReadsSFA " << filename << endl;

	FILE * fp = xfopen(filename, "r");

	char rbuffer[BUFFER_SIZE];
	char sbuffer[BUFFER_SIZE];
	char qbuffer[BUFFER_SIZE];

	int readid;

	while (fscanf(fp, "%s\t%s", rbuffer, sbuffer) == 2)
	{
		readid = addRead("rd", rbuffer, sbuffer, 'L');

		int l = strlen(sbuffer);
		for (int i = 0; i < l; i++) { qbuffer[i] = MIN_QUAL; } qbuffer[l] = '\0';

		trimAndLoad(readid, sbuffer, qbuffer, false);
	}

	xfclose(fp);
}



// dfs to detect cycles
//////////////////////////////////////////////////////////////
bool Graph_t::hasCycle() {
	
	//cout << "Check for cycles (kmer = " << K << ")..." << endl;
	
	bool ans1 = false; 
	bool ans2 = false; 
	bool ans = false;
	
	if ( (source_m != NULL) && (sink_m != NULL) ) {
		
		MerTable_t::iterator mi;
	
		for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++) {
			Node_t * node = mi->second;

			if (node->isRef_m)    { continue; }
			//if (node->touchRef_m) { continue; }
		
			node->setColor(WHITE);
		}
		
		hasCycleRec(source_m, F, &ans1);
		hasCycleRec(source_m, R, &ans2);
		ans = ans1 || ans2;
	}
	
	if(ans) {
		cout << "Cycle found in the graph (kmer = " << K << ")!" << endl;
	}
	
	return ans;
}

void Graph_t::hasCycleRec(Node_t * node, Ori_t dir, bool *ans) {
		
	if(node != NULL && !(*ans)) {
		
		node->setColor(GREY);
			
		for (unsigned int i = 0; i < node->edges_m.size(); i++) {
						
			Edge_t & edge = node->edges_m[i];
			if (edge.isDir(dir)) {

				Node_t * other = getNode(edge);
				
				if(other->isRef_m) { continue; }
				
				//if (DFS_VERBOSE) { cerr << "     ==> " << other->nodeid_m << endl; }

				if(other->getColor() == GREY) { // cycle!
					//cout << "cycle found!" << endl;
					*ans = true;
					break;
				}
				if(other->getColor() == WHITE) { 
					hasCycleRec(other,  edge.destdir(), ans);
				}
			}
		}
		node->setColor(BLACK);
	}
}

// dfs
//////////////////////////////////////////////////////////////

bool Graph_t::findRepeatsInGraphPaths(Node_t * source, Node_t * sink, Ori_t dir)
{
	//cerr << "Check for perfect and near perfect repeats in the graph" << endl;
	
	// exit if source or sink are not defined
	if ( (source_m == NULL) || (sink_m == NULL) ) { return false; }

	bool answer = false;
	int complete = 0;
	int deadend = 0;
	int visit = 0;

	deque<Path_t *> Q;

	Path_t * start = new Path_t(K);
	start->nodes_m.push_back(source);
	start->dir_m = dir;
	start->len_m = K;

	Q.push_back(start);
	
	while (!Q.empty())
	{
		visit++;

		if ((DFS_LIMIT) && (visit > DFS_LIMIT))
		{
			//cerr << "WARNING: DFS_LIMIT (" << DFS_LIMIT << ") exceeded" << endl;
			break;
		}

		Path_t * path = Q.front();
		Q.pop_front();

		Node_t * cur = path->curNode();

		if (cur == sink)
		{
			// success!
			complete++;
			
			if(isAlmostRepeat(path->str(), K, MAX_MISMATCH)) {
				answer = true;
				cerr << "Near-perfect repeat in assembled sequence for kmer " << K << endl;
				break;
			}
		}
		else {
			int tried = 0;

			for (unsigned int i = 0; i < cur->edges_m.size(); i++)
			{
				Edge_t & edge = cur->edges_m[i];
				if (edge.isDir(path->dir_m))
				{
					tried++;

					Node_t * other = getNode(edge);

					Path_t * newpath = new Path_t(path,K);

					newpath->nodes_m.push_back(other);
					newpath->edgedir_m.push_back(edge.dir_m);
					newpath->dir_m = edge.destdir();
					newpath->len_m = path->len_m + other->strlen() - K + 1;

					Q.push_back(newpath);
				}
			}

			if (tried == 0)
			{
				deadend++;
				//cerr << "deadend: " <<  cur->nodeid_m << endl;
			}
		}

		delete path;
	}

	while (!Q.empty())
	{
		Path_t * path = Q.front();
		delete path;
		Q.pop_front();
	}
	
	return answer;
}



// dfs
//////////////////////////////////////////////////////////////

void Graph_t::dfs(Node_t * source, Node_t * sink, Ori_t dir, 
	Ref_t * ref, FILE * fp, bool printPathsToFile)
{
	cerr << endl << "searching from " << source->nodeid_m << " to " << sink->nodeid_m << " dir: " << dir << endl;

	const string & refseq = ref->seq;

	int complete = 0;
	int toolong = 0;
	int deadend = 0;
	int visit = 0;
	int shortpaths = 0;
	int allcycles = 0;

	int perfect   = 0;
	int withsnps  = 0;
	int withindel = 0;
	int withmix   = 0;

	int reflen = refseq.length();

	deque<Path_t *> Q;

	Path_t * start = new Path_t(K);
	start->nodes_m.push_back(source);
	start->dir_m = dir;
	start->len_m = K;

	Q.push_back(start);

	bool DFS_VERBOSE = 0;
	bool OLD_VERBOSE = VERBOSE;
	
	while (!Q.empty())
	{
		visit++;

		if ((DFS_LIMIT) && (visit > DFS_LIMIT))
		{
			cerr << "WARNING: DFS_LIMIT (" << DFS_LIMIT << ") exceeded" << endl;
			break;
		}

		Path_t * path = Q.front();
		Q.pop_front();

		Node_t * cur = path->curNode();

		if (DFS_VERBOSE) { cerr << " --> " << cur->nodeid_m << " [" << Q.size() << "]" << endl; }

		if (cur == sink)
		{
			// success!
			complete++;
			
			//if(isAlmostRepeat(path->str(), K, MAX_MISMATCH)) {
			//	cerr << "Near-perfect repeat in assembled sequence for kmer " << K << "... skip!" << endl;
			//	continue;
			//}
			
			path->match_bp = 0;
			path->snp_bp = 0;
			path->ins_bp = 0;
			path->del_bp = 0;

			if (VERBOSE) { cerr << "alignment" << endl; }
			if (VERBOSE) { cerr << "r:  " << refseq << endl; }
			if (VERBOSE) { cerr << "p:  " << path->str() << endl; }

			// Get alignment
			string ref_aln;
			string path_aln;
			//string cov_aln;
			
			global_align_aff(refseq, path->str(), ref_aln, path_aln, 0, 0);
			//global_cov_align_aff(refseq, path->str(), path->covDistr(), ref_aln, path_aln, cov_aln, 0, 0);

			assert(ref_aln.length() == path_aln.length());

			cerr << "r':" << ref_aln << endl;  
			cerr << "p':" << path_aln << " " << path->cov() << " [" << path->mincov() << " - " << path->maxcov() << "]" << endl; 
			cerr << "d':"; 
			for (unsigned int i = 0; i < ref_aln.length(); i++)
			{
				if (ref_aln[i] == path_aln[i]) { path->match_bp++;  cerr << ' '; }
				else if (ref_aln[i] == '-')    { path->ins_bp++;    cerr << '^'; }
				else if (path_aln[i] == '-')   { path->del_bp++;    cerr << 'v'; }
				else                           { path->snp_bp++;    cerr << 'x'; }
			}
			cerr << "\n";
			//cerr << "c':" << cov_aln << endl;
			
			//print coverage distribution along the sequence path
			cerr << "c':" << path->covstr() << endl;
			vector<int> coverage = path->covDistr();
			
			//print coverage distribution along the edges of the path
			/*
			coverage = path->readCovNodes();
			num = new char[30];
			cerr << "e':";
			for (unsigned int i = 0; i < coverage.size(); i++) {
				sprintf(num, "%.1f", coverage[i]);
				cerr << num << " ";
			}
			cerr << endl;
			*/
			
			try {
				// scan aligned sequences for differences

				unsigned int pos_in_ref = 0; 
				unsigned int refpos  = 0; 
				unsigned int pathpos = 0; 

				Node_t * spanner;

				char code;

				vector<Transcript_t> transcript;
			
				// cov_window keeps track of the minimum coverage in a window of size K
				multiset<int> cov_window;
			 			
				int end = min( (int)(K-1), (int)(coverage.size()-1) );
				assert(end >= 0);
				assert(end < (int)coverage.size());
				for (int t=0; t<end; t++) { cov_window.insert(coverage[t]); }
			
				for (unsigned int i = 0; i < ref_aln.length(); i++) {	
					int toadd = min( (int)(pathpos+K-1), (int)(coverage.size()-1));
					assert(toadd >= 0);
					assert(toadd < (int)coverage.size());
					cov_window.insert(coverage[toadd]);
							
					unsigned int old_pathpos = pathpos;
					if (ref_aln[i] == '-') {
						code = '^'; // insertion
						pos_in_ref = refpos; // save value of position in reference before increment
						pathpos++;           
					}
					else if (path_aln[i] == '-') { 
						code = 'v'; // deletion
						pos_in_ref = refpos; // save value of position in reference before increment
						refpos++;            
					}
					else { 
						code = '=';
						if (ref_aln[i] != path_aln[i]) { code = 'x'; }
						pos_in_ref = refpos; // save value of position in reference before increment
						refpos++; 
						pathpos++; 
					}

					// pathpos is a 1-based coordinate, need to substruc 1 to correctly search 
					// into the sequence path which is indexed using 0-based coordinates
					spanner = path->pathcontig(pathpos);
					spanner->setRead2InfoList(&readid2info);

					// get min coverage in the window
					int cov_at_pos = *(cov_window.begin()); // need to be adjusted for possible left normalization after alignment				
				
					//float cov_at_pos = path->covAt(pathpos+(K/2)); // need to be adjusted for possible left normalization after alignment
					//float cov_at_pos = path->covAt(pathpos); 
				
					if(old_pathpos != pathpos) { // remove old coverage only if change in path position
						assert(pathpos > 0);
						assert(pathpos <= coverage.size());
						multiset<int>:: iterator it = cov_window.find(coverage[pathpos-1]);
						cov_window.erase(it);
					}
					//multiset<float>:: iterator it = cov_window.find(path->covAt(i));
					//cov_window.erase(path->covAt(i));
				
					if (code != '=')
					{ 
						cerr << (ref_aln[i] == '-' ? '*' : ref_aln[i]) << " " << (path_aln[i] == '-' ? '*' : path_aln[i]) << " " << code 
							<< " " << pos_in_ref + ref->refstart + ref->trim5 << " " << pathpos 
							<< " " << spanner->nodeid_m << " " << spanner->cov_m << " " << cov_at_pos << " " << spanner->reads_m.size() << " " << spanner->cntReadCode(CODE_BASTARD)
							<< endl;

						unsigned int rrpos = pos_in_ref+ref->refstart+ref->trim5;
						unsigned int ts = transcript.size();

						//if (ts > 0)
						//{
						//  cerr << "== " << ts << " "
						//       << transcript[ts-1].code << " "
						//       << transcript[ts-1].pos << " " 
						//       << transcript[ts-1].ref << " "
						//       << transcript[ts-1].qry << endl;
						//}
						
						// compute previous base to the event for both reference and alternative sequences
						// [required for VCF output format]
						int pr=i-1; // referecne index
						assert(pr >= 0);
						int pa=i-1; // alternative index
						assert(pa >= 0);
						while( (ref_aln[pr] != 'A') && (ref_aln[pr] != 'C') && (ref_aln[pr]) != 'G' && (ref_aln[pr] != 'T') ) { pr--; }
						while( (path_aln[pa] != 'A') && (path_aln[pa] != 'C') && (path_aln[pa] != 'G') && (path_aln[pa] != 'T') ) { pa--; }
						

						if (((code == '^') && (ts > 0) && (transcript[ts-1].code == code) && (transcript[ts-1].pos == rrpos)) ||
							((code == 'v') && (ts > 0) && (transcript[ts-1].code == code) && ((transcript[ts-1].pos + transcript[ts-1].ref.length()) == rrpos)))
						{
							// extend the indel
							transcript[ts-1].ref += ref_aln[i];
							transcript[ts-1].qry += path_aln[i];
							(transcript[ts-1].cov_distr).push_back(cov_at_pos);
						
							//if(code == '^') { // insertion 
							//	if(cov_at_pos < min_cov) { min_cov = cov_at_pos; } 
							//}
						}
						else
						{
							//transcript.push_back(Transcript_t(rrpos, code, ref_aln[i], path_aln[i], spanner->cov_m));
							transcript.push_back(Transcript_t(rrpos, code, ref_aln[i], path_aln[i], cov_at_pos, ref_aln[pr], path_aln[pa]));
							//min_cov = 1000000;
						}
					}
				}
		
				cerr << ">p_" << ref->refchr << ":" << ref->refstart << "-" << ref->refend << "_" << complete
					<< " cycle: " << path->hasCycle_m
					<< " match: " << path->match_bp
					<< " snp: "   << path->snp_bp
					<< " ins: "   << path->ins_bp
					<< " del: "   << path->del_bp;

				for (unsigned int ti = 0; ti < transcript.size(); ti++)
				{
					//cerr << " " << transcript[ti].pos << ":" << transcript[ti].ref << "|" << transcript[ti].qry << "|" << transcript[ti].cov;
					cerr << " " << transcript[ti].pos << ":" << transcript[ti].ref << "|" << transcript[ti].qry << "|" << transcript[ti].getAvgCov() << "|" << transcript[ti].getMinCov() << "|" << transcript[ti].prev_bp_ref << "|" << transcript[ti].prev_bp_alt;
				}

				cerr << endl;

				if      ((path->snp_bp + path->ins_bp + path->del_bp) == 0) { perfect++;   }
				else if ((path->snp_bp) == 0)                               { withindel++; }
				else if ((path->ins_bp + path->del_bp) == 0)                 { withsnps++;  }
				else                                                        { withmix++;   }

				if(printPathsToFile) {
					fprintf(fp,  ">p_%s:%d-%d_%d len=%d cov=%0.2f mincov=%0.2f maxcov=%0.2f pathlen=%d hasCycle=%d match=%d snp=%d ins=%d del=%d pathstr=%s\n",
						ref->refchr.c_str(), ref->refstart, ref->refend, complete, 
						path->strlen(), path->cov(), path->mincov(), path->maxcov(), path->pathlen(), 
						path->hasCycle_m, path->match_bp, path->snp_bp, path->ins_bp, path->del_bp, path->pathstr().c_str());
					
					fprintf(fp, "%s\n", path->str().c_str());
				}

				for (unsigned int i = 0; i < path->nodes_m.size(); i++)
				{
					Node_t * cur = path->nodes_m[i];
					cur->onRefPath_m++;
				}

				if ((PATH_LIMIT) && (complete > PATH_LIMIT))
				{
					cerr << "WARNING: PATH_LIMIT reached: " << PATH_LIMIT << endl;
					break;
				}
			
			}
			catch(std::out_of_range& e) {
		   		cerr << "An exception occurred: " << e.what( ) << endl;
		 	}
			catch (...) { 
				cout << "default exception" << endl; 
			}
			
		}
		else if (path->len_m > reflen + MAX_INDEL_LEN)
		{
			// abort
			toolong++;
			//cerr << "too long: " << path->pathstr() << " " << path->str() << endl;
		}
		else
		{
			int tried = 0;

			for (unsigned int i = 0; i < cur->edges_m.size(); i++)
			{
				Edge_t & edge = cur->edges_m[i];
				if (edge.isDir(path->dir_m))
				{
					tried++;

					Node_t * other = getNode(edge);
					if (DFS_VERBOSE) { cerr << "     ==> " << other->nodeid_m << endl; }

					if (!path->hasCycle_m && path->hasCycle(other))
					{
						allcycles++;
					}

					Path_t * newpath = new Path_t(path,K);

					newpath->nodes_m.push_back(other);
					newpath->edgedir_m.push_back(edge.dir_m);
					newpath->dir_m = edge.destdir();
					newpath->len_m = path->len_m + other->strlen() - K + 1;

					Q.push_back(newpath);
				}
			}

			if (tried == 0)
			{
				deadend++;
				//cerr << "deadend: " <<  cur->nodeid_m << endl;
			}
		}

		delete path;
	}

	while (!Q.empty())
	{
		Path_t * path = Q.front();
		delete path;
		Q.pop_front();
	}

	if (complete == 0)
	{
		// didn't find an end-to-end path

		if (visit == 2)
		{
			// source to single node
			assert(source->edges_m.size() == 1);

			complete++;
			shortpaths++;

			Node_t * node = getNode(source->edges_m[0]);
			node->onRefPath_m++;

			string str = node->str_m;

			if (source->edges_m[0].destdir() == R)
			{
				str = CanonicalMer_t::rc(str);
			}


			VERBOSE = 1;

			if (VERBOSE) { cerr << "partial align" << endl; }

			if (VERBOSE) { cerr << "alignment" << endl; }
			if (VERBOSE) { cerr << "r:  " << refseq << endl; }
			if (VERBOSE) { cerr << "p:  " << str << endl; }

			// Get alignment
			string ref_aln;
			string path_aln;
			string cov_aln;

			global_align_aff(refseq, str, ref_aln, path_aln, 1, 0);

			if (VERBOSE) 
			{ 
				cerr << "r': " << ref_aln << endl;  
				cerr << "p': " << path_aln << endl; 
			}

			VERBOSE = OLD_VERBOSE;

			int match_bp = 0;
			int snp_bp = 0;
			int ins_bp = 0;
			int del_bp = 0;

			assert(ref_aln.length() == path_aln.length());

			for (unsigned int i = 0; i < ref_aln.length(); i++)
			{
				if (ref_aln[i] == path_aln[i]) { match_bp++; }
				else if (ref_aln[i] == '-')    { ins_bp++; }
				else if (path_aln[i] == '-')   { del_bp++; }
				else                           { snp_bp++; }
			}

			cerr << ">sp" << complete
				<< " match: " << match_bp
				<< " snp: " << snp_bp
				<< " ins: " << ins_bp
				<< " del: " << del_bp
				<< endl;

			if      ((snp_bp + ins_bp + del_bp) == 0) { perfect++;   }
			else if ((snp_bp) == 0)                   { withindel++; }
			else if ((ins_bp + del_bp) == 0)          { withsnps++;  }
			else                                      { withmix++;   }

			if(printPathsToFile) {
				fprintf(fp,  ">p_%d len=%d cov=%0.2f mincov=%0.2f maxcov=%0.2f pathlen=%d match=%d snp=%d ins=%d del=%d pathstr=%s\n",
					complete, (int) str.length(), node->cov_m, node->cov_m, node->cov_m, 1, match_bp, snp_bp, ins_bp, del_bp, "shortmatch");
				
				fprintf(fp, "%s\n", str.c_str());
			}
		}
	}

	VERBOSE = OLD_VERBOSE;

	int withmixindel = withmix + withindel;
	int withmixsnp   = withmix + withsnps;
	int withvar      = withsnps + withindel + withmix;

	cerr << " refcomp: "    << ref_m->refcomp
		<< " refnodes: "   << ref_m->refnodes-2
		<< " visit: "      << visit
		<< " complete: "   << complete
		<< " allcycles: "  << allcycles
		<< " shortpaths: " << shortpaths
		<< " toolong: "    << toolong
		<< " deadend: "    << deadend << endl;

	cerr << " perfect: "      << perfect
		<< " withsnps: "     << withsnps
		<< " withindel: "    << withindel
		<< " withmix: "      << withmix 
		<< " withmixindel: " << withmixindel
		<< endl;

	if(printPathsToFile) {
		fprintf(fp, ">stats\treflen=%d\tnumreads=%d\tcov=%0.02f\ttrim5=%d\ttrim3=%d\tnodes=%d\trefnodes=%d\tcomp=%d\trefcomp=%d\tvisit=%d\tcomplete=%d\tallcycles=%d\tshortpath=%d\ttoolong=%d\tdeadend=%d\tperfect=%d\twithsnps=%d\twithindel=%d\twithmix=%d\twithmixindel=%d\twithmixsnp=%d\twithvar=%d\n",
			(int) ref_m->seq.length(), (int) readid2info.size(), ((float)totalreadbp_m / (float) ref_m->seq.length()),
			ref_m->trim5, ref_m->trim3, (int) nodes_m.size()-2, ref_m->refnodes-2, ref_m->allcomp, ref_m->refcomp,
			visit, complete, allcycles, shortpaths, toolong, deadend, perfect, withsnps, withindel, withmix, withmixindel, withmixsnp, withvar);
	}
}

string Graph_t::nodeColor(Node_t * cur, string & who)
{
	string color = COLOR_ALL;

	if      (cur->nodeid_m == "source") { color = COLOR_SOURCE; }
	else if (cur->nodeid_m == "sink")   { color = COLOR_SINK; }
	else if (cur->touchRef_m)           { color = COLOR_TOUCH; }

	stringstream whostr;

	map<string, int> whocnt;

	set<ReadId_t>::iterator si;
	for (si = cur->reads_m.begin(); si != cur->reads_m.end(); si++)
	{
		whocnt[readid2info[*si].set_m]++;
	}

	bool isChl = false;
	bool isParent = false;

	map<string, int>::iterator mi;
	for (mi = whocnt.begin(); mi != whocnt.end(); mi++)
	{
		if (mi != whocnt.begin()) { whostr << " "; }
		whostr << mi->first << ":" << mi->second;

		if (mi->second >= COV_THRESHOLD)
		{
			if ( (mi->first == "self") || (mi->first == "sib") )
			{
				isChl = true;
			}
			else
			{
				isParent = true;
			}
		}
	}


	if ((color == COLOR_ALL) && (cur->cov_m < COV_THRESHOLD))
	{
		color = COLOR_LOW;
	}
	else if (isChl && !isParent)
	{
		color = COLOR_NOVO;
	}

	who = whostr.str();

	return color;
}

string Graph_t::edgeColor(Node_t * cur, Edge_t & e)
{
	string color = COLOR_ALL;

	string w;
	string c1 = nodeColor(cur, w);
	string c2 = nodeColor(getNode(e), w);

	if ((c1 == COLOR_LOW) || (c2 == COLOR_LOW))
	{
		color = COLOR_LOW;
	}
	else if ((c1 == COLOR_NOVO) && (c2 == COLOR_NOVO))
	{
		color = COLOR_NOVO;
	}

	return color;
}



// printDot
//////////////////////////////////////////////////////////////
void Graph_t::printDot(const string & filename)
{
	cerr << "Saving graph: " << filename << endl;

	FILE * fp = xfopen(filename, "w");

	if (PRINT_DOT_READS)
	{
		for (unsigned int i = 0; i < readid2info.size(); i++)
		{
			fprintf(fp, "// %s %d %s -> %d (%s)\n",
				readid2info[i].set_m.c_str(),
				i, 
				readid2info[i].readname_m.c_str(),
				readid2info[i].mateid_m, 
				readid2info[i].contigid_m.c_str());
		}

		fprintf(fp, "\n\n");
	}
	
	fprintf(fp, "digraph structs{\n");

    int printstrings = 1;

    if (printstrings)
    {
	  fprintf(fp, "  graph [bgcolor=white,size=\"8.5,11\",ratio=fill,center=true]\n");
	  fprintf(fp, "  node [shape=record];\n");
	  fprintf(fp, "  rankdir=LR\n");
    }
    else
    {
	  fprintf(fp, "  graph [bgcolor=black,layout=neato,rankdir=LR]\n");
	  fprintf(fp, "  node [shape=circle,style=filled,fontsize=1,fixedsize=true,hight=1,width=1];\n");
	  fprintf(fp, "  edge [fixedsize=true,len=1.2];\n");
    }

	int nodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		nodes++;

		Node_t * cur = mi->second;
		cur->setRead2InfoList(&readid2info);

		string who;
		string color = nodeColor(cur, who);

		if (NODE_STRLEN == 0)
		{
			fprintf(fp, "  %s [label=\". | <F> . | <R> .\" color=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				color.c_str());
		}
		else if (cur->strlen() > NODE_STRLEN)
		{
			fprintf(fp, "  %s [label=\"%d:%s | <F> %s... | <R> len=%d cov=%0.02f rd:%d B:%d %s\" color=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				nodes,
				cur->nodeid_m.c_str(),
				cur->str_m.substr(0, NODE_STRLEN).c_str(),
				cur->strlen(),
				cur->cov_m,
				(int) cur->reads_m.size(),
				cur->cntReadCode(CODE_BASTARD),
				who.c_str(),
				color.c_str());
		}
		else
		{
			fprintf(fp, "  %s [label=\"%d:%s | <F> %s | <R> len=%d cov=%0.02f rd:%d B:%d %s\" color=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				nodes,
				cur->nodeid_m.c_str(),
				cur->str_m.substr(0, NODE_STRLEN).c_str(),
				cur->strlen(),
				cur->cov_m,
				(int) cur->reads_m.size(),
				cur->cntReadCode(CODE_BASTARD),
				who.c_str(),
				color.c_str());
		}

		if (PRINT_DOT_READS)
		{
			fprintf(fp, "  //reads:");

			set<ReadId_t>::iterator ri;

			for (ri = cur->reads_m.begin(); ri != cur->reads_m.end(); ri++)
			{
				//fprintf(fp, " %s", readid2info[e.readids_m[j]].readname_m.c_str());
				fprintf(fp, " %d", *ri);
			}

			fprintf(fp, "\n");

			fprintf(fp, "  //readstarts:");

			for (unsigned int i = 0; i < cur->readstarts_m.size(); i++)
			{
				fprintf(fp, " %d:%c%d",
					cur->readstarts_m[i].readid_m,
					(cur->readstarts_m[i].ori_m == F) ? '+' : '-',
					cur->readstarts_m[i].nodeoffset_m);
			}

			fprintf(fp, "\n");

			fprintf(fp, "  //links:");

			ContigLinkMap_t::iterator li;
			for (li = cur->contiglinks_m.begin();
			li != cur->contiglinks_m.end();
			li++)
			{
				fprintf(fp, " %s(%d)", li->first.c_str(), li->second->linkCnt());
			}

			fprintf(fp, "\n");
		}

		for (unsigned int i = 0; i < cur->edges_m.size(); i++)
		{
			Edge_t & e = cur->edges_m[i];

			bool printed = 0;

			if (e.dir_m != RR)
			{
				if ((e.dir_m == FF) || (cur->nodeid_m <= e.nodeid_m))
				{
					string ecolor = edgeColor(cur, e);

					printed = 1;
					fprintf(fp, "    %s:%c -> %s:%c [arrowhead=\"normal\" arrowtail=\"normal\" color=\"%s\"]\n",
						cur->nodeid_m.c_str(), Edge_t::toString(e.startdir()),
						e.nodeid_m.c_str(),    Edge_t::toString(e.destdir()),
						ecolor.c_str());
				}
			}

			if (!printed)
			{
				fprintf(fp, "    //%s:%c -> %s:%c\n",
					cur->nodeid_m.c_str(), Edge_t::toString(e.startdir()),
					e.nodeid_m.c_str(),    Edge_t::toString(e.destdir()));
			}


			if (PRINT_DOT_READS)
			{
				fprintf(fp, "    //reads:");

				for (unsigned int j = 0; j < e.readids_m.size(); j++)
				{
					//fprintf(fp, " %s", readid2info[e.readids_m[j]].readname_m.c_str());
					fprintf(fp, " %d", e.readids_m[j]);
				}

				fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n\n");
	}

	fprintf(fp, "}\n");

	xfclose(fp);
}


// printFasta: print all contigs
//////////////////////////////////////////////////////////////

void Graph_t::printFasta(const string & filename)
{
	cerr << "Saving fasta: " << filename << endl;

	FILE * fp = xfopen(filename, "w");

	int nodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		nodes++;
		Node_t * cur = mi->second;

        int fdeg = 0;
        int rdeg = 0;

        for (unsigned int i = 0; i < cur->edges_m.size(); i++)
        {
            Edge_t & edge = cur->edges_m[i];
            if (edge.isDir(F)) { fdeg++; }
            else               { rdeg++; }
        }

		fprintf(fp, ">%d:%s len=%d cov=%0.2f fdeg=%d rdeg=%d\n", nodes, cur->nodeid_m.c_str(), cur->strlen(), cur->cov_m, fdeg, rdeg);
		fprintf(fp, "%s\n",  cur->str_m.c_str());
	}

	xfclose(fp);
}

// printPairs: print all pairs of neighboring contigs
//////////////////////////////////////////////////////////////

void Graph_t::printPairs(const string & filename)
{
	cerr << "Saving pairs fasta: " << filename << endl;

	FILE * fp = xfopen(filename, "w");

	int nodes = 0;

	Path_t pairpath(K);

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		if (cur->edges_m.size() == 0)
		{
		// print isolated contigs

			nodes++;
			fprintf(fp, ">%d:%s len=%d cov=%0.2f\n", nodes, cur->nodeid_m.c_str(), cur->strlen(), cur->cov_m);
			fprintf(fp, "%s\n",  cur->str_m.c_str());
		}
		else
		{
			// print contig pairs, making sure a given pair is only printed once

			for (unsigned int i = 0; i < cur->edges_m.size(); i++)
			{
				Edge_t & e = cur->edges_m[i];

				if (e.dir_m != RR) 
				{
					if ((e.dir_m == FF) || (cur->nodeid_m < e.nodeid_m))
					{
						nodes++;

						Node_t * other = getNode(e);

						pairpath.reset();
						pairpath.nodes_m.push_back(cur);
						pairpath.edgedir_m.push_back(e.dir_m);
						pairpath.nodes_m.push_back(other);

						string pathstr = pairpath.pathstr();
						string path    = pairpath.str();
						float  cov     = pairpath.cov();

						int    strlen  = path.length();

						fprintf(fp, ">%d:%s len=%d cov=%0.2f\n", nodes, pathstr.c_str(), strlen, cov);
						fprintf(fp, "%s\n",  path.c_str());
					}
				}
			}
		}
	}

	xfclose(fp);
}


// markRefEnds
//////////////////////////////////////////////////////////////

void Graph_t::markRefEnds(Ref_t * refinfo)
{
	if (VERBOSE) 
	{ 
		cerr << "Looking at " << refinfo->hdr << " " << refinfo->seq.length() 
			<< " " << refinfo->seq << endl; 
	}

	ref_m = refinfo;
	int refid; 
	if(!is_ref_added) {
		refid = addRead("ref", ref_m->hdr, ref_m->rawseq, 'R');
		is_ref_added = true;
		if (VERBOSE) { cerr << "refid: " << refid << endl; }
	}

	//loadSequence(refid, ref_m->seq, 1, 0);

	CanonicalMer_t source_mer;
	CanonicalMer_t sink_mer;
	CanonicalMer_t tmp_mer;
	
	Node_t * source_tmp;
	Node_t * sink_tmp;

	int source_offset = -1;
	int sink_offset = -1;
	int offset = -1;

	ref_m->trim5 = -1;
	ref_m->trim3 = -1;

	source_m = NULL;
	sink_m = NULL;
	
	bool ambiguous_match;
	
	// Find the first matching mer with sufficient coverage
	source_tmp = NULL;
	ambiguous_match = false;
	for (offset = 0; offset < (int) ref_m->rawseq.length(); offset++)
	{	
		tmp_mer.set(ref_m->rawseq.substr(offset, K));
		source_tmp = getNode(tmp_mer);

		if ((source_tmp) && (source_tmp->cov_m >= COV_THRESHOLD))
		{ 
			if(source_m == NULL) { // found 1st match
				source_m = source_tmp;
				source_mer = tmp_mer;
				source_offset = offset;
			}
			else { // check if there is another match for the same mer
				if (source_m == source_tmp) { //found another identical match
					source_m = NULL;
					ambiguous_match = true;
					break;
				}
			}
			//source_m = source_tmp;
			//break; 
		}
	}

	if(ambiguous_match) {
		cerr << "Ambiguous match to reference for source" << endl;
		return;
	}
	
	if (!source_m) {
		cerr << "No match to reference for source" << endl;
		return;
	}

	// Find the last matching mer
	sink_tmp = NULL;
	ambiguous_match = false;
	for (offset = ref_m->rawseq.length()-K; offset >= 0; offset--)
	{
		tmp_mer.set(ref_m->rawseq.substr(offset, K));
		sink_tmp = getNode(tmp_mer);

		if ((sink_tmp) && (sink_tmp->cov_m >= COV_THRESHOLD))
		{
			if(sink_m == NULL) { // found 1st macth
				sink_m = sink_tmp;
				sink_mer = tmp_mer;
				sink_offset = offset;
			}
			else { // check if there is another match for the same mer
				if (sink_m == sink_tmp) { //found another identical match
					sink_m = NULL;
					ambiguous_match = true;
					break;
				}
			}
			//sink_m = sink_tmp;
			//break; 
		}
	}
	
	if(ambiguous_match) {
		cerr << "Ambiguous match to reference for sink" << endl;
		return;
	}
	
	if (!sink_m) {
		cerr << "No match to reference for sink" << endl;
		return;
	}

	int ref_dist = sink_offset - source_offset + K;
	sink_offset = ref_m->rawseq.length() - sink_offset - K;
	ref_m->seq = ref_m->rawseq.substr(source_offset, ref_dist);

	cerr << "ref trim5: "     << source_offset 
		<< " trim3: "     << sink_offset
		<< " uncovered: " << source_offset + sink_offset
		<< " ref_dist: "  << ref_dist << endl;

	ref_m->trim5 = source_offset;
	ref_m->trim3 = sink_offset;

	//cerr << " searching " << source_m->nodeid_m << " to " << sink_m->nodeid_m << endl;

	// Add the fake source node
	Node_t * newsource = new Node_t("source");

	Edgedir_t sourcedir = FF;
	if (source_mer.ori_m == R) { sourcedir = FR; }

	bool CLIP_REF_ENDS = true;
	if (CLIP_REF_ENDS)
	{
		if (VERBOSE) { cerr << "checking for other source edges" << endl; }
		for (int i = source_m->edges_m.size()-1; i >= 0; i--)
		{
			if (Edge_t::edgedir_start(source_m->edges_m[i].dir_m) == Edge_t::flipdir(source_mer.ori_m))
			{
				Node_t * other = getNode(source_m->edges_m[i]);
				
				if ( (other != NULL) && (other != source_m) ) { 

					if (VERBOSE) { cerr << "  removing node before source: " << other->nodeid_m << endl; }

					other->removeEdge(source_m->nodeid_m, Edge_t::fliplink(source_m->edges_m[i].dir_m));
					source_m->edges_m.erase(source_m->edges_m.begin() + i);
				}
			}
		}
	}

	newsource->addEdge(source_mer.mer_m, sourcedir, refid);
	newsource->isRef_m = true;
	source_m->addEdge(newsource->nodeid_m, Edge_t::fliplink(sourcedir), refid);
	source_m = newsource;

	nodes_m.insert(make_pair(newsource->nodeid_m, newsource));

	// Add the fake sink node
	Node_t * newsink = new Node_t("sink");

	Edgedir_t sinkdir = RR;
	if (sink_mer.ori_m == R) { sinkdir = FF; } 

	if (CLIP_REF_ENDS) 
	{
		if (VERBOSE) { cerr << "checking for other sink edges" << endl; }
		for (int i = sink_m->edges_m.size()-1; i >= 0; i--)
		{
			if (Edge_t::edgedir_start(sink_m->edges_m[i].dir_m) == sink_mer.ori_m)
			{
				Node_t * other = getNode(sink_m->edges_m[i]);
				
				if ( (other != NULL) && (other != sink_m) ) { 
				
					if (VERBOSE) { cerr << "  removing node after sink: " << other->nodeid_m << endl; }

					other->removeEdge(sink_m->nodeid_m, Edge_t::fliplink(sink_m->edges_m[i].dir_m));
					sink_m->edges_m.erase(sink_m->edges_m.begin() + i);
				}
			}
		}
	}

	newsink->addEdge(sink_mer.mer_m, sinkdir, refid);
	newsink->isRef_m = true;
	sink_m->addEdge(newsink->nodeid_m, Edge_t::fliplink(sinkdir), refid);
	sink_m = newsink;

	nodes_m.insert(make_pair(newsink->nodeid_m, newsink));
}

// markRefNodes
//////////////////////////////////////////////////////////////

void Graph_t::markRefNodes()
{
	cerr << endl << "mark refnodes" << endl;
	int nodes = 0;
	int refnodes = 0;

	ref_m->refcompids.clear();

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		nodes++;
		refnodes += mi->second->markRef(ref_m, K);
		mi->second->component_m = 0;
	}

	int comp = 0;
	int refcomp = 0;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		if (mi->second->component_m != 0) { continue; }

		comp++;

		deque<Node_t *> Q;
		Q.push_back(mi->second);

		int touches = 0;

		while (!Q.empty())
		{
			Node_t * cur = Q.front();
			Q.pop_front();

			if (cur->component_m == 0)
			{
				cur->component_m = comp;

				if (cur->touchRef_m) { touches++; }

				for (unsigned int i = 0; i < cur->edges_m.size(); i++)
				{
					Node_t * next = getNode(cur->edges_m[i]);
					Q.push_back(next);
				}
			}
			else
			{
				if (cur->component_m != comp)
				{
					cerr << "ERROR: cur->comp: " << cur->component_m << " != comp: " << comp << endl;
					exit(1);
				}
			}
		}

		if (touches)
		{
			refcomp++;
			ref_m->refcompids.insert(comp);
		}
	}


	ref_m->refnodes = refnodes;
	ref_m->refcomp  = refcomp;
	ref_m->allcomp  = comp;

	cerr << " nodes: "    << nodes
		<< " refnodes: " << refnodes
		<< " comp: "     << comp
		<< " refcomp: "  << refcomp
		<< " refcompids: ";

	set<int>::iterator si;
	for (si = ref_m->refcompids.begin(); si != ref_m->refcompids.end(); si++)
	{
		cerr << " " << *si;
	}
	cerr << endl;
}

// denovoNodes
//////////////////////////////////////////////////////////////

void Graph_t::denovoNodes(const string & filename, const string & refname)
{
	cerr << endl << "finding denovo nodes" << endl;

	markRefNodes();

	FILE * fp = xfopen(filename, "w");

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		map<string, int> who;

		set<ReadId_t>::iterator si;
		for (si = cur->reads_m.begin(); si != cur->reads_m.end(); si++)
		{
			string & set = readid2info[*si].set_m;
			if (set != "ref")
			{
				who[set]++;
			}
		}

		if (who.size() == 1)
		{
			map<string, int>::iterator mi = who.begin();

			//char isRef = ref_m->isRefComp(cur->component_m) ? 'R' : 'N';

			if (VERBOSE)
			{
				fprintf(stderr, ">%s_%s len=%d cov=%0.02f comp=%d\n%s\n",
					mi->first.c_str(), 
					cur->nodeid_m.c_str(),
					cur->strlen(), 
					cur->cov_m,
					cur->component_m,
					cur->str_m.c_str());

			}

			fprintf(fp, ">%s_%s len=%d cov=%0.02f\n%s\n",
				mi->first.c_str(), 
				cur->nodeid_m.c_str(),
				cur->strlen(), 
				cur->cov_m,
				cur->str_m.c_str());
		}
	}

	xfclose(fp);
}

// alignRefNodes
//////////////////////////////////////////////////////////////

void Graph_t::alignRefNodes()
{
	int refpathnodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		if (cur->onRefPath_m)
		{
			refpathnodes++;
		}
	}

	cerr << " Found " << refpathnodes << " on ref path" << endl;
}

// countRefPath
//////////////////////////////////////////////////////////////

void Graph_t::countRefPath(const string & filename, const string & refname, bool printPathsToFile)
{
	if (source_m)
	{
		FILE * fp = NULL;
		
		if(printPathsToFile) {
			FILE * fp = xfopen(filename, "w");
		
			fprintf(fp, ">ref_raw\n%s\n",
				ref_m->rawseq.c_str());

			fprintf(fp, ">ref_trim %s trim5:%d trim3:%d\n%s\n", 
				refname.c_str(), ref_m->trim5, ref_m->trim3, ref_m->seq.c_str());
		}

		if (source_m != NULL && sink_m != NULL) {
			dfs(source_m, sink_m, F, ref_m, fp, printPathsToFile);
		}
		alignRefNodes();

		if(printPathsToFile) { xfclose(fp); }
		
	}
}



// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(Mer_t nodeid)
{
	MerTable_t::iterator ni = nodes_m.find(nodeid);

	if (ni == nodes_m.end()) { return false; }

	return ni->second;
}

// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(CanonicalMer_t mer)
{
	return getNode(mer.mer_m);
}


// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(Edge_t & edge)
{
	MerTable_t::iterator ni = nodes_m.find(edge.nodeid_m);

	if (ni == nodes_m.end()) { return false; }

	return ni->second;
}


// compressNode
//////////////////////////////////////////////////////////////

void Graph_t::compressNode(Node_t * node, Ori_t dir)
{
	bool cnVERBOSE = 0;

	if (cnVERBOSE) { cerr << "Compress " << node->nodeid_m << ":" << dir << endl; }

	while (true)
	{
		if (cnVERBOSE) { cerr << node << endl; }

		int uniqueid = node->getBuddy(dir);
		if (uniqueid == -1) { return; }
		if (node->isTandem()) { return; }

		// make sure they are mutual buddies
		Edgedir_t edir = node->edges_m[uniqueid].dir_m;

		if (cnVERBOSE) { cerr << " --> " << node->edges_m[uniqueid] << endl; }

		Ori_t bdir = F;

		if (edir == FF || edir == RF) { bdir = R; }

		Node_t * buddy = getNode(node->edges_m[uniqueid]);
		if (buddy->isTandem()) { return; }

		if (!buddy)
		{
			cerr << "couldn't get " << node->nodeid_m << " - buddyid: " << node->edges_m[uniqueid] << endl;
		}

		assert(buddy);
		assert(!buddy->dead_m);

		if (cnVERBOSE)
		{
			cerr << "    " << buddy << endl;
		}

		int buniqueid = buddy->getBuddy(bdir);

		if (buniqueid == -1) 
		{ 
			if (cnVERBOSE) { cerr << "node's buddy's buddy is not me" << endl; }
			return; 
		}

		assert(buddy->edges_m[buniqueid].nodeid_m == node->nodeid_m);

		if (cnVERBOSE) { cerr << " --> " << node->edges_m[uniqueid] << endl; }

		// str
		string astr = node->str_m;
		if (dir == R)
		{ 
			astr = CanonicalMer_t::rc(astr); 
			node->revreads();
			node->revCovDistr();
		}

		string bstr = buddy->str_m;
		if (Edge_t::edgedir_dest(edir) == R)
		{ 
			bstr = CanonicalMer_t::rc(bstr); 
			buddy->revreads();
			buddy->revCovDistr();
		}

		if (cnVERBOSE)
		{
			string shift;
			shift.append(astr.length()-K+1, ' ');
			cerr << "     " << astr << endl; 
			cerr << "     " << shift << bstr << endl;
		}


		assert(astr.substr(astr.length()-K+1, K-1) == bstr.substr(0, K-1));

		string mstr = astr + bstr.substr(K-1);

		if (cnVERBOSE) { cerr << "     " << mstr << endl; }

		if (dir == R) 
		{ 
			mstr = CanonicalMer_t::rc(mstr); 
		}

		node->str_m = mstr;

		// coverage
		int   amerlen = astr.length() - K + 1;
		int   bmerlen = bstr.length() - K + 1;
		float ncov = node->cov_m;
		float ccov = buddy->cov_m;
		
		// update the coverage of the overlapping region
		//for (unsigned int j = (node->cov_distr.size())-K+1; j < node->cov_distr.size(); j++) {
		//	node->cov_distr[j] += 1;
		//}
		// add coverage info for the new base-pairs 	
		for (unsigned int j = (K-1); j < buddy->cov_distr.size(); j++) {
			node->cov_distr.push_back(buddy->cov_distr[j]);
		}
		node->cov_m = ((ncov * amerlen) + (ccov * bmerlen)) / (amerlen + bmerlen);

		// reads
		set<ReadId_t>::iterator ri;
		for (ri = buddy->reads_m.begin(); ri != buddy->reads_m.end(); ri++)
		{
			node->reads_m.insert(*ri);
		}

		// add buddy read starts
		int shift = amerlen;
		for (unsigned int i = 0; i < buddy->readstarts_m.size(); i++)
		{
			ReadStart_t & rs = buddy->readstarts_m[i];
			node->readstarts_m.push_back(ReadStart_t(rs.readid_m, rs.nodeoffset_m+shift, rs.trim5_m, rs.ori_m));
		}

		// flip read starts
		if (dir == R)
		{
			node->revreads();
			node->revCovDistr();
		}

		node->sortReadStarts();

		// dead flag
		buddy->dead_m = true;

		// isRef
		node->isRef_m |= buddy->isRef_m;

		// node edges
		node->edges_m.erase(node->edges_m.begin()+uniqueid);

		// now move over the buddy edges
		for (int i = 0; i < (int) buddy->edges_m.size(); i++)
		{
			if (i == buniqueid) { continue; }

			Edge_t ne(buddy->edges_m[i]);

			if (edir == FR || edir == RF)
			{
				ne.dir_m = Edge_t::flipme(ne.dir_m);
			}

			Node_t * other = getNode(ne);

			if (cnVERBOSE) { cerr << "Keeping: " << ne << endl; }

			// NOTE: this comprassion can generate the error: 
			// scalpel: Path.cc:66: std::string Path_t::str(): Assertion `retval.substr(retval.length()-K+1) == nstr.substr(0, K-1)' failed.  
			if (other == buddy) 
			{
				cerr << "circle to buddy" << endl;
				ne.nodeid_m = node->nodeid_m;
				node->edges_m.push_back(ne);
			}
			else
			{
				node->edges_m.push_back(ne);
				other->updateEdge(buddy->nodeid_m, Edge_t::fliplink(buddy->edges_m[i].dir_m),
					node->nodeid_m, Edge_t::fliplink(ne.dir_m));
			}
		}
	}
}


// compress
//////////////////////////////////////////////////////////////

void Graph_t::compress()
{
	cerr << "compressing graph:";

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		if (mi->second->dead_m)  { continue; }
		if (mi->second->isRef_m) { continue; }

		compressNode(mi->second, F);
		compressNode(mi->second, R);
	}


	cleanDead();
}

// cleanDead
///////////////////////////////////////////////////////////////

void Graph_t::cleanDead()
{
	set<Mer_t> deadnodes;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		if (mi->second->dead_m) 
		{
			deadnodes.insert(mi->second->nodeid_m);
		}
	}

	cerr << "  removing " << deadnodes.size() << " dead nodes" << endl;

	set<Mer_t>::iterator di;
	for (di = deadnodes.begin(); di != deadnodes.end(); di++)
	{
		mi = nodes_m.find(*di);
		assert(mi != nodes_m.end());
		assert(mi->second->dead_m);

		delete mi->second;
		nodes_m.erase(mi);
	}
}


// removeNode
//////////////////////////////////////////////////////////////

void Graph_t::removeNode(Node_t * node)
{
	assert(node);
	assert(!node->dead_m);

	node->dead_m = true;

	for(unsigned int i = 0; i < node->edges_m.size(); i++)
	{
		Node_t * nn = getNode(node->edges_m[i]);

		if ((nn) && (nn != node))
		{
			nn->removeEdge(node->nodeid_m, Edge_t::fliplink(node->edges_m[i].dir_m));
		}
	}
}


// removeLowCov
//////////////////////////////////////////////////////////////

void Graph_t::removeLowCov()
{
	cerr << endl << "removing low coverage:";

	int lowcovnodes = 0;
	
	double avgcov = ((double) totalreadbp_m) / ((double) ref_m->rawseq.length());
	//cerr << "avgcov: " << avgcov << endl;

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * node = mi->second;

		if (node->isRef_m)    { continue; }
		//if (node->touchRef_m) { continue; }

		//if (node->cov_m <= TIP_COV_THRESHOLD)
		//if (node->minCov() <= TIP_COV_THRESHOLD)
		if ( (node->minCov() <= TIP_COV_THRESHOLD) || (node->minCov() <= (MIN_COV_RATIO*avgcov)) )
		{
			lowcovnodes++;
			removeNode(node);
		}
	}

	cerr << " found " << lowcovnodes;

	cleanDead();
	compress();

	printStats();
}


// removeTips
//////////////////////////////////////////////////////////////

void Graph_t::removeTips()
{
	int tips = 0;
	int round = 0;

	do
	{
		round++;
		tips = 0;

		cerr << endl << "remove tips round: " << round;

		MerTable_t::iterator mi;

		for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
		{
			Node_t * cur = mi->second;

			if (cur->isRef_m) { continue; }

			int deg = cur->edges_m.size();
			int len = cur->strlen() - K + 1;

			if ((deg <= 1) && (len < MAX_TIP_LEN))
			{
				removeNode(cur);
				tips++;
			}
		}

		cerr << " removed: " << tips << endl;

		if (tips) { compress(); }
	}
	while (tips);

	printStats();
}

// greedyTrim
//////////////////////////////////////////////////////////////

class CovCmp
{
public:

  CovCmp(Graph_t * g) : _g(g) {}

  bool operator() (const string & a, const string & b)
  {
    MerTable_t::iterator ai = _g->nodes_m.find(a);
    MerTable_t::iterator bi = _g->nodes_m.find(b);

    return ai->second->cov_m > bi->second->cov_m;
  }
  
  Graph_t * _g;
};

void Graph_t::greedyTrim()
{
	int branches = 0;

	cerr << endl << "greedy trim" << endl;

	MerTable_t::iterator mi;

  vector<string> nodelist;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
    nodelist.push_back(mi->first);
  }

  CovCmp covcmp(this);
  sort (nodelist.begin(), nodelist.end(), covcmp);

  for (unsigned int i = 0; i < nodelist.size(); i++)
  {
    mi = nodes_m.find(nodelist[i]);

		Node_t * cur = mi->second;

    if (i < 10)
    {
      cerr << i << " " << nodelist[i] << " " << cur->cov_m << endl;
    }

    // figure out the best edges in the forward and reverse

    if (cur->edges_m.size() < 2) { continue; }

    Edge_t bestf(cur->edges_m[0]); float covf = -1; int degf = 0;
    Edge_t bestr(cur->edges_m[0]); float covr = -1; int degr = 0;

    for (unsigned int j = 0; j < cur->edges_m.size(); j++)
    {
      Edge_t & edge = cur->edges_m[j];
      Node_t * other = getNode(edge);

      if (edge.isDir(F)) { degf++; if (other->cov_m > covf) { bestf = edge; covf = other->cov_m; } }
      else               { degr++; if (other->cov_m > covr) { bestr = edge; covr = other->cov_m; } }
    }

    // prune away all the other edges

    if (degf > 1 || degr > 1)
    {
      branches++;

      for (unsigned int i = 0; i < cur->edges_m.size(); i++)
      {
        Edge_t & edge = cur->edges_m[i];
        Node_t * other = getNode(edge);

        bool removeEdge = true;

        if (edge.isDir(F)) { if (edge.nodeid_m == bestf.nodeid_m && edge.dir_m == bestf.dir_m) { removeEdge = false; } }
        else               { if (edge.nodeid_m == bestr.nodeid_m && edge.dir_m == bestr.dir_m) { removeEdge = false; } }

        if (removeEdge)
        {
          other->removeEdge(cur->nodeid_m, Edge_t::fliplink(edge.dir_m));
        }
      }

      cur->edges_m.clear();

      if (covf != -1) { cur->edges_m.push_back(bestf); }
      if (covr != -1) { cur->edges_m.push_back(bestr); }
    }
	}

	cerr << " removed: " << branches << endl;

	if (branches) { compress(); }

	printStats();
}



// threadReads
//////////////////////////////////////////////////////////////
/* 
**
**         1\         /3
**           \       /
**            ---N---
**           /       \
**         2/         \4
**
**  Check for individual reads that span 1-N-3, 1-N-4, 2-N-3, 2-N-4
**
*/

void Graph_t::threadReads()
{
    if (MIN_THREAD_READS == -1)
    {
      cerr << "Skipping threading reads" << endl;
      return;
    }

	MerTable_t::iterator mi;

	bool oldverbose = VERBOSE;

	//VERBOSE = true;

	int thread = 1;
	int threadround = 0;

	while (thread)
	{
		thread = 0;
		threadround++;

		cerr << endl << "threading round: " << threadround << " " << nodes_m.size() << " nodes" << endl;

		if (VERBOSE)
		{
			cerr << "================================================================" << endl;
			printGraph();

			char buffer [1024];
			sprintf(buffer, "thread_%d.dot", threadround);
			printDot(buffer);

			cerr << "================================================================" << endl;
		}

		for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
		{
			Node_t * cur = mi->second;

			if (cur->isRef_m) { continue; }

			//if ((cur->degree(F) > 1) || (cur->degree(R) > 1))
			if ((cur->degree(F) > 1) && (cur->degree(R) > 1))
			{
				unsigned int numedges = cur->edges_m.size();

				vector<unsigned int> threaded;
				threaded.resize(numedges);

				for (unsigned int ti = 0; ti < numedges; ti++)
				{
					threaded[ti] = 0;
				}

				if (VERBOSE)
				{
					cerr << endl;
					cerr << "=Checking: " << cur->nodeid_m << " " << cur->strlen() << "bp" << endl;
					cerr << "================================================================" << endl;

					for (unsigned int i = 0; i < numedges; i++)
					{
						cerr << "  " << i << ": " << cur->nodeid_m << ":" << cur->edges_m[i] << endl;
					}
				}

				if (cur->isTandem()) 
				{ 
					if (VERBOSE) { cerr << " tandem, skipping" << endl << endl; }

					continue; 
				}

				// Find reads that span across the node

				for (unsigned int e1i = 0; e1i < numedges; e1i++)
				{
					Edge_t & e1 = cur->edges_m[e1i];
					Node_t * n1 = getNode(e1);

					if (e1.startdir() != F) { continue; }

					set<ReadId_t> e1reads;
					for (unsigned int j = 0; j < e1.readids_m.size(); j++)
					{
						e1reads.insert(e1.readids_m[j]);
					}

					for (unsigned int e2i = 0; e2i < numedges; e2i++)
					{
						Edge_t & e2 = cur->edges_m[e2i];
						Node_t * n2 = getNode(e2);
						

						if (e2.startdir() != R) { continue; }

						// check if any reads span edges from e1 to e2
						//////////////////////////////////////////////

						set<ReadId_t> overlap;

						if (VERBOSE) cerr << e1i << ":" << e1.label() << " == " << e2i << ":" << e2.label() << " :";

						for (unsigned int j = 0; j < e2.readids_m.size(); j++)
						{
							int jj = e2.readids_m[j];

							if (e1reads.find(jj) != e1reads.end())
							{
								if (VERBOSE) cerr << " " << jj;
								overlap.insert(jj);
							}
						}

						if (VERBOSE) cerr << endl;


						// check if any mates of n1 span to n2
						//////////////////////////////////////

						if (VERBOSE) cerr << "mates: " << endl;

						set<ReadId_t> mateoverlap;
						set<ReadId_t>::iterator s1, s2;

						for (s1 =  n1->reads_m.begin();
						s1 != n1->reads_m.end();
						s1++)
						{
							s2 = n2->reads_m.find(readid2info[*s1].mateid_m);

							if (s2 != n2->reads_m.end())
							{
								// mates are in n1 and n2
								// TODO: check orientation and spacing
								if (VERBOSE) cerr << " " << *s1 << " " << *s2 << endl;
							}
						}

						// overlap has reads that span from e1 (in F dir) to e2 (in R dir)

						if ( (int)overlap.size() >= MIN_THREAD_READS)
						{
							threaded[e1i]++;
							threaded[e2i]++;
						}
					}
				}

				// see what fraction of nodes were completely threaded

				unsigned int threadcnt = 0;
				for (unsigned int ti = 0; ti < numedges; ti++)
				{
					if (threaded[ti] > 0)
					{
						threadcnt++;
					}
				}
				

				if (VERBOSE)
				{ 
					cerr << "threads through " << threadcnt << " of " << numedges << endl;
				}

				// if the node was completely resolved, go ahead and resolve it

				if (threadcnt == numedges)
				{
					if (VERBOSE) { cerr << endl << "All edges threaded, applying..." << endl; }

					int copy = 0;
					vector<Node_t*> newnodes;

					for (unsigned int e1i = 0; e1i < numedges; e1i++)
					{
						Edge_t & e1 = cur->edges_m[e1i];

						if (e1.startdir() != F) { continue; }

						set<ReadId_t> e1reads;
						for (unsigned int j = 0; j < e1.readids_m.size(); j++)
						{
							e1reads.insert(e1.readids_m[j]);
						}

						for (unsigned int e2i = 0; e2i < numedges; e2i++)
						{
							Edge_t & e2 = cur->edges_m[e2i];

							if (e2.startdir() != R) { continue; }

							set<ReadId_t> overlap;

							for (unsigned int j = 0; j < e2.readids_m.size(); j++)
							{
								int jj = e2.readids_m[j];

								if (e1reads.find(jj) != e1reads.end())
								{
									overlap.insert(jj);
								}
							}

							// overlap has reads that span from e1 (in F dir) to e2 (in R dir)

							if ((int)overlap.size() >= MIN_THREAD_READS)
							{
								copy++;

								char buffer [1024];
								sprintf(buffer, "%s_%d", cur->nodeid_m.c_str(), copy);

								if (VERBOSE)
								{
									cerr << "thread " << cur->nodeid_m << " " << cur->strlen() << "bp" << endl;
									cerr << "  1: " << e1 << endl
										<< "  2: " << e2 << endl
										<< "  r[" << overlap.size() << "]:";

									set<ReadId_t>::iterator si;
									for (si = overlap.begin(); si != overlap.end(); si++) { cerr << " " << *si; }
									cerr << endl;

									cerr << "  making " << buffer << endl << endl;
								}

								Node_t * copy = new Node_t(cur->nodeid_m);

								copy->nodeid_m = buffer;      // nodeid
								copy->str_m = cur->str_m;     // sequence
								copy->cov_m = overlap.size(); // coverage
								copy->cov_distr.resize(copy->str_m.size());
								copy->updateCovDistr(overlap.size()); // coverage distribution

								// reads and edges
								set<ReadId_t>::iterator si;
								for (si = overlap.begin(); si != overlap.end(); si++)
								{
									copy->addEdge(e1.nodeid_m, e1.dir_m, *si);
									copy->addEdge(e2.nodeid_m, e2.dir_m, *si);
								}

								// TODO: update read starts?

								copy->touchRef_m = cur->touchRef_m;

								newnodes.push_back(copy);
							}
						}
					}

					if (!newnodes.empty())
					{
						if (VERBOSE) 
						{ 
							cerr << "cleaning up: " << cur->nodeid_m << endl; 

							for (unsigned int e1i = 0; e1i < cur->edges_m.size(); e1i++)
							{
								Edge_t & e1 = cur->edges_m[e1i];
								cerr << "   " << e1 << endl;
							} 
						} 

						removeNode(cur);

						thread++;

						for (unsigned int j = 0; j < newnodes.size(); j++)
						{
							Node_t * nn = newnodes[j];
							nodes_m.insert(make_pair(nn->nodeid_m, nn));

							if (VERBOSE) { cerr << "  swapping in: " << nn->nodeid_m << endl; }

							for (unsigned int k = 0; k < nn->edges_m.size(); k++)
							{
								Edge_t & e = nn->edges_m[k];
								Node_t * other = getNode(e.nodeid_m);

								if (VERBOSE) { cerr << "    edge: " << e << endl; }

								for (unsigned int r = 0; r < e.readids_m.size(); r++)
								{
									other->addEdge(nn->nodeid_m, Edge_t::fliplink(e.dir_m), e.readids_m[r]);
								}

								if (VERBOSE) { cerr << "      " << other << endl; }
							}
						}

						if (VERBOSE) { cerr << endl << endl; }
					}
				}
				else
				{
					if (VERBOSE) { cerr << "Not all edges threaded, skipping" << endl << endl; }
				}
			}
		}

		cerr << "  round: " << threadround << " threaded: " << thread << endl;

		if (thread > 0)
		{
			cleanDead();
			compress();
		}
	}

	VERBOSE=oldverbose;

	printStats();
}

// checkReadStarts
//////////////////////////////////////////////////////////////

void Graph_t::checkReadStarts()
{
	cerr << "checking read starts.... ";

	int all = 0;
	int bad = 0;

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		for (unsigned int i = 0; i < cur->readstarts_m.size(); i++)
		{
			ReadStart_t & rstart = cur->readstarts_m[i];
			ReadId_t rid         = rstart.readid_m;
			ReadInfo_t & rinfo   = readid2info[rid];

			string ckmer;
			string rkmer = rinfo.seq_m.substr(rstart.trim5_m, K);

			all++;

			if (rstart.ori_m == R)
			{
				ckmer = cur->str_m.substr(rstart.nodeoffset_m-K+1, K);
				ckmer = CanonicalMer_t::rc(ckmer);
			}
			else
			{
				ckmer = cur->str_m.substr(rstart.nodeoffset_m, K);
			}

			if ((rkmer != ckmer)) // || VERBOSE)
			{
				cerr << "Checking " << rid << " " << rinfo.readname_m 
					<< " " << rstart.ori_m 
					<< " offset:" << rstart.nodeoffset_m 
					<< " trim5:" << rstart.trim5_m << endl;
				cerr << "  " << rkmer << endl;
				cerr << "  " << ckmer << endl;

				cur->print(cerr) << endl;

				if (rkmer != ckmer)
				{
					bad++;
					cerr << "mismatch: " << cur->str_m << endl;
				}
				else
				{
					cerr << "ok" << endl;
				}
			}
		}
	}

	cerr << " found " << bad << " bad starts out of " << all << endl;
}

// updateContigReadStarts
//////////////////////////////////////////////////////////////

void Graph_t::updateContigReadStarts()
{
	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		for (unsigned int i = 0; i < cur->readstarts_m.size(); i++)
		{
			ReadId_t rid = cur->readstarts_m[i].readid_m;

			readid2info[rid].contigid_m = cur->nodeid_m;
			readid2info[rid].readstartidx_m = i;
		}
	}
}

// bundleMates
//////////////////////////////////////////////////////////////

void Graph_t::bundleMates()
{
	MerTable_t::iterator mi;

	// compute contig links

	int links = 0;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		for (unsigned int i = 0; i < cur->readstarts_m.size(); i++)
		{
			ReadId_t mateid = readid2info[cur->readstarts_m[i].readid_m].mateid_m;

			if (mateid != -1)
			{
				Mer_t matecontig = readid2info[mateid].contigid_m;

				if (matecontig != "")
				{
					cur->addContigLink(matecontig, cur->readstarts_m[i].readid_m);
					links++;
				}
			}
		}
	}

	cerr << "created " << links << " contig links" << endl;


	// bundle consistent links

	for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++)
	{
		Node_t * cur = mi->second;

		if (cur->contiglinks_m.size() != 0)
		{
			// cur is linked to at least 1 other contig

			ContigLinkMap_t::iterator mi;

			for (mi =  cur->contiglinks_m.begin();
			mi != cur->contiglinks_m.end();
			mi++)
			{
				ContigLinkList_t * list = mi->second;

				Node_t * other = getNode(mi->first);

				if (0)
				{
					cerr << cur->nodeid_m << "(" << cur->strlen() <<  ") -- " 
						<< list->linkCnt() << " -- " << mi->first
						<< "(" << other->strlen() << ")" << endl;
				}

				ContigLinkList_t bundles [4];

				int internal = 0;
				int invalid = 0;

				int lastpos = -1;
				int lastmpos = -1;

				for (unsigned int i = 0; i < list->linkCnt(); i++)
				{
					ContigLink_t & link = list->linklist_m[i];

					ReadId_t rid         = link.rid_m;
					ReadInfo_t & rinfo   = readid2info[rid];
					ReadStart_t & rstart = cur->readstarts_m[rinfo.readstartidx_m];

					ReadId_t mid         = rinfo.mateid_m;
					ReadInfo_t & minfo   = readid2info[mid];
					ReadStart_t & mstart = other->readstarts_m[minfo.readstartidx_m];

					Edgedir_t linkdir;

					if      ((rstart.ori_m == F) && (mstart.ori_m == F)) { linkdir = FF; }
					else if ((rstart.ori_m == F) && (mstart.ori_m == R)) { linkdir = FR; }
					else if ((rstart.ori_m == R) && (mstart.ori_m == F)) { linkdir = RF; }
					else if ((rstart.ori_m == R) && (mstart.ori_m == R)) { linkdir = RR; }

					int linkdist = 0;

					if (cur == other)
					{
						internal++;

						if      (linkdir == FR) { linkdist = mstart.nodeoffset_m - rstart.nodeoffset_m; }
						else if (linkdir == RF) { linkdist = rstart.nodeoffset_m - mstart.nodeoffset_m; }
						else
						{
							invalid++;
						}
					}
					else
					{
						int adist = cur->strlen() - rstart.nodeoffset_m;
						if (rstart.ori_m == R) { adist = rstart.nodeoffset_m; }

						int bdist = other->strlen() - mstart.nodeoffset_m;
						if (mstart.ori_m == R) { bdist = mstart.nodeoffset_m; }

						linkdist = INSERT_SIZE - (adist + bdist);
					}

					int dup = 0;

					if ((rstart.nodeoffset_m == lastpos) && (mstart.nodeoffset_m == lastmpos))
					{
						dup = 1;
					}

					lastpos  = rstart.nodeoffset_m;
					lastmpos = mstart.nodeoffset_m;

					if (0)
					{
						int lo = linkdist - 2 * INSERT_STDEV;
						int hi = linkdist + 2 * INSERT_STDEV;

						cerr << Edge_t::toString(linkdir) << ":" << linkdist << "\t"
							<< lo << "\t" << hi << "\t"
							<< rinfo.code_m << "\t" << dup << "\t"
							<< rid << "\t"  << rinfo.readname_m <<  "\t" << rstart.nodeoffset_m << "\t" << rstart.ori_m << "\t"
							<< mid << "\t"  << minfo.readname_m <<  "\t" << mstart.nodeoffset_m << "\t" << mstart.ori_m << endl;
					}

					if (!dup)
					{
						bundles[linkdir].addLink(rid, linkdist);
					}
					else
					{
						bundles[linkdir].addDup();
					}
				}

				//cerr << "--" << endl;

				for (int i = 0; i < 4; i++)
				{
					int cnt = bundles[i].linkCnt();
					int dups = bundles[i].dupCnt();

					if (cnt > 0)
					{
						float mean  = bundles[i].mean();
						float stdev = bundles[i].stdev(mean);

						cerr << cur->nodeid_m << "(" << cur->strlen() <<  ") "
							<< Edge_t::toString((Edgedir_t) i) << " "
							<< mi->first << "(" << other->strlen() << ")\t" 
							<< mean << "\t+/-\t" << stdev << "\t"
							<< cnt << "\t(" << dups << ")" << endl;
					}
				}

				//cerr << endl;
			}
		}
	}
}


// scaffold
//////////////////////////////////////////////////////////////

void Graph_t::scaffoldContigs()
{
	cerr << endl << "== scaffolding ==" << endl;

	checkReadStarts();
	updateContigReadStarts();
	bundleMates();

	cerr << endl << "==" << endl << endl;
}

// printGraph
//////////////////////////////////////////////////////////////

void Graph_t::printGraph()
{
	MerTable_t::iterator gi;
	for (gi = nodes_m.begin(); gi != nodes_m.end(); gi++)
	{
		cout << gi->second << endl;
	}
}

// printStats
//////////////////////////////////////////////////////////////

void Graph_t::printStats()
{
	int edgecnt = 0;
	int span = 0;

	MerTable_t::iterator gi;
	for (gi = nodes_m.begin(); gi != nodes_m.end(); gi++)
	{
		edgecnt += gi->second->edges_m.size();
		span += gi->second->strlen();
	}

	cerr << "  nodes: " << nodes_m.size() 
		<< " edges: " << edgecnt
		<< " span: " << span << endl;
}
