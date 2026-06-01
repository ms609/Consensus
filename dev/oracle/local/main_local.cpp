// Patched copy of dev/reference/FDCT_new/main_new.cpp for the LOCAL oracle.
// ONLY CHANGE: the algorithm is taken from argv[1] (defaulting to "freq")
// instead of being hard-wired to "freq", while keeping the LOCAL_TEST
// inp.txt / oup.txt I/O so the freqdiff R bridge can be reused verbatim.
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

#include "nex_parser.h"

// #include "freqdiff.h"
#include "freqdiff2.h"
// #include "freqdiff2_backup.h"
#include "local_consensus.h"
#include "Tree.cpp"

using namespace std::chrono;

#define LOCAL_TEST

#define Now high_resolution_clock::now
#define duration_sec(start,end) (float)(duration_cast<milliseconds>(end-start).count() / 1000.)

using namespace std;

int main(int argc, char** argv) {

	#ifdef LOCAL_TEST

	string algo(argc > 1 ? argv[1] : "freq");   // <-- ONLY CHANGE

	ifstream fin("inp.txt");
	vector<Tree*> trees = parse_nex_short(fin);

	#else

	if (argc != 3) {
		cout << "Usage: exec [freq|minrs|minis] filename" << endl;
		cout << "More info to come" << endl;
		return 0;
	}

	string algo(argv[1]);
	ifstream fin(argv[2]);
	vector<Tree*> trees = parse_nex(fin);

	#endif



	auto st = Now();

	Tree* consensus;
	if (algo == "freq") {
		consensus = freqdiff(trees);
	} else if (algo == "minrlc_exact") {
		consensus = minRLC_exact(trees);
	} else if (algo == "minilc_exact") {
		consensus = minILC_exact(trees);
	} else if (algo == "aho-build") {
		consensus = ahoBuild(trees);
	} else {
		std::cerr << "Algorithm not supported" << std::endl;
		return 1;
	}

	auto ed = Now();
	cout << duration_sec(st, ed) << endl;

	if (consensus == NULL) {
		std::cout << "No valid consensus found." << std::endl;
	} else {
		freopen("oup.txt", "w", stdout);
		cout << consensus->to_newick() << endl;
		fclose(stdout);
	}

	delete consensus;

	return 0;
}
