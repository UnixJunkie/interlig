/*
 * InterLig.cpp
 *
 *  Created on: Feb 10, 2016
 *      Author: Claudio Mirabello <claudio.mirabello@liu.se>
 */

#include "InterLig.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#define FASTRAND_MAX 32767.00

#include <stdio.h>  /* defines FILENAME_MAX */
#include <unistd.h>

using namespace std;

bool verbose = false;
bool superimpose = false;
bool map = false;

double best = -1;

int rounds = 1;
double chainMult = 10.0;

double percent = 0;
double T0 = 0.2;
double d0 = 0.5;
double epsilon = 0.5;
double dW = 0.5; //DEFAULT VALUES
long double pvalues[30][1001];
long double pvaluesMol2[36][1001];

std::string blosum;
std::string superfile;

int nulls;
RunningStat rs;

Molecule* molA;
Molecule* molB;
InterLig* compare;

static unsigned int g_seed;
int seed = 0;
//srand (time(NULL));

//Used to seed the generator.
inline void fast_srand(int seed) {

	g_seed = seed;

}

//fastrand routine returns one integer, similar output value range as C lib.
inline int fastrand() {

	g_seed = (214013 * g_seed + 2531011);

	return (g_seed >> 16) & 0x7FFF;

}

std::string GetCurrentWorkingDir(){

	//char buff[FILENAME_MAX];
	//GetCurrentDir(buff, FILENAME_MAX);
	//std::string current_working_dir(buff);

	std::string path = "";
	pid_t pid = getpid();
	char buf[20] = {0};
	sprintf(buf,"%d",pid);

	std::string _link = "/proc/";
	_link.append( buf );
	_link.append( "/exe");
	char buff[1024];

	int ch = readlink(_link.c_str(), buff, 1024);

	if (ch != -1) {
		buff[ch] = 0;
		path = buff;
		std::string::size_type t = path.find_last_of("/");
		path = path.substr(0, t);
	}

	return path;
}

void p_values_mol2() {

	std::ifstream file((GetCurrentWorkingDir() + "//p_values_mol2size_step1_0.6_0.25_0.75").c_str());

	int intl = 0;

	for (int row = 0; row < 36; row++) {
		std::string line;
		std::getline(file, line);
		if (!file.good())
			break;

		std::stringstream iss(line);

		for (int col = 0; col < 1002; ++col) {

			std::string val = "";
			std::getline(iss, val, ' ');
			if (!iss.good())
				break;
			std::stringstream convertor(val);

			if (col == 0) {
				convertor >> intl;
			} else {
				convertor >> pvaluesMol2[int(intl) - 5][col - 1];
			}
		}
	}
	return;
}

std::vector<Molecule*> read_list_of_mol2s(char* file){

	std::vector<Molecule*> molecules;

	ifstream fin(file); //to count stuff
	ifstream fin2(file);	//to read stuff

	string line;
	string molecule("@<TRIPOS>MOLECULE");

	int nMolecules = 0;

	getline(fin, line);

	while (fin.good()) {

		if (line.compare(0, molecule.length(), molecule) == 0)
			nMolecules++;

		getline(fin, line);
	}

	fin.clear();
	fin.seekg(0, std::ios::beg);


	for (int i = 0; i < nMolecules; i++) {

		Molecule* thisMol = new Molecule(&fin, &fin2, "mol");
		molecules.push_back(thisMol);

		if(i % 10000 == 0 && verbose){
			cout << "Loading molecules: " << i << "\n" << flush;
		}
	}

	return molecules;

}

std::vector<Molecule*> copy_list_of_mol2s(char* file){

	ifstream fin(file); //to count stuff
	ifstream fin2(file);	//to read stuff

	std::vector<Molecule*> moleculesCopy;

	string line;
	string molecule("@<TRIPOS>MOLECULE");

	int nMolecules = 0;

	getline(fin, line);

	while (fin.good()) {

		if (line.compare(0, molecule.length(), molecule) == 0)
			nMolecules++;

		getline(fin, line);
	}

	fin.clear();
	fin.seekg(0, std::ios::beg);


	for (int i = 0; i < nMolecules; i++) {

		Molecule* thisMolCopy = new Molecule(&fin, &fin2, "fullmol");
		moleculesCopy.push_back(thisMolCopy);

	}

	return moleculesCopy;

}

void parse_options(int argc, char* argv[]){

	for (int arg = 4; arg < argc; arg++){

		if (!strcmp(argv[arg], "-nullP")){

			percent = atoi(argv[arg + 1]);

		}else if (!strcmp(argv[arg], "-anneal")){

			rounds = atoi(argv[arg + 1]);
			if (verbose)
				cout << "Annealing rounds: " << rounds << "\n";

		} else if (!strcmp(argv[arg], "-d0")){

			d0 = atof(argv[arg + 1]);

			if (verbose)
				cout << "d0: " << d0 << "\n";

		}else if(!strcmp(argv[arg], "-v")){

			cout << "Verbose" << "\n" << flush;
			verbose = true;

		}else if(!strcmp(argv[arg], "-super")){

			cout << "Superimpose" << "\n" << flush;
			superfile = argv[arg + 1];
			superimpose = true;

		}else if(!strcmp(argv[arg], "-eps")){

			epsilon = atof(argv[arg + 1]);
			if (verbose)
				cout << "Epsilon: " << epsilon << "\n" << flush;

		}else if(!strcmp(argv[arg], "-matrix")){

			blosum = argv[arg + 1];

			if (verbose)
				cout << "Blosum matrix: " << blosum << "\n" << flush;

		}else if(!strcmp(argv[arg], "-dW")){

			dW = atof(argv[arg + 1]);

			if (verbose)
				cout << "Sequence weight: " << dW << "\n" << flush;

		}else if(!strcmp(argv[arg], "-seed")){

			seed = atoi(argv[arg + 1]);

			if (verbose)
				cout << "New seed: " << atoi(argv[arg + 1]) << "\n" << flush;

		}else if(!strcmp(argv[arg], "-ch")){

			chainMult = atof(argv[arg + 1]);
			if (verbose)
				cout << "New chainMult: " << atof(argv[arg + 1]) << "\n" << flush;

		}else if(!strcmp(argv[arg], "-T0")){

			T0 = atof(argv[arg + 1]);
			if (verbose)
				cout << "New initial Temperature: " << atof(argv[arg + 1]) << "\n" << flush;

		}else if(!strcmp(argv[arg], "-map")){

			map = true;
			if (verbose)
				cout << "Output mapping of residues\n" << flush;

		}

	}
}

double anneal() {

	//shuffle B, reset null correspondences
	//molB->shuffle(1000);
	molA->resetNull();

	//calculate the difference between the two (only where A and B overlap and A has no null correspondences)

	//rs.Clear();
	int r;
	int c;

	int loop = 1;

	double p = 0;
	double de = 0;


	double T = T0;

	double dice = 0;

	double mean = 0;
	double stdev = 0;

	double accRate = 1;
	int accepted = 1;
	int rejected = 0;

	int maxChainLength = max(400,
			int(chainMult * molB->getLength() * (molB->getLength() - 1))); //2*lengthMin*(lengthMin-1);
	int maxAccept = molB->getLength() * (molB->getLength() - 1);

	int toNull = -1;

	for (int i = 0; i < nulls; i++) {

		while (toNull < 0 || molA->getNull()[toNull]) {
			toNull = fastrand() % molA->getLength();
		}

		molA->setNull(toNull);
	}
	if (verbose)
		cout << "Caching difference..." << "\n" << flush;
	compare->cacheDifference();
	best = compare->getScore();

	if (verbose)
		cout << "Difference cached" << "\n" << flush;

	//a or b selection parameters
	char AorB = 'A';

	while (1) {

		if (nulls == 0 || loop % 3) {

			//swapping on B
			AorB = 'B';

			r = 0;
			c = 0;

			/*illegal moves are:				(for c, r in [0; lengthB), c > r)
			 1) swap two residues aligned to null correspondences in A (nullA[c] and nullA[r])
			 2) swap a residue in a null correspondence in A with a residue in the outer subset in B (nullA[r] and c >= lengthA)
			 3) swap two residues in the outer subset of B (r >= lengthA, c will be too then)
			 */
			while (r >= c || (molA->getNull()[c] and molA->getNull()[r])
					|| (molA->getNull()[r] and c >= molA->getLength())
					|| r >= molA->getLength()) {

				r = fastrand() % molB->getLength();
				c = fastrand() % molB->getLength();

			}

		} else {
			//null correspondences on A

			AorB = 'A';

			r = 0;
			c = 0;

			/*illegal moves are:				(for c, r in [0; lengthA), c > r)
			 1) swap two residues aligned to null correspondences in A (nullA[c] and nullA[r])
			 2) swap two non null residues
			 */
			while (r >= c || (molA->getNull()[c] and molA->getNull()[r])
					|| (!(molA->getNull()[c]) and !(molA->getNull()[r]))) {
				r = fastrand() % molA->getLength();
				c = fastrand() % molA->getLength();
			}

		}

		if (AorB == 'A') {

			compare->deltaA(r, c);
			compare->deltaSeqA(r, c);
		}

		if (AorB == 'B') {

			compare->deltaB(r, c);
			compare->deltaSeqB(r, c);

		}

		de = epsilon * -compare->getDeltaN()
										+ dW * (double) compare->getDeltaS();// / compare->getDeltaD();
		rs.Push(compare->getScore() + de);

		if (accepted > maxAccept || accepted + rejected > maxChainLength) {

			//starting a new markov chain
			mean = rs.Mean();
			stdev = rs.StandardDeviation();

			T = compare->temperatureSD(mean, stdev, T);

			accepted = 1;
			rejected = 0;

			rs.Clear();
		}

		if (de < 0) {

			p = 1;

		} else {

			p = exp(-de / T);

		}

		dice = (((double) fastrand() / (FASTRAND_MAX)));

		if (p > dice) {	//probability

			accepted++;

			if (AorB == 'A') {

				molA->swapA(r, c);
				compare->updateDifferenceA2(r, c);

			}

			if (AorB == 'B') {

				molB->swapB(r, c);
				compare->updateDifferenceB2(r, c);
			}

			if (compare->getScore() > best) {

				best = compare->getScore();

			}

		} else {

			rejected++;

		}

		accRate = (double) accepted / (double) rejected;

		if (accRate < 0.003) {

			if (verbose) {
				molA->printSeq(molA->getLength());
				molB->printSeq(molB->getLength());
			}

			break;

		}
		loop++;
	}

	return best;
}

void run_mol(int argc, char* argv[]) {

	blosum = GetCurrentWorkingDir() + "/MOL2";
	p_values_mol2();

	ifstream fin(argv[2]); //to count stuff
	ifstream fin2(argv[2]);
	//according to algorithm with null correspondences, B molecule is always the largest of the two
	Molecule* mol1 = new Molecule(&fin, &fin2, "mol");

	mol1->cacheDist();
	int length1, length2;

	length1 = mol1->getLength();

	if (verbose)
		cout << "reading molecules...\n";

	int molN = 0;
	std::vector<Molecule*> mols2 = read_list_of_mol2s(argv[3]);
	std::vector<Molecule*> mols3;
	std::ofstream out;

	if(superimpose){

		mols3 = copy_list_of_mol2s(argv[3]);
		out.open(superfile.c_str());

	}

	cout << "TARGET_NAME\tDATABASE_NAME\tSTR_SCORE\tP_VALUE\tSEQ_SCORE\tTARGET_LENGTH\tDB_LENGTH\n" << flush;

	int moleculeIndex = 0;

	for(std::vector<Molecule*>::iterator mol2 = mols2.begin(); mol2 != mols2.end(); ++mol2){

		srand(seed);
		fast_srand(seed);
		(*mol2)->cacheDist();

		length2 = (*mol2)->getLength();

		nulls = round(double(min(length1, length2)) / 100.0 * percent);

		if (verbose)
			cout << "Nulls: " << nulls << "\n";

		if (length1 > length2) {

			molA = (*mol2);
			molB = mol1;

		} else {

			molA = mol1;
			molB = (*mol2);

		}

		if (verbose) {

			cout << "Molecule A length: " << molA->getLength() << "\n";
			cout << "Molecule B length: " << molB->getLength() << "\n";

		}

		rs.Clear();

		if(molN == 0){
			if (verbose)
				cout << "Creating InterLig object..." << "\n" << flush;

			compare = new InterLig(molA, molB);
			compare->loadMOL2(blosum.c_str());
			compare->initializeLgmap();
			compare->setD0(d0);

		}else{

			if(verbose)
				cout << "Resetting InterLig object...\n" << flush;
			compare->reset(molA, molB);
			compare->setD0(d0);
		}
		molN++;

		if (verbose)
			cout << "Done." << "\n";

		for (int r = 0; r < rounds; r++) {

			anneal();

		}

		int pvalue_xindex = int(molA->getLength() - 5);

		pvalue_xindex = min(pvalue_xindex, 35);
		pvalue_xindex = max(pvalue_xindex, 0);

		int pvalue_zindex = round(best * 1000);

		double best_pvalue = pvaluesMol2[pvalue_xindex][pvalue_zindex];

		printf("%11s\t%13s\t%9.2f\t%7.2f\t%9.0f\t%13d\t%9d\t\n", (mol1->getName()).c_str(), ((*mol2)->getName()).c_str(), best, best_pvalue, compare->getsScore(), mol1->getLength(), (*mol2)->getLength());

		best = 0;

		if (superimpose) {

			//fp1 = fopen(filename_superA.c_str(), "w+");

			//rotates mol2 onto mol1
			compare->superimpose(mol1, (*mol2));

			double* firstTrans = compare->getTranslation2();

			for (int i = 0; i < 3; i++) {
				firstTrans[i] = -firstTrans[i];
			}

			mols3[moleculeIndex]->rototranslate(compare->getI(), firstTrans);
			mols3[moleculeIndex]->rototranslate(compare->getRotation(),
												compare->getTranslation1());
			mols3[moleculeIndex]->printMolecule(&out);

		}

		(*mol2)->destroy();
		moleculeIndex++;
	}
	if(superimpose)
		out.close();

}

int main(int argc, char *argv[]) {

	fast_srand(time(NULL));
	srand(time(NULL));

	if(argc == 1){
		cout << "Usage: ./InterLig -mol2 target.mol2 database.mol2 <options>\n"
				<< "\nFor more help: InterLig -h\n";
		return EXIT_FAILURE;
	}
	if (!strcmp(argv[1], "-h")) {

		cout << "\nUsage: ./InterLig  -mol2 target.mol2 database.mol2 <options>\n\n" <<
				"   -mol2  : Compare a mol2 file with a mol2 database (multiple molecules in the database)\n" <<
				"          (./InterLig -mol2 examples/target.mol2 examples/database.mol2)\n" <<

				"\nGeneral options:\n\n" <<

				"   -v             : verbose output\n" <<
				"   -h             : print this help\n" <<
				"   -matrix <path> : optional path to the sequence similarity matrix (default: <current_path>/MOL2)\n"  <<

				"\nAlignment options:\n\n" <<

				"   -d0            : d0 parameter in the Levitt-Gerstein score (default: 0.5)\n" <<
				"   -dW            : weight of sequence similarity (e.g. BLOSUM62) in the optimization procedure (default: 0.5)\n" <<
				"   -eps           : weight of structural similarity in the optimization procedure (default: 0.5)\n" <<
				"   -nullP         : percentage of ignored atoms in the smaller molecule (default: 0)\n" <<
				"   -super <path>  : calculate optimal superposition between aligned molecules\n" <<

				"\nAnnealing options:\n\n" <<
				"   -anneal <int>  : number of annealing rounds (default: 1)\n" <<
				"   -ch <int>      : Markov chain length multiplier (default: 10)\n" <<
				"   -T0 <float>    : Initial annealing temperature (default: 0.2)\n" <<
				"   -seed <int>    : seed for the stochastic process (random number generation)\n";

	}else if (!strcmp(argv[1], "-mol2")) {

		parse_options(argc, argv);
		run_mol(argc, argv);

	} else {

		cout << "Usage: ./InterLig -mol2 target.mol2 database.mol2 <options>\n"
				<< "For more help: InterLig -h\n";
		return EXIT_FAILURE;
	}

}

