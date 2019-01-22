/*
 * Molecule.h
 *
 *  Created on: Feb 10, 2016
 *      Author: claudio
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string>

#include <malloc.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>        // std::abs

using namespace std;

class Molecule{

private:

	char Emap(string AA);
	void get_xyz(string line, double *x, double *y, double *z, char *atomname);
	void RemoveSpaces(char* source);
	int readMol(char *filename);
	int readMol(ifstream *fin);
	int readFullMol(ifstream* fin);

	//private sets
	void setLength(char *filename);
	void setLength(ifstream* fin);
	void setFullLength(ifstream* fin);

	void setDist();
	void setDistSqrt();
	void setDiff();

	int length;

	double **coord;		//coordinates[length][3]
	int    *resnumber;
	char  **atomname;
	char   *seq;		//atom sequences
	char   *chain;		//protein sequences
	std::string fn;
	std::string name;
	const char   *type;

	bool *null;

	double **dist;  //distance matrix

	std::vector<std::string> textLines;

public:

	Molecule(char* filename);
	Molecule(ifstream* fin, ifstream* fin2, const char* t);

	void cacheDist();
	void destroy();

	int getLength();

	char* getSeq();
	double** getDist();
	double* getCoord(int i);
	int getResN(int i);
	void setCoord(int i, int j, double newCoord);
	void setName(std::string n);
	std::string getName();
	bool* getNull();

	void printSeq(int size);

	void printSeq2();
	void printDist();
	void printMolecule(std::ofstream* os);

	void swapA(int i, int j);
	void swapB(int i, int j);

	void resetNull();
	void setNull(int i);
	void shuffle();
	void shuffle(int rounds);
	void rototranslate(double** rot, double* trans);
	std::string getFileName();
	double* center_molecule(int length, double* distances);


};

Molecule::Molecule(char* filename){

	fn = filename;
	setLength(filename);

	coord = new double*[this->length];
	resnumber = new int[this->length];
	atomname = new char*[this->length];

	for(int i=0; i<this->length; i++){
		coord[i] = new double[3];
		resnumber[i] = 0;
		atomname[i] = new char[5];
	}

	this->seq   = new char[this->length+1];
	this->chain   = new char[this->length+1];
	null  = new bool[this->length];

	seq[this->length] = '\0';
	chain[this->length] = '\0';
	//load coordinates
	this->readMol(filename);

}

Molecule::Molecule(ifstream* fin, ifstream* fin2, const char* t){

	type = t;

	if(!strcmp(type, "mol")){
		setLength(fin);
	}else{
		setFullLength(fin);
	}

	coord = new double*[this->length];
	resnumber = new int[this->length];
	atomname = new char*[this->length];

	for(int i=0; i<this->length; i++){
		coord[i] = new double[3];
		resnumber[i] = 0;
		atomname[i] = new char[5];
	}

	this->seq   = new char[this->length+1];
	this->chain   = new char[this->length+1];
	null  = new bool[this->length];

	seq[this->length] = '\0';
	chain[this->length] = '\0';

	//load
	if(!strcmp(type, "mol")){
		this->readMol(fin2);
	}else{
		this->readFullMol(fin2);
	}

	dist = new double*[this->length];

	for(int x=0; x<this->length; x++){

		dist[x] = new double[this->length];
		null[x] = false;
	}

	this->setDist();

}

void Molecule::cacheDist(){
	dist = new double*[this->length];


	for(int x=0; x<this->length; x++){

		dist[x] = new double[this->length];
		null[x] = false;
	}

//		this->setDistSqrt();
		this->setDist();

}

void Molecule::rototranslate(double** rot, double* trans){

	for(int i=0; i<this->length; i++){

		double x = rot[0][0] * this->getCoord(i)[0] + rot[0][1] * this->getCoord(i)[1]  + rot[0][2] * this->getCoord(i)[2] + trans[0];
		double y = rot[1][0] * this->getCoord(i)[0] + rot[1][1] * this->getCoord(i)[1]	+ rot[1][2] * this->getCoord(i)[2] + trans[1];
		double z = rot[2][0] * this->getCoord(i)[0] + rot[2][1] * this->getCoord(i)[1]	+ rot[2][2] * this->getCoord(i)[2] + trans[2];

		this->setCoord(i, 0, x);
		this->setCoord(i, 1, y);
		this->setCoord(i, 2, z);

	}
}

double* Molecule::center_molecule(int length, double* distances){

	int	i;	/* Counter variables */
	int   natoms;  // Number of selected atoms
	double	xcen,ycen,zcen;		/* Temporary coordinates values */
	double* vec = new double[3];

	natoms=0;
	xcen=ycen=zcen=0;
	for (i=0; i<length; i++){
		if(distances[i] < 5){
			xcen+=this->getCoord(i)[0];
			ycen+=this->getCoord(i)[1];
			zcen+=this->getCoord(i)[2];
			natoms++;
		}
	}
	/* Now center molecule */
	xcen/=(double)natoms;
	ycen/=(double)natoms;
	zcen/=(double)natoms;

	vec[0]=xcen;
	vec[1]=ycen;
	vec[2]=zcen;

	for (i=0;i<this->getLength();i++)
	{

		this->setCoord(i, 0, this->getCoord(i)[0]-xcen);
		this->setCoord(i, 1, this->getCoord(i)[1]-ycen);
		this->setCoord(i, 2, this->getCoord(i)[2]-zcen);

	}

	return vec;
}

char Molecule::Emap(string AA){

	char A=' ';
	if(     AA.compare("Br")==0)   A='A';
	else if(AA.compare("C.1")==0)   A='B';
	else if(AA.compare("C.2")==0)   A='C';
	else if(AA.compare("C.3")==0)   A='D';
	else if(AA.compare("C.ar")==0)   A='E';
	else if(AA.compare("Cl")==0)   A='F';
	else if(AA.compare("F")==0)   A='G';
	else if(AA.compare("F.3")==0)   A='H';
	else if(AA.compare("H")==0)   A='I';
	else if(AA.compare("N.1")==0)   A='L';
	else if(AA.compare("N.2")==0)   A='M';
	else if(AA.compare("N.3")==0)   A='N';
	else if(AA.compare("N.4")==0)   A='O';
	else if(AA.compare("N.am")==0)   A='P';
	else if(AA.compare("N.ar")==0)   A='Q';
	else if(AA.compare("N.pl3")==0)   A='R';
	else if(AA.compare("O.2")==0)   A='S';
	else if(AA.compare("O.3")==0)   A='T';
	else if(AA.compare("O.co2")==0)   A='U';
	else if(AA.compare("S.2")==0)   A='V';
	else if(AA.compare("S.3")==0)   A='W';
	else if(AA.compare("S.02")==0)   A='X';
	else A='Z';

	return A;
}

void Molecule::get_xyz(string line, double *x, double *y, double *z, char *atomname){

	char cstr[50];

	strcpy(cstr, (line.substr(17, 9)).c_str());
	sscanf(cstr, "%lf", x);

	strcpy(cstr, (line.substr(27, 9)).c_str());
	sscanf(cstr, "%lf", y);

	strcpy(cstr, (line.substr(37, 9)).c_str());
	sscanf(cstr, "%lf", z);

	strcpy(cstr, (line.substr(47, 9)).c_str());
	RemoveSpaces(cstr);
	*atomname=Emap(cstr);

}

void Molecule::RemoveSpaces(char* source){
	char* i = source;
	char* j = source;
	while(*j != 0)
	{
		*i = *j++;
		if(*i != ' ')
			i++;
	}
	*i = 0;
}

void Molecule::setLength(char *filename){

	int i=0;
	string line="";

	ifstream fin(filename);

	string atom("@<TRIPOS>ATOM");
	string bond("@<TRIPOS>BOND");

	if (fin.is_open()){

		while (fin.good()){

			getline(fin, line);

			if(line.compare(0, atom.length(), atom)==0){

				getline(fin, line);

				while(line.compare(0, bond.length(), bond)){

					if(line.at(47) != 'H'){
						i++;
					}
					getline(fin, line);

				}
			}
		}
		fin.close();
	}else{
		char error[5000];
		sprintf(error, "Can not open file: %s", filename);
		cerr << error << endl;
		exit(1);

	}

	if(i==0){
		char error[5000];
		sprintf(error, "Can not find atoms in file: %s", filename);
		cerr << error << endl;
		exit(1);

	}

	this->length = i;

}

void Molecule::setLength(ifstream* fin){

	int i=0;
	string line="";
	string n="";
	bool readname = false;

	string molecule("@<TRIPOS>MOLECULE");
	string atom ("@<TRIPOS>ATOM");
	string bond ("@<TRIPOS>BOND");

	if ((*fin).is_open()){

		while (line.compare(0, bond.length(), bond) && (*fin).good()){

			getline(*fin, line);

			if(readname){
				n = line;
				readname = false;
			}

			if(line.compare(0, molecule.length(), molecule)==0){
				readname = true;
			}

			if(line.compare(0, atom.length(), atom)==0){

				getline(*fin, line);

				while(line.compare(0, bond.length(), bond)){

					if(line.at(47) != 'H'){
						i++;
					}

					getline(*fin, line);

				}
			}
		}

	}else{

		char error[5000];
		sprintf(error, "Read error");
		cerr << error << endl;

	}

	if(i==0){

		char error[5000];
		sprintf(error, "Read error, size zero");
		cerr << error << endl;

	}

	this->length = i;
	this->setName(n);

}

void Molecule::setFullLength(ifstream* fin){

	int i=0;
	string line="";
	string n="";
	bool readname = false;

	string molecule("@<TRIPOS>MOLECULE");
	string atom ("@<TRIPOS>ATOM");
	string bond ("@<TRIPOS>BOND");

	if ((*fin).is_open()){

		while (line.compare(0, bond.length(), bond) && (*fin).good()){

			getline(*fin, line);

			if(readname){
				n = line;
				readname = false;
			}

			if(line.compare(0, molecule.length(), molecule)==0){
				readname = true;
			}

			if(line.compare(0, atom.length(), atom)==0){

				getline(*fin, line);

				while(line.compare(0, bond.length(), bond)){

					//if(line.at(47) != 'H'){
					i++;
					//}

					getline(*fin, line);

				}
			}
		}

	}else{

		char error[5000];
		sprintf(error, "Read error");
		cerr << error << endl;

	}

	if(i==0){

		char error[5000];
		sprintf(error, "Read error, size zero");
		cerr << error << endl;

	}

	this->length = i;
	this->setName(n);

}

int Molecule::readMol(char *filename){
	int i=0;
	string line, str;

	ifstream fin (filename);

	string atom ("@<TRIPOS>ATOM");
	string bond ("@<TRIPOS>BOND");

	if (fin.is_open()){

		while (fin.good()){

			getline(fin, line);
			if(line.compare(0, atom.length(), atom)==0){

				getline(fin, line);

				while(line.compare(0, bond.length(), bond)){

					if(line.at(47) != 'H'){

						get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i]);
						i++;
					}
					getline(fin, line);

				}
			}
		}
		fin.close();

	}

	return 0;
}

int Molecule::readMol(ifstream* fin){

	int i=0;
	string line, str;

	getline(*fin, line);

	string molecule("@<TRIPOS>MOLECULE");
	string atom("@<TRIPOS>ATOM");
	string bond("@<TRIPOS>BOND");

	if ((*fin).is_open()){

		getline(*fin, line);


		while (line.compare(0, bond.length(), bond)!=0 && (*fin).good()){//(*fin).good()){

			getline(*fin, line);

			if(line.compare(0, atom.length(), atom)==0){

				getline(*fin, line);


				while(line.compare(0, bond.length(), bond)){

					if(line.at(47) != 'H'){
						//cout << "Loading coord " << i << " " << line.at(47) << "\n";
						get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i]);
						i++;
					}
					getline(*fin, line);

				}
			}
		}

	}

	return 0;
}

int Molecule::readFullMol(ifstream* fin){

	int i=0;
	string line, str;

	getline(*fin, line);
	textLines.push_back(line);
	string molecule("@<TRIPOS>MOLECULE");
	string atom("@<TRIPOS>ATOM");
	string bond("@<TRIPOS>BOND");
	streampos oldpos;

	if ((*fin).is_open()){

		getline(*fin, line);
		textLines.push_back(line);

		while (line.compare(0, bond.length(), bond)!=0 && (*fin).good()){

			getline(*fin, line);

			textLines.push_back(line);

			if(line.compare(0, atom.length(), atom)==0){

				getline(*fin, line);
				textLines.push_back(line);

				while(line.compare(0, bond.length(), bond)){

					//if(line.at(47) != 'H'){
					get_xyz(line, &coord[i][0], &coord[i][1], &coord[i][2], &seq[i]);
					i++;
					//}
					getline(*fin, line);
					textLines.push_back(line);
				}
			}
		}

		getline(*fin, line);

		while (line.compare(0, molecule.length(), molecule)!=0 && (*fin).good()){

			textLines.push_back(line);
			oldpos = (*fin).tellg();
			getline(*fin, line);

		}

		(*fin).seekg(oldpos);
	}

	return 0;
}

void Molecule::printMolecule(std::ofstream* os){

	string atom("@<TRIPOS>ATOM");
	string bond("@<TRIPOS>BOND");

	bool pastAtom = false;
	bool pastBond = false;
	int thisatom = 0;
	for(std::vector<string>::iterator text = textLines.begin(); text != textLines.end(); ++text){

		if((*text).compare(0, bond.length(), bond)==0){
			pastBond = true;
		}

		if(pastAtom && !pastBond){

			std::stringstream ss;

			ss << std::fixed << std::setw(9) << std::setprecision(4) << coord[thisatom][0];
			string newx = ss.str();
			ss.str(std::string());

			ss << std::fixed << std::setw(9) << std::setprecision(4) << coord[thisatom][1];
			string newy = ss.str();
			ss.str(std::string());

			ss << std::fixed << std::setw(9) << std::setprecision(4) << coord[thisatom][2];
			string newz = ss.str();
			ss.str(std::string());

			std::string newtext = (*text).replace(17, 9, newx);
						newtext = (*text).replace(27, 9, newy);
						newtext = (*text).replace(37, 9, newz);
			(*os) << newtext << "\n";

			thisatom++;

		}else{
			(*os) << *text << "\n";
		}

		if((*text).compare(0, atom.length(), atom)==0){
			pastAtom = true;
		}

	}
}

//fills the distance matrices the first time
void Molecule::setDist(){

	for(int i=0; i<this->length; i++){

		for(int j=i; j<this->length; j++){

			if(coord[i][0] != 0 && coord[j][0] != 0 && i!=j){
				this->dist[i][j] = this->dist[j][i] = (coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0])+(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1])+(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]);
			}else{
				this->dist[i][j] = this->dist[j][i] = 0;
			}

		}

	}

}

//fills the distance matrices the first time
void Molecule::setDistSqrt(){

	for(int i=0; i<this->length; i++){

		for(int j=i; j<this->length; j++){

			if(coord[i][0] != 0 && coord[j][0] != 0 && i!=j){
				this->dist[i][j] = this->dist[j][i] = sqrt((coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0])+(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1])+(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]));
			}else{
				this->dist[i][j] = this->dist[j][i] = 0;
			}

		}

	}

}

int Molecule::getLength(){

	return this->length;

}

std::string Molecule::getFileName(){
	return this->fn;
}

char* Molecule::getSeq(){

	return this->seq;

}

double** Molecule::getDist(){

	return this->dist;

}

double* Molecule::getCoord(int i){

	return this->coord[i];

}

void Molecule::setCoord(int i, int j, double newCoord){

	this->coord[i][j] = newCoord;

}

int Molecule::getResN(int i){

	return this->resnumber[i];

}

bool* Molecule::getNull(){

	return this->null;

}

//swapping A is much simpler, it is only about setting a new null point
void Molecule::swapA(int c, int r){

	if(!null[r] and null[c]){
		null[r] = true;
		null[c]= false;
	}else if(null[r] and !null[c]){
		null[r] = false;
		null[c] = true;
	}else{
		cout << "Something's wrong in swapA\n";
	}
}

//swap two columns/rows on the distance matrix
//(use only when you have accepted the move)
void Molecule::swapB(int i, int j){

	int m;

	//swap coordinates
	std::swap(this->coord[i], this->coord[j]);

	//normal D swapping, no need to recalculate distances
	//swap rows in the distance matrix
	std::swap(this->dist[i], this->dist[j]);

	std::swap(this->dist[i][j], this->dist[j][j]);
	std::swap(this->dist[j][i], this->dist[i][i]);

	for(m=0; m<this->length; m++){

		if(m!=i && m!=j){

			this->dist[m][i] = this->dist[i][m];
			this->dist[m][j] = this->dist[j][m];
		}
	}

	//swap letters in the sequence
	std::swap(this->seq[i], this->seq[j]);
	std::swap(this->resnumber[i], this->resnumber[j]);

}

void Molecule::shuffle(){

	int i = 0;
	int j = 0;

	for(int r = 0; r<10000; r++){

		do{
			i = rand() % this->length;
			j = rand() % this->length;
		}while(i == j);

		this->swapB(i, j);

	}

}

void Molecule::shuffle(int rounds){

	int i = 0;
	int j = 0;

	for(int r = 0; r<rounds; r++){

		do{
			i = rand() % this->length;
			j = rand() % this->length;
		}while(i == j);

		this->swapB(i, j);

	}

}

void Molecule::resetNull(){


	for(int i = 0; i< this->length; i++){

		this->null[i] = false;

	}

}

void Molecule::setNull(int i){

	this->null[i] = true;

}

void Molecule::setName(std::string n){

	this->name = n;

}

std::string Molecule::getName(){

	return this->name;

}

void Molecule::printSeq(int size){

	cout << "> Molecule";

	cout << "\n";

	for(int a=0; a<size; a++){
		if(!this->null[a]){
			cout << seq[a] << "";
		}else{
			cout << "-" << "";
		}
	}

	cout << "\n";
}

void Molecule::printSeq2(){

	for(int a=0; a<this->length; a++){
		if(this->null[a]){
			cout << " ";
		}else{
			cout << "-";
		}
	}
	cout << "\n";

	for(int a=0; a<this->length; a++){

		cout << seq[a] << "";
	}

	cout << "\n";
}

void Molecule::printDist(){

	for(int a=0; a<this->length; a++)
		cout << seq[a] << "\t";

	cout << "\n";

	for(int i=0; i<this->length; i++){
		for(int j=0; j<this->length; j++){

			if(!this->null[i] && !this->null[j]){
				cout << dist[i][j] << "\t";
			}else{
				cout << "0\t";
			}
		}
		cout << "\n";
	}
}

void Molecule::destroy(){

	//fix this
	delete[] seq;
	delete[] null;
	delete[] chain;
	delete[] resnumber;
	delete[] atomname;


	for(int i=0; i<this->length; i++){
		//cout << "deleting row: " << i << " /" << length << "\n" << flush;
		delete[] dist[i];
		delete[] coord[i];
	}

	delete[] dist;
	delete[] coord;

}

#endif /*  */
