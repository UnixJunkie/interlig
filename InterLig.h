/*
 * InterLig.h
 *
 *  Created on: Feb 10, 2016
 *      Author: claudio
 */

#ifndef INTERLIG_H_
#define INTERLIG_H_

#define max(a,b) a>b?a:b
#define min(a,b) a>b?b:a

#define neighborhood 1000000
#define DBL_MAX 100000000.00
#define PI 3.14159265358979323846  /* pi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include "RunningStat.h"
#include <sstream>
#include <iostream>
#include <fstream>

//#include <algorithm>
#include <cmath>        // std::fabs

#include "Molecule.h"

class InterLig{

	private:

		double **diff;  //distance matrix
		double *diffCumulative;
		double score;
		double dScore;
		double sScore;

		double de;


		double d0;

		int length;

		int contributes; //the number of contributes to the score

		Molecule* molA;
		Molecule* molB;

		double deltaD;
		double deltaN;
		double deltaS;

		int** MOL2;
		int* AAtoi;

		double** I;
        double** s;
        double* trans1;
        double* trans2;

        int lgmapsize; //10 angstrom at 0.0001 steps
        double* lgmap;


	public:

		InterLig(Molecule* molA, Molecule* molB);
		void reset(Molecule* molA, Molecule* molB);
		void destroy();

		void cacheDifference();
		void differencePerResidue();
		void updateDifference(int i, int j);

		void updateDifferenceA(int i, int j);
		void updateDifferenceA2(int i, int j);

		void updateDifferenceB(int i, int j);
		void updateDifferenceB2(int i, int j);

		void deltaA(int i, int j);
		void deltaB(int i, int j);

		void initializeLgmap();
		double levittGerstein(double d);

		void deltaSeqA(int i, int j);
		void deltaSeqB(int i, int j);

		double temperatureSD(double mi, double sigma, double Told);

		double getScore();
		double getsScore();
		double getdScore();
		double getDeltaD();
		double getDeltaN();
		double getDeltaS();
                double** getRotation();
                double** getI();
                double* getTranslation1();
                double* getTranslation2();
		int getContributes();
		double* getDiffCumulative();
                //double* center_molecule(Molecule *m);
		double superimpose(Molecule* A, Molecule* B);
                
		void transpose_matrix(double** m);
		void multiply_matrix(double** a, double** b, double**);
		void copy_matrix(double** f, double** t);

		void loadMOL2(const char* file);

		void setD0(double d);
		void printDiff();

};

InterLig::InterLig(Molecule* A, Molecule* B){

	this->molA = A;
	this->molB = B;

	this->d0 = 1;

	this->length = A->getLength();
	this->score = 0;
	this->dScore = DBL_MAX;
	this->contributes = 0;

    this->diff = new double*[this->length];
    this->diffCumulative = new double[this->length];

    this->lgmapsize = 10000000;
    this->lgmap = new double[lgmapsize];

	for(int x=0; x<this->length; x++){
		diffCumulative[x] = 0;
		diff[x]  = new double[this->length];

		for(int y=0; y<this->length; y++){

			diff[x][y] = 0;

		}
	}


}

void InterLig::reset(Molecule* A, Molecule* B){

	this->destroy();
	this->molA = A;
	this->molB = B;

	this->d0 = 1;

	this->length = A->getLength();

	this->score = 0;
	this->dScore = DBL_MAX;
	this->contributes = 0;

    this->diff = new double*[this->length];
    this->diffCumulative = new double[this->length];

	for(int x=0; x<this->length; x++){
		diffCumulative[x] = 0;
		diff[x]  = new double[this->length];
		for(int y=0; y<this->length; y++){

			diff[x][y] = 0;

		}
	}

}


double InterLig::superimpose(Molecule* m1, Molecule* m2) {

	int i, j, k; /* Counter variables */
	int natoms = min(m1->getLength(), m2->getLength());
	double** u = new double*[3]; /* direct product matrix */
	I = new double*[3]; /* identity matrix */
	double** t = new double*[3]; /* Temporary storage matrix */
	double** ma = new double*[3]; /* x axis rotation matrix */
	double** mb = new double*[3]; /* y axis rotation matrix */
	double** mg = new double*[3]; /* z axis rotation matrix */
	double *d1, *d2; /* useful pointers */
	double error; /* Final error */
	double error2; /* Final error */
	double alpha = 0.0; /* Angle of rotation around x axis */
	double beta = 0.0; /* Angle of rotation around y axis */
	double gamma = 0.0; /* Angle of rotation around z axis */
	double x, y, z; /* Temporary coordinate variables */
	s = new double*[3]; /* Final transformation matrix */
	double error_cut = 0.0000001;

	this->differencePerResidue();

	for (i = 0; i < 3; i++){ /* Initialize matrices */
		s[i] = new double[3];
		u[i] = new double[3];
		t[i] = new double[3];
		for (j = 0; j < 3; j++){
			s[i][j] = u[i][j] = t[i][j] = 0.0;
		}
	}

	s[0][0] = s[1][1] = s[2][2] = 1.0; /* Initialize S matrix to I */
	for (i = 0; i < 3; i++){ /* Initialize rotation matrices to I */
		ma[i] = new double[3];
		mb[i] = new double[3];
		mg[i] = new double[3];
                I[i] = new double[3];
		for (j = 0; j < 3; j++){
			ma[i][j] = mb[i][j] = mg[i][j] = I[i][j] = s[i][j];
		}
	}
	
	trans1 = m1->center_molecule(natoms, this->diffCumulative);
    trans2 = m2->center_molecule(natoms, this->diffCumulative);
        
	for (i = 0; i < natoms; i++){ /* Construct U matrix */

		if (this->diffCumulative[i] < 5.0) {
			d1 = m1->getCoord(i); //m1->atm[i].x

			for (j = 0; j < 3; j++){
                            d2 = m2->getCoord(i); //m2->atm[i].x

				for (k = 0; k < 3; k++){

					u[j][k] += (*d1) * (*d2);
                                    d2++;

				}
				d1++;
			}
		}
	}

	do{
		error = 0.0;
		/* Calculate x axis rotation */
		alpha = atan((u[2][1] - u[1][2]) / (u[1][1] + u[2][2]));

		/* Insure we are heading for a minimum, not a maximum */
		if (cos(alpha) * (u[1][1] + u[2][2]) + sin(alpha) * (u[2][1] - u[1][2])
				< 0.0)
			alpha += PI;
		ma[1][1] = ma[2][2] = cos(alpha);
		ma[2][1] = sin(alpha);
		ma[1][2] = -ma[2][1];
		transpose_matrix(ma);
		multiply_matrix(u, ma, t);
		transpose_matrix(ma);
		copy_matrix(t, u);
		multiply_matrix(ma, s, t);
		copy_matrix(t, s);
		/* Calculate y axis rotation */
		beta = atan((u[0][2] - u[2][0]) / (u[0][0] + u[2][2]));
		/* Insure we are heading for a minimum, not a maximum */
		if (cos(beta) * (u[0][0] + u[2][2]) + sin(beta) * (u[0][2] - u[2][0])
				< 0.0)
			beta += PI;
		mb[0][0] = mb[2][2] = cos(beta);
		mb[0][2] = sin(beta);
		mb[2][0] = -mb[0][2];
		transpose_matrix(mb);
		multiply_matrix(u, mb, t);
		transpose_matrix(mb);
		copy_matrix(t, u);
		multiply_matrix(mb, s, t);
		copy_matrix(t, s);
		/* Calculate z axis rotation */
		gamma = atan((u[1][0] - u[0][1]) / (u[0][0] + u[1][1]));
		/* Insure we are heading for a minimum, not a maximum */
		if (cos(gamma) * (u[0][0] + u[1][1]) + sin(gamma) * (u[1][0] - u[0][1])
				< 0.0)
			gamma += PI;
		mg[0][0] = mg[1][1] = cos(gamma);
		mg[1][0] = sin(gamma);
		mg[0][1] = -mg[1][0];
		transpose_matrix(mg);
		multiply_matrix(u, mg, t);
		transpose_matrix(mg);
		copy_matrix(t, u);
		multiply_matrix(mg, s, t);
		copy_matrix(t, s);
		error = fabs(alpha) + fabs(beta) + fabs(gamma);
	} while (error > error_cut); /* was 0.0001 before optimization Is error low enough to stop? */
	/* Now calculate final RMS superimposition */
	error = 0.0;
	error2 = 0.0;
	//natoms = m1->getLength();
	//printf("testing7b %f %f\n",m2->atm[1].x,error);
	//printf ("s1:\t%f\t%f\t%f\n",s[0][0],s[0][1],s[0][2]);
	//printf ("s2:\t%f\t%f\t%f\n",s[1][0],s[1][1],s[1][2]);
	//printf ("s3:\t%f\t%f\t%f\n",s[2][0],s[2][1],s[2][2]);
	if (!(s[0][0] > 0 || s[0][0] < 0)) {
		//      printf ("bar\n");
		for (i = 0; i < 3; i++) /* Initialize matrices */
			for (j = 0; j < 3; j++)
				s[i][j] = u[i][j] = 0.0;
		s[0][0] = s[1][1] = s[2][2] = 1.0; /* Initialize S matrix to I */
	}
	//for (i = 0; i < m1->getLength(); i++) {

		//dist=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
		//printf("TEST1-m2\t %d\t%6.3f %6.3f %6.3f\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist);
		//x = s[0][0] * m2->getCoord(i)[0] + s[0][1] * m2->getCoord(i)[1] + s[0][2] * m2->getCoord(i)[2] + trans1[0];
		//y = s[1][0] * m2->getCoord(i)[0] + s[1][1] * m2->getCoord(i)[1]	+ s[1][2] * m2->getCoord(i)[2] + trans1[1];
		//z = s[2][0] * m2->getCoord(i)[0] + s[2][1] * m2->getCoord(i)[1]	+ s[2][2] * m2->getCoord(i)[2] + trans1[2];

		//m2->setCoord(i, 0, x);
		//m2->setCoord(i, 1, y);
		//m2->setCoord(i, 2, z);
        m2->rototranslate(s, trans1);
        m1->rototranslate(I, trans1);
                //bring back m1 to its previous position
                //m1->setCoord(i, 0, m1->getCoord(i)[0]+trans1[0]);
                //m1->setCoord(i, 1, m1->getCoord(i)[1]+trans1[1]);
                //m1->setCoord(i, 2, m1->getCoord(i)[2]+trans1[2]);
        

		//dist2=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
		//printf("TEST2-m2\t %d\t%6.3f %6.3f %6.3f\t%e\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist2,dist2-dist);
        for (i = 0; i < natoms; i++) {
        	if (this->diffCumulative[i] < 5.0) {
        		x = m1->getCoord(i)[0] - m2->getCoord(i)[0];
        		y = m1->getCoord(i)[1] - m2->getCoord(i)[1];
        		z = m1->getCoord(i)[2] - m2->getCoord(i)[2];

				error += x * x + y * y + z * z;
        	}
		//printf("TEST:\t %d\t%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n",i,m1->getCoord(i)[0],m2->getCoord(i)[0],m1->getCoord(i)[1],m2->getCoord(i)[1],m1->getCoord(i)[2],m2->getCoord(i)[2]);
		//printf("TEST1:\t %d\t%f %f %f\n",i,x*x,y*y,z*z);
		//printf("TEST2:\t %d %d\t%f\n",i,m1->atm[i].selected,m1->atm[i].rms);
		//if (m1->atm[i].selected) {
		//	error2 += m1->atm[i].rms;
		//	natoms++;
		//}
		/*printf ("TEST1: %f %f %f\n",m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
		 printf ("TEST2: %f %f %f\n",m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
		 printf ("TEST3: %d %f %f\n",i,error,m1->atm[i].rms); */
	}
	//    printf("testing8 %f\n",m2->atm[1].x);
	error /= (double) (natoms);
	//error2 /= (double) (natoms);
	if (natoms == 0)
		error2 = 99999999;
	return (sqrt(error));
	return 0;
}

void InterLig::transpose_matrix(double** m) {/* Transpose a 3x3 matrix */
	double dummy;
	dummy = m[0][1];
	m[0][1] = m[1][0];
	m[1][0] = dummy;
	dummy = m[0][2];
	m[0][2] = m[2][0];
	m[2][0] = dummy;
	dummy = m[1][2];
	m[1][2] = m[2][1];
	m[2][1] = dummy;
}

void InterLig::multiply_matrix(double** a, double** b,
		double** c) { /* computes C=AB */

	int i, j, k;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			c[i][j] = 0.0;
			for (k = 0; k < 3; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
}

void InterLig::copy_matrix(double** f, double** t){ /* copy matrix f into matrix t */

	int i, j;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			t[i][j] = f[i][j];
}

//calculates the full difference between the two matrices
void InterLig::cacheDifference(){

	//cout << dScore << " ";

	this->score = 0;
	this->dScore = 0;
	this->sScore = 0;

	this->contributes = 0;

	double tempSeqScore = 0;

	for(int i=0; i<this->length; i++){
		//cout << i << "\n"<<flush;
		this->diff[i][i] = 0;

		//cout << molA->getSeq()[i] << " " << molB->getSeq()[i] << "\n" << flush;
		if(!molA->getNull()[i])
			tempSeqScore += MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']];

		for(int j=i; j<this->length; j++){

				if(!molA->getNull()[i] && !molA->getNull()[j]){

					this->diff[i][j] = this->diff[j][i] = fabs(molA->getDist()[i][j]-molB->getDist()[i][j]);

					this->dScore += this->diff[i][j];
					this->score +=  1.0/(1.0+(this->diff[i][j]/d0)*(this->diff[i][j]/d0));

					this->contributes++;
				}else{
					this->diff[i][j] = this->diff[j][i] = 0;
				}

		}

	}

	this->sScore = tempSeqScore;

}


void InterLig::differencePerResidue(){


	for(int i=0; i<this->length; i++){

		this->diffCumulative[i] = 0;

		for(int j=i; j<this->length; j++){


				if(!molA->getNull()[i] && !molA->getNull()[j]){

					this->diffCumulative[i] += fabs(molA->getDist()[i][j]-molB->getDist()[i][j]);


				}else{
					this->diffCumulative[i] += 0;
				}

		}
		this->diffCumulative[i] /= (double) this->length;

	}

}


void InterLig::updateDifferenceA(int i, int j){


	//if i is going to be enabled (j disabled)
	if(molA->getNull()[i]){


		//going down the i-j columns for a and b
		for(int r=0; r<this->length; r++){


			if(r!=i && r!=j){//if the difference is lower, delta has to be negative

				           //  - 1*(MOL2[AAtoi[seqA[i]-'A']][AAtoi[seqB[i]-'A']] + MOL2[AAtoi[seqA[j]-'A']][AAtoi[seqB[j]-'A']]);

					if(!molA->getNull()[r]){
						this->dScore += 1*(fabs(molA->getDist()[j][r] - molB->getDist()[j][r]));
						this->score  += 1.0/(1.0+(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0)*(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0));
					}

					if(!molA->getNull()[r]){
						this->dScore -= 1*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r]));
						this->score  -= 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0));
					}


			}


		}

	//if j is going to be enabled (i disabled)
	}else if(molA->getNull()[j]){
	  //de += MOL2[AAtoi[seqA[i]-'A']][AAtoi[seqB[i]-'A']]  - MOL2[AAtoi[seqA[j]-'A']][AAtoi[seqB[j]-'A']];

		for(int r=0; r<this->length; r++){



			if(r!=i && r!=j){//if the difference is lower, delta has to be negative

					if(!molA->getNull()[r]){
						this->dScore += 1*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r]));
						this->score  += 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0));
					}
					if(!molA->getNull()[r]){
						this->dScore -= 1*(fabs(molA->getDist()[j][r] - molB->getDist()[j][r]));
						this->score  -= 1.0/(1.0+(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0)*(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0));
					}

			}
		}



	}

	for(int m=0; m<this->length; m++){

	//if i is going to be enabled (j disabled)
		if(molA->getNull()[i]){
			if(m!=i && m!=j){

				if(!molA->getNull()[m]){

					diff[i][m] = diff[m][i] = fabs(molA->getDist()[i][m]-molB->getDist()[i][m]);

				}



				if(!molA->getNull()[m]){

					diff[j][m] = diff[m][j] = 0;

				}

			}

		//if j is going to be enabled (i disabled)
		}else  if(molA->getNull()[j]){
				if(!molA->getNull()[m]){

					diff[i][m] = diff[m][i] = 0;

				}



				if(!molA->getNull()[m]){

					diff[j][m] = diff[m][j] = fabs(molA->getDist()[j][m]-molB->getDist()[j][m]);

				}
		}
	}

}

void InterLig::updateDifferenceA2(int i, int j){


	this->dScore += deltaD;
	this->score += deltaN;
	this->sScore -= deltaS;

	for(int m=0; m<this->length; m++){

	//if i is going to be enabled (j disabled)
		if(molA->getNull()[i]){
			if(m!=i && m!=j){

				if(!molA->getNull()[m]){

					diff[i][m] = diff[m][i] = fabs(molA->getDist()[i][m]-molB->getDist()[i][m]);

				}



				if(!molA->getNull()[m]){

					diff[j][m] = diff[m][j] = 0;

				}

			}

		//if j is going to be enabled (i disabled)
		}else  if(molA->getNull()[j]){
				if(!molA->getNull()[m]){

					diff[i][m] = diff[m][i] = 0;

				}



				if(!molA->getNull()[m]){

					diff[j][m] = diff[m][j] = fabs(molA->getDist()[j][m]-molB->getDist()[j][m]);

				}
		}
	}

}

void InterLig::updateDifferenceB(int i, int j){


	if(j < this->length){

		//going down the i-j columns for a and b
		for(int r=0; r<this->length; r++){


			if(r!=i && r!=j && !molA->getNull()[i] && !molA->getNull()[r]){//if the difference is lower, delta has to be negative   //



				             	 //old difference   									//new difference
				this->dScore -= 1*(fabs(molA->getDist()[i][r] - molB->getDist()[j][r]) - fabs(molA->getDist()[i][r] - molB->getDist()[i][r]) );///(size);//*size);

				this->score -= 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[j][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[j][r])/d0));
				this->score += 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0));



			}

			if(r!=i && r!=j && !molA->getNull()[j] && !molA->getNull()[r]){//if the difference is lower, delta has to be negative //


				             	 //old difference   	 								//new difference
				this->dScore -= 1*(fabs(molA->getDist()[j][r] - molB->getDist()[i][r]) - fabs(molA->getDist()[j][r] - molB->getDist()[j][r]));///(size);//*size);

				this->score -= 1.0/(1.0+(fabs(molA->getDist()[j][r] - molB->getDist()[i][r])/d0)*(fabs(molA->getDist()[j][r] - molB->getDist()[i][r])/d0));
				this->score += 1.0/(1.0+(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0)*(fabs(molA->getDist()[j][r] - molB->getDist()[j][r])/d0));


			}
		}



	//j is in the outer subset, molA->getNull()[i] can't be true
	}else{

		for(int r=0; r<this->length; r++){



			if(r!=i && r!=j && !molA->getNull()[r]){//

				             //old difference   	 									//new difference
				this->dScore -= 1*(fabs(molA->getDist()[i][r] - molB->getDist()[j][r]) - fabs(molA->getDist()[i][r] - molB->getDist()[i][r]) );///(size);//*size);
				this->score -= 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[j][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[j][r])/d0));
				this->score += 1.0/(1.0+(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0)*(fabs(molA->getDist()[i][r] - molB->getDist()[i][r])/d0));
			}
		}


	}

	for(int m=0; m<this->length; m++){


		if(m!=i && m!=j){



			if(i < this->length){

				if(!molA->getNull()[i] && !molA->getNull()[m]){

					diff[i][m] = diff[m][i] = fabs(molA->getDist()[i][m]-molB->getDist()[i][m]);

				}else{
					diff[i][m] = diff[m][i] = 0;
				}
			}


			if(j < this->length){
				if(!molA->getNull()[j] & !molA->getNull()[m]){

					diff[j][m] = diff[m][j] = fabs(molA->getDist()[j][m]-molB->getDist()[j][m]);

				}else{
					diff[j][m] = diff[j][i] = 0;
				}
			}
		}
	}

	//re-normalise scores with new count
	//*score /= *count;//(sizeMax*sizeMax/4);
	//*dScore /= *count;//(sizeMax*sizeMax/4);
	//cout << dScore << "\n";

}

void InterLig::updateDifferenceB2(int i, int j){


	this->dScore += deltaD;
	this->score += deltaN;
	this->sScore -= deltaS;

	/*illegal moves are:				(for i, j in [0; lengthB), j > i)
		1) swap two residues aligned to null correspondences in A (nullA[i] and nullA[j])
		2) swap a residue in a null correspondence in A with a residue in the outer subset in B (nullA[i] and j >= lengthA)
		3) swap two residues in the outer subset of B (i >= lengthA, j will be too then)
	*/

	for(int m=0; m<this->length; m++){




			//UPDATE THE I-TH COLUMN
			if(!molA->getNull()[i] && !molA->getNull()[m]){


				diff[i][m] = diff[m][i] = fabs(molA->getDist()[i][m]-molB->getDist()[i][m]);


			}else{

				diff[i][m] = diff[m][i] = 0;

			}


			//UPDATE THE J-TH COLUMN
			if(j < this->length){//i and j both in the subset
				if(!molA->getNull()[j] & !molA->getNull()[m]){

					diff[j][m] = diff[m][j] = fabs(molA->getDist()[j][m]-molB->getDist()[j][m]);

				}else{
					diff[j][m] = diff[j][i] = 0;
				}
			}

	}

}



void InterLig::deltaA(int i, int j){

	de=0;

	deltaD = 0;
	deltaN = 0;
	//if i is going to be enabled (j disabled)
	if(molA->getNull()[i]){


		//going down the i-j columns for a and b
		for(int r=0; r<this->length; r++){


			if(r!=i && r!=j){//if the difference is lower, delta has to be negative


					if(!(molA->getNull()[r])){

						de = (fabs(molA->getDist()[i][r] - molB->getDist()[i][r]));///(size*size);

						deltaD += de;
						deltaN += 1.0/(1.0+(de/d0)*(de/d0));

						de = (fabs(molA->getDist()[j][r] - molB->getDist()[j][r]));///(size);//*size);

						//deltaD -= de;
						deltaN -= 1.0/(1.0+(de/d0)*(de/d0));

					}


			}


		}

	//if j is going to be enabled (i disabled)
	}else if(molA->getNull()[j]){

		for(int r=0; r<this->length; r++){



			if(r!=i && r!=j){//if the difference is lower, delta has to be negative

					if(!(molA->getNull()[r])){

						de = (fabs(molA->getDist()[i][r] - molB->getDist()[i][r]));///(size);//*size);
						//deltaD -= de;
						deltaN -= 1.0/(1.0+(de/d0)*(de/d0));

						de = (fabs(molA->getDist()[j][r] - molB->getDist()[j][r]));///(size);//*size);
						deltaD += de;
						deltaN += 1.0/(1.0+(de/d0)*(de/d0));

					}

			}
		}



	}


}

void InterLig::deltaSeqA(int i, int j){

	double de=0;


	//if i is going to be enabled (j disabled)
	if(molA->getNull()[i]){
	  de += - (MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']]  - MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[j]-'A']]);
		//de = MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[j]-'A']] - MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']];

		//if j is going to be enabled (i disabled)
	}else if(molA->getNull()[j]){

	  de +=    MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']]  - MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[j]-'A']];
		//de =    MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']]  - MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[j]-'A']];
	}

	this->deltaS = de;

}

double InterLig::getDeltaD(){

	return deltaD;///(double)this->contributes;

}

double InterLig::getDeltaN(){

	return deltaN;///((double)this->contributes);

}

double InterLig::getDeltaS(){

	return deltaS;///((double)this->contributes);

}

inline double InterLig::levittGerstein(double d){

	//int dindex = int(d*(1000)); //(100000/100) if the distance diff is 10 angstroms, d*10000 is 100,000 and is capped
	//return lgmap[dindex];

	if(d0 == 0.5){
		return 1.0/(1+(d*d*4));
	}else if(d0 == 1){
		return 1.0/(1+(d*d));
	}else{
		return 1.0/(1+((d/d0)*(d/d0)));
	}
}

void InterLig::initializeLgmap(){
	double d = 0;
	double d0mult = 1.0/(d0*d0);
	for(int i = 0; i<lgmapsize; i++){

		d = (double)i/1000.0;
		lgmap[i] = 1.0/(1+(d*d*d0mult));


	}
}

void InterLig::deltaB(int i, int j){

	double de=0;
	//double count=0;
	deltaD = 0;
	deltaN = 0;


	//b[j] and b[i] to be switched. one of them might be a gap (only j, actually. so the b[j] positions shouldn't compared with the a[i] positions, in this case)

	//both in the inner subset
	if(j < this->length){

		//going down the i-j columns for a and b
		for(int r=0; r<this->length; r++){


			if(r!=i && r!=j && !molA->getNull()[i] && !molA->getNull()[r]){//if the difference is lower, delta has to be negative   //


				    de = fabs(molA->getDist()[i][r] - molB->getDist()[j][r]);
				    deltaD += de;
				    deltaN += levittGerstein(de);

					de = fabs(molA->getDist()[i][r] - molB->getDist()[i][r]);
					//deltaD -= de;
					deltaN -= levittGerstein(de);

			}

			if(r!=i && r!=j && !molA->getNull()[j]  && !molA->getNull()[r]){//if the difference is lower, delta has to be negative //


					de = fabs(molA->getDist()[j][r] - molB->getDist()[i][r]);
				    deltaD += de;
				    deltaN += levittGerstein(de);

					de = fabs(molA->getDist()[j][r] - molB->getDist()[j][r]);
				    //deltaD -= de;
				    deltaN -= levittGerstein(de);

			}
		}



	//j is in the outer subset, molA->getNull()[i] can't be true
	}else{

		for(int r=0; r<this->length; r++){



			if(r!=i && r!=j && !molA->getNull()[r]){//

					de = fabs(molA->getDist()[i][r] - molB->getDist()[j][r]);
				    deltaD += de;
				    deltaN += levittGerstein(de);

					de = fabs(molA->getDist()[i][r] - molB->getDist()[i][r]);
				    //deltaD -= de;
				    deltaN -= levittGerstein(de);

			}
		}


	}


}

void InterLig::deltaSeqB(int i, int j){

	double de=0;



	//b[j] and b[i] to be switched. one of them might be a gap (only j, actually. so the b[j] positions shouldn't compared with the a[i] positions, in this case)

	//both in the inner subset
	if(j < this->length){




			if(!molA->getNull()[i]){//if the difference is lower, delta has to be negative   //



				             //new difference   	 //old difference
				de += - (MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[j]-'A']] - MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']]);


			}

			if(!molA->getNull()[j]){//if the difference is lower, delta has to be negative //

				             //new difference   	 //old difference
				de+= - (MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[i]-'A']] - MOL2[AAtoi[molA->getSeq()[j]-'A']][AAtoi[molB->getSeq()[j]-'A']]);


			}




	//j is in the outer subset, molA->getNull()[i] can't be true
	}else{



				        //new difference   	 													//old difference
				de += - (MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[j]-'A']] - MOL2[AAtoi[molA->getSeq()[i]-'A']][AAtoi[molB->getSeq()[i]-'A']]);



	}

	this->deltaS = de;

}

void InterLig::setD0(double d){

	this->d0 = d;

}

double InterLig::getScore(){

	return this->score/this->contributes;

}

double InterLig::getdScore(){

	return this->dScore/this->contributes;

}

double InterLig::getsScore(){

	return this->sScore;

}

double** InterLig::getRotation(){

	return this->s;

}

double** InterLig::getI(){

	return this->I;

}

double* InterLig::getTranslation1(){

	return this->trans1;

}

double* InterLig::getTranslation2(){

	return this->trans2;

}

int InterLig::getContributes(){

	return this->contributes;

}

/*double InterLig::temperature(int i, double T0){

	double a;
	//a=a*a;
	a=T0*exp(-i/2000);
	return a;

}*/

double InterLig::temperatureSD(double mi, double sigma, double Told){

	double Tnew;
	double fdyn;

	//fdyn = 1/(1 + Told*mi*log(1+0.0001)/sigma);
	fdyn = exp(-(0.8*sigma)/(sigma*sigma));
	Tnew=Told*fdyn;

	return Tnew;

}

void InterLig::destroy()
{

	//fix this


	delete[] diffCumulative;
	for(int x=0; x<this->length; x++){

		delete[] diff[x];

	}

	delete[] diff;

}

void InterLig::printDiff(){

	for(int i=0; i<this->length; i++){
		for(int j=0; j<this->length; j++){


			//cout << printf("%1.2f\t", this->diff[i][j]);
			cout << this->diff[i][j] << "\t";
		}
		cout << "\n";
	}
}

double* InterLig::getDiffCumulative(){
	return diffCumulative;
}

void InterLig::loadMOL2(const char *fileA)
{

	int l=24;

    MOL2 = new int*[l];
	AAtoi = new int[l];

	for(int i=0; i<l; i++){
		MOL2[i] = new int[l];
	}
    	string line;
	string word;
	//load PDB coordinates
	ifstream fin (fileA);

    	if (fin.is_open()){
        	for(int m=0; m<l; m++){

			getline(fin, line);
			istringstream iss(line, istringstream::in);
			for(int n=0; n<l+1; n++){

				iss >> word;

				if(n>0){
					MOL2[m][n-1] = atoi(word.c_str());
				}else if(word.c_str()[0] != '-'){

					AAtoi[word.c_str()[0]-'A'] = m;
					//cout << word.c_str()[0] << " " << word.c_str()[0]-'A' << " " << m << "\n" << flush;
				}

			}

        	}
        fin.close();
    }
    else
    {
	char error[5000];
	sprintf(error, "Can not open file: %s\n", fileA);
	cerr << error << endl;
	exit(1);
    }



}

#endif /* INTERLIG_H_ */
