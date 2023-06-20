//2D interpolation
#include <cmath>
#include "polyinterp2.h"
#include <iostream>
#include <fstream>


#include <sys/types.h> //getpid()
#include <unistd.h> //getpid()
#include <stdlib.h> //labs
#include "prng2.h"

long int seed=seedgen();
Ran myran(seed);
Ran * pmyran = &myran;
Normaldev_BM_prng myrang(0.,1.,pmyran);


#include "Hydroinfo_h5.h"
// Read in hydro data
HydroinfoH5 testreadin("data/JetData0010.h5", 500, 0);
fluidCell* fluidCellptr = new fluidCell();


//code from gubtest.cpp
// Useful constants
const double hbarc = 0.197, PI = 3.14159, mc = 1.5, mb = 4.18;
const int Nc = 3; //nd is the number of spatial dimensions
// double maxLi[2][2] = {{-10.,10.},{-10.,10.}}; //NOTE HARD CODING OF BOX'S DIMENSIONS!!!!


// g in terms of alpha
double gYM (double alf) {
	double val;
	val = sqrt(4. * PI * alf);
	return val;
}


// 't Hooft constant in terms of g and Nc
double lam (double g, int Nc) {
	double val;
	val = g * g * Nc;
	return val;
}


// update the local fluid rest frame temperature and the velocity of the local fluid cell in the lab frame at time t and position x
// function is passed a pointer to the fluid temperature (and flow) array; function updates the values the pointer points to
void tempupdate (double t, double x[], double hytemp[], double Tmult, int nd) {
	testreadin.getHydroinfo(t, x[0], x[1], fluidCellptr);
	hytemp[0] = Tmult * (fluidCellptr->temperature);
	hytemp[1] = fluidCellptr->vx;
	hytemp[2] = fluidCellptr->vy;
	return;
}


// energy = sqrt( m^2 + p^2)
inline double energy (double p[], double m, int nd) {
	double en2 = m * m;
	int i;
	for (i = 0; i < nd; ++i) {
		en2 += p[i] * p[i];
	}
	return sqrt(en2);
}

bool isinbox (double x[], double maxLi[][2], int nd){
	int i; //iterators
	for (i = 0; i < nd; ++i){
		if (maxLi[i][0] >= x[i])
			return false;
	}
	for (i = 0; i < nd; ++i){
		if (x[i] >= maxLi[i][1])
			return false;
	}
	return true;
}

inline double mypospow (const double base, const int exp) {
	double res=base;
	for (int i = 1; i < exp; ++i){
		res *= res;
	}
	return res;
}


// void gubdrag (const double xi[], const long int lxi, const double pi[], const long int lpi, const double maxl, const double t0,
// 		const double dt, const double m, const double lamb, const double Tmult, const double maxt, const double minT, const int nd, const int runnumber) {
int gubdrag (const double xi[], const long int lxi, const double pi[], const long int lpi, const double maxl, const double t0,
		const double dt, const double m, const double lamb, const double Tmult, const double maxt, const double minT, const double Epmax, const int nd, double pf[]) {
	//define the position, momentum, and time.  Particle starts at initial position xi with initial momentum pi.
	//We assume particle creation is at t = 0.  The particle propagates unmodified until thermalization, which occurs at t = t0.
	double t = t0, x[nd], p[nd], maxLi[nd][2];
	int i, j; //iterators

	//set the size of the box the particle is in
	for (i = 0; i < nd; ++i){
		maxLi[i][0] = -maxl;
		maxLi[i][1] = maxl;
	}

	//set the position and momenta equal to their initial values
	for (i = 0; i < nd; ++i){
		x[i]=xi[i];
		p[i]=pi[i];
	}


	// double hytemp[nd+1]; //define an array that will hold the local fluid temperature and flow velocity (in the lab frame)
	// tempupdate(t,x,hytemp,Tmult,nd);

	double hytemp[nd+1]; //define an array that will hold the local fluid temperature and flow velocity (in the lab frame)
	// for (i = 0; i < nd + 1; ++i){
	// 	hytemp[i] = ftemp[i];
	// }
	tempupdate(t,x,hytemp,Tmult,nd);

	double Ep=energy(p,m,nd); //energy of the heavy quark given its momentum p and mass m in the lab frame
	
	double mu; //define and initialize the drag
	double pmult; //define the momentum multiplier for the drag process

	double fbeta2, fgamma; //define fluid flow quantities used for boosting

	double Lambda[nd+1][nd+1], invLambda[nd+1][nd+1]; //define the boost and inverse boost matrices

	double boostedp[nd], boostedEp, boosteddt; //define the hq momentum and energy and the timestep in the boosted frame 

	double Cij[nd][nd], kc, cc; //define the diffusion correlation matrix and related quantities
	double pgamma, pgamma12; //these will hold the gamma factor and square root of the particle in the local fluid rest frame, respectively

	//move the particle along without altering its direction or momentum for the first t0 fm
	for (i = 0; i < nd; ++i) {
		x[i] += p[i] * t / Ep;
	}

	double dW[nd]; //define the differential Wiener process; this will hold Normal(0.,1.) distributed pseudo random numbers

	//update the particle position and momentum so long as: the maximum time has not elapsed, the local temperature is above the minimum value
	//and the particle remains inside a defined box. Note hard coding of 2 spatial dimensions in the conditional here.
	while (t < maxt && hytemp[0] > (Tmult * minT) && isinbox(x,maxLi,nd) && Ep < Epmax){
		tempupdate(t,x,hytemp,Tmult,nd); //update the local fluid temp & velocities

		//compute beta^2 = beta.beta, where beta is the velocity of the local fluid cell in the lab frame
		fbeta2 = 0.;
		for (i = 1; i < nd + 1; ++i) {
			fbeta2 += hytemp[i] * hytemp[i];
		}

		fgamma = 1./sqrt(1.-fbeta2); //compute the gamma factor for the local fluid cell

		//load up the boost and inverse boost matrices
		Lambda[0][0] = fgamma;
		invLambda[0][0] = fgamma;

		for (i = 1; i < nd + 1; ++i) {
			Lambda[0][i] = -fgamma * hytemp[i];
			Lambda[i][0] = -fgamma * hytemp[i];
			invLambda[0][i] = fgamma * hytemp[i];
			invLambda[i][0] = fgamma * hytemp[i];
		}

		for (i = 1; i < nd + 1; ++i) {
			for (j = 1; j < nd + 1; ++j) {
				Lambda[i][j] = hytemp[i] * hytemp[j] * fgamma * fgamma / (1. + fgamma);
				if (i == j)
					Lambda[i][j] += 1.;
				invLambda[i][j] = Lambda[i][j];
			}
		}

		//boost the momentum into the local rest frame of the fluid
		for (i = 0; i < nd; ++i) {
			boostedp[i] = Lambda[i+1][0] * Ep;
			for (j = 0; j < nd; ++j) {
				boostedp[i] += Lambda[i+1][j+1] * p[j];
			}
		}

		boosteddt = dt / fgamma; //boost the time step

		boostedEp = energy(boostedp,m,nd); //compute the energy of the hq in the fluid rest frame
		pgamma = boostedEp / m; //gamma factor of the HQ in the local fluid rest frame
		pgamma12 = sqrt(pgamma);

		mu = PI * sqrt(lamb) * hytemp[0] * hytemp[0] / (2. * m); //update the drag coef. for the local fluid cell according to the AdS HQ drag calcs
		kc = PI * sqrt(lamb) * hytemp[0] * hytemp[0] * hytemp[0] / hbarc;
		//update the momentum multiplier; extra terms come from interpreting the Langevin SDE solution as a Stratonovich integral
		pmult = 1. - (mu * boosteddt / hbarc)
			+ (kc * boosteddt * ((5. * pgamma12 * pgamma12 * pgamma12 * pgamma12 * pgamma12 / (4. * boostedEp * boostedEp)) + ((nd - 1.) * pgamma12 / (m * m * (pgamma + 1.)))) / 2.); 
		cc = sqrt(boosteddt * kc * pgamma12) / (m * m * (pgamma + 1.));
		for (i = 0; i < nd; ++i){
			for (j = 0; j < nd; ++j){
				if (i==j)
					Cij[i][j] = cc * (((pgamma + 1.) * m * m) + (boostedp[i] * boostedp[i]));
				else
					Cij[i][j] = cc * boostedp[i] * boostedp[j];
			}
		}

		//load up the normally distributed (mean 0., variance 1.) Wiener random variable
		for (i = 0; i < nd; ++i){
			dW[i] = myrang.dev();
		}
		
		//update the momentum of the hq in the fluid rest frame
		for (i = 0; i < nd; ++i) {
			boostedp[i] *= pmult;
			for (j = 0; j < nd; ++j) {
				boostedp[i] += Cij[i][j] * dW[j];
			}
		}

		boostedEp = energy(boostedp,m,nd); //compute the new energy of the hq in the fluid rest frame

		//boost back to the lab frame
		for (i = 0; i < nd; ++i) {
			p[i] = invLambda[i+1][0] * boostedEp;
			for (j = 0; j < nd; ++j) {
				p[i] += invLambda[i+1][j+1] * boostedp[j];
			}
		}

		Ep = energy(p,m,nd); //compute the hq energy in the lab frame

		//update the position of the hq
		for (i = 0; i < nd; ++i) {
			x[i] += p[i] * dt / Ep;
		}

		t += dt; //update the time
	}
	if (Ep < Epmax) {
		// std::cout.precision(15);
		// for (i = 0; i < nd; ++i){
		// 	std::cout << xi[i] << "\t";
		// }
		// for (i = 0; i < nd; ++i){
		// 	std::cout << pi[i] << "\t";
		// }
		// for (i = 0; i < nd; ++i){
		// 	std::cout << x[i] << "\t";
		// }
		// for (i = 0; i < nd; ++i){
		// 	std::cout << p[i] << "\t";
		// }
		// std::cout << std::endl;
		for (i = 0; i < nd; ++i) {
			pf[i] = p[i];
		}
		return 0;
	}
	else
		return 1;
}

int main () {
	make1dprng pT("data/invbCDF5.dat");
	make1dprng x("data/Pb0010mcprodx.dat");
	make2dprng y("data/Pb0010mcprody.dat");

	int i, j, k, goodreturn, goodreturn2;
	const int jmax = 5000000, jprint = jmax/10, tnd=2;
	double tempx, tempangle, temppT, maxL = 10., t0 = 0.6, tmax = 13., tTf = 0.16, tEpmax = 10000., tTmult = 0.759835686, tm = mb, x00[2] = {0.,0.};
	double tlam, tempxi[2], temppi[2], temppi2[2];

	// paramter must be adjusted to choose which lambda value to use
	const std::string LAMBDA = "G";
	// determining whether to use gubser or "reasonable" lambda
	if(LAMBDA == "R"){
		// reasonable
		std::cout << "'t Hooft's constant determined using \"Reasonable\"\n";
		tlam = gYM(0.3) * gYM(0.3) * Nc;
	} else if (LAMBDA == "G") {
		// gubser
		std::cout << "'t Hooft's constant determined using Gubser\n";
		tlam = 5.5;
	} else {
	// invalid argument
		std::cout << "Invalid Lambda Selected\n";
		return 0;
	}

	double thytemp[tnd+1]; //define an array that will hold the local fluid temperature and flow velocity (in the lab frame)
	tempupdate(t0,x00,thytemp,tTmult,tnd);

	//compute the largest possible drag coefficient and use it to determine the appropriate dt.
	//from code testing we found that a dt <~ mu / 150 was needed to get results to converge to a thermal distribution when the diffusion was given by the FD thm
	const double tmu = PI * sqrt(tlam) * thytemp[0] * thytemp[0] / (2. * tm), tdt = hbarc / (150. * tmu);

	//define variables for storing data produced by the drag/diffusion code
	double pfi[2], pTf, angle; 
	double pfi2[2], pTf2, angle2;

	//set up the bins for the data; make sure they're sensible
	const double startx = 0., endx = 600., dx = .5, tol = 0.01; //x = pT
	if (!(startx<endx)) {
		std::cout << "startx is not less than endx" << std::endl;
		return 1;
	}
	if ( std::abs( ((endx - startx)/dx) - round((endx - startx)/dx)) > tol ){
		std::cout << "dx does not create equal sized bins within tolerance" << tol << std::endl;
	}

	/*
	const int numang = 6;
	int numpts = (int) ((endx - startx)/dx);
	unsigned long int counts [numpts][numang]; //determine the size of the binning array
	//initialize the binning array
	for (i = 0; i < numpts; ++i) {
		for (j = 0; j < numang; ++j){
			counts[i][j] = 0;
		}
	}

	int countsbin, angbin; //declare variables for the correct bin to count the result from the drag code
	double dphi = PI / (2 * numang); //size of the delta phi bins
	*/

	clock_t startTime = clock();

	std::cout << "Seed is " << seed << "\n";

	std::cout.precision(15);
	ofstream resultFile("out/output_"+LAMBDA+"b.csv");
	// write headers to file
	resultFile << "quarkLead_PT,quarkSecondary_PT,deltaAngle\n";

	int trials = 0;
	char esc_char = 27;
	char sep = ',';
	double deltaAngle;
	for (j = 1; j < jmax + 1; ++j){
		// myran.state();
		tempx = x.val(myran.doub());
		tempangle = 2 * PI * myran.doub();
		temppT = pT.val(myran.doub());
		tempxi[0] = tempx;
		tempxi[1] = y.val(myran.doub(),tempx);
		temppi[0] = temppT * cos(tempangle);
		temppi[1] = temppT * sin(tempangle);
		// std::cout << tempxi[0] << "\t" << tempxi[1] << "\t" << temppi[0] << "\t" << temppi[1] << std::endl;
		// for original particle
		goodreturn = gubdrag (tempxi, 1L, temppi, 1L, maxL, t0, tdt, tm, tlam, tTmult, tmax, tTf, tEpmax, tnd, pfi);
		
		// for 180 phase shifted particle
		// alter intial momentum to be 180 out of phase
		for(int z = 0; z < tnd; z++){
			temppi2[z] = -1 * temppi[z];
		}

		goodreturn2 = gubdrag (tempxi, 1L, temppi2, 1L, maxL, t0, tdt, tm, tlam, tTmult, tmax, tTf, tEpmax, tnd, pfi2);
		if ((goodreturn == 0) && (goodreturn2 == 0)){
			// FOR ORIGINAL PARTICLE
			pTf = sqrt( (pfi[0]*pfi[0]) + (pfi[1]*pfi[1]) );
			//countsbin = (int) ceil( (pTf - startx)/dx ) - 1; //determine the bin in the array the data point falls in
			angle = atan2(pfi[1],pfi[0]);
			

			// FOR OPPOSITE PARTICLE
			pTf2 = sqrt( (pfi2[0]*pfi2[0]) + (pfi2[1]*pfi2[1]) );
			angle2 = atan2(pfi2[1],pfi2[0]);

			// Doing some calculations with both particles
			deltaAngle = angle2 - angle;
			// Want angles between 0 and 2*PI
			if (deltaAngle < 0){deltaAngle+=2*PI;}
			//std:cout << deltaAngle << "\n";

			// write output to File
			if (pTf >= pTf2) {
				// first particle is leading
				resultFile << pTf << sep << pTf2 << sep << deltaAngle << "\n";
			} else {
				// second particle is leading
				resultFile << pTf2 << sep << pTf << sep << deltaAngle << "\n";
			}
			
			/*
			//determine the angular bin the data point falls in
			i = 0;
			while (! ((i * dphi) - PI <= angle && angle < ((i + 1) * dphi) - PI)) {++i;}
			// std::cout << "i = " << i << std::endl;
			if ((-PI <= angle && angle < -PI/2.) || (0. <= angle && angle < PI/2.)){
				angbin = i % numang;
			}
			else {
				angbin = numang - 1 - (i % numang);
			}
			//count the point
			if (0 <= countsbin && countsbin < numpts){ //make sure the bin is within the array
				++counts[countsbin][angbin]; //add the point to the bin count
			}
			//print out the results so far if the particle number is a multiple of jprint, currently set at 5 * 10^6
			if ((j % jprint) == 0){
			//	print out the binned counts
				for (i = 0; i < (endx - startx)/dx; ++i){
					std::cout << startx + ((i + 0.5) * dx) << "\t";
					for (k = 0; k < numang; ++k) {
						std::cout << counts[i][k] << "\t";
					}
					std::cout << std::endl;
				}
			
				
			}
			*/
			// prints how many trials have occured thus far
			if ((j % jprint) == 0){
				trials += 10;
				std::cout << esc_char << "[1m" << trials << "%" << esc_char << "[0m\n";
			}


			if (j == jmax) {
				//print out the end state of the prng generator
				std::cout << "Done!" << "\n";
				std::cout << "The state of the generator is" << "\n";
				myran.state();
				std::cout << "Runtime of " << double(clock() - startTime ) / double(CLOCKS_PER_SEC) << " seconds.\n";
			}
			
		}
	}
	resultFile.close();

//	print out the binned counts if the final number of particles is not a multiple of jprint (since the final result would already have been printed in this case)
	/*
	if ((jmax % jprint) != 0){
		for (i = 0; i < (endx - startx)/dx; ++i){
			std::cout << startx + ((i + 0.5) * dx) << "\t";
			for (j = 0; j < numang; ++j) {
				std::cout << counts[i][j] << "\t";
			}
			std::cout << std::endl;
		}
		//print out the end state of the prng generator
		std::cout << "The state of the generator is" << std::endl;
		myran.state();
		std::cout << "Runtime of " << double(clock() - startTime ) / double(CLOCKS_PER_SEC) << " seconds." << std::endl;
	}
	*/
	return 0;
}