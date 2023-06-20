
struct Ran {
	//program (slightly altered) from Numerical Recipes 3rd Edition
private:
	unsigned long long int u,v,w;
public:
	Ran(unsigned long long int j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline unsigned long long int int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline unsigned int int32() { return (unsigned int) int64(); }
	//my own addition
	inline void state() {std::cout << "(u, v, w) = " << "(" << u << ", " << v << ", " << w << ")" << std::endl; return;}
	inline void setstate(unsigned long long int uin, unsigned long long int vin, unsigned long long int win){
		u=uin;
		v=vin;
		w=win;
		return;
	} 
};

struct Normaldev_BM : Ran {
private:
	double mu,sig;
	double storedval;
public:
	Normaldev_BM(double mmu, double ssig, unsigned long long int i) : Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
	double dev() {
		double v1,v2,rsq,fac;
		if (storedval == 0.) {
			do {
				v1=2.0*doub()-1.0;
				v2=2.0*doub()-1.0;
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=sqrt(-2.0*log(rsq)/rsq);
			storedval = v1*fac;
			return mu + sig*v2*fac;
		}
		else {
			fac = storedval;
			storedval = 0.;
			return mu + sig*fac;
		}
	}
};


struct Normaldev_BM_prng {
	private:
		double mu,sig;
		double storedval;
		Ran * myran;
	public:
		Normaldev_BM_prng(double mmu, double ssig, Ran * pmyran) : mu(mmu), sig(ssig), storedval(0.), myran(pmyran) {}
		double dev() {
			double v1,v2,rsq,fac;
			if (storedval == 0.) {
				do {
					v1=2.0*myran->doub()-1.0;
					v2=2.0*myran->doub()-1.0;
					rsq=v1*v1+v2*v2;
				} while (rsq >= 1.0 || rsq == 0.0);
				fac=sqrt(-2.0*log(rsq)/rsq);
				storedval = v1*fac;
				return mu + sig*v2*fac;
			}
			else {
				fac = storedval;
				storedval = 0.;
				return mu + sig*fac;
		}
	}
};

long seedgen(void) {
	//seed generation from 1005.4117, compatible with parallelization on multiple processors for the generation of different seeds despite very similar initial run times.  
	//note that this method of computing random numbers for parallel processes is discouraged in 0905.4238 (e.g. there might be sequence overlap)
	//code from 1005.4117
	long s, seed, pid;

	pid = getpid(); //get process PID
	s = time (NULL); //get CPU seconds from 01/01/1970

	seed = labs(((s*181)*((pid-83)*359))%104729);
	return seed;
}

//create a non-uniform pseudo random number generator from a specific input inverse CDF specified by filename fn.  the inverse CDF must be of the form of
//a column of double reals; the first n values are the array of points x_i at which the inverse CDF is evaluated, and the second n points are the evaluations
//of the inverse CDF
struct make1dprng {
	private:
		int i, numpts;
		double temp, *x, *y;
		Poly_interp * pinvcdf;
	public:
		make1dprng (const char * fn) {
			//read in the inverse CDF data
			std::ifstream qproddat (fn); //open the file
			//determine the length of the file
			i = 0;
			while (!qproddat.eof()){
				qproddat >> temp;
				++i;
			}
			numpts = i/2; //determine the number of x and y values to be read in
			qproddat.seekg(0);
			x = new double [numpts]; //make an appropriately sized array to hold the x values
			y = new double [numpts]; //make an appropriately sized array to hold the y values
			for (i = 0; i < numpts; ++i) {
				qproddat >> x[i]; //read in x
			}
			for (i = 0; i < numpts; ++i) {
				qproddat >> y[i]; //read in y
			}
			// for (i = 0; i < numpts; ++i){
			// 	std::cout << x[i] << "\t" << y[i] << std::endl;
			// // }
			// Poly_interp invcdf(x, numpts, y, 4); //create a 3rd order polynomial interpolating function of the inverse CDF data
			// pinvcdf = &invcdf;
//			Poly_interp invcdf(x, numpts, y, 4); //create a 3rd order polynomial interpolating function of the inverse CDF data
			pinvcdf = new Poly_interp(x, numpts, y, 4); //create a 3rd order polynomial interpolating function of the inverse CDF data;
		}

		//feed val a random number uniformly distributed between 0. and 1. and it will return a val
		//with the val's distributed according to the inverse CDF input by the constructor
		double val (double prn) {return pinvcdf->interp(prn);}
};

//create an object that will return a non-uniform pseudo random number for a given x and a prng between 0. and 1. with the distribution of prn's determined
//by the input file
	//input file is assumed to be of the following form: a single column of numbers, the first two are integers, the rest are double reals.
	//The first integer is nx, the number of x points in the grid
	//The second integer is ny, the number of y points in the grid
	//The first nx values are the x1 positions
	//The following ny values are the x2 positions
	//The last nx * ny values are the values of the function y(x1,x2).  
	//The matrix y is scanned along for nx values of x1, then the next set of nx x1 values for the next step in x2
struct make2dprng {
	private:
		int i, j, nx, ny;
		double *x1, *x2, *y;
		Bilin_interp * pinvCDF;
	public:
		make2dprng (const char * fn) {
			std::ifstream input (fn); //open the file
			input >> nx; //grab the number of x1 points
			input >> ny; //grab the number of x2 points

			x1 = new double [nx];
			x2 = new double [ny];
			y = new double [nx * ny];

			for (i = 0; i < nx; ++i) {
				input >> x1[i]; //fill an array with the x1 points
			}
			for (i = 0; i < ny; ++i) {
				input >> x2[i]; //fill an array with the x2 points
			}
			for (j = 0; j < ny; ++j) {
				for (i = 0; i < nx; ++i) {
					input >> y[(nx * j) + i]; //fill the matrix y[x1,x2]
				}
			}

			pinvCDF = new Bilin_interp(x1,nx,x2,ny,y); //linearly interpolate between calculated values

		}

		//feed val a pseudo random number between 0. and 1. and a given value of x, and val returns a non-uniformly distributed pseudo
		//random number with distribution given by the input file fn.
		double val (double prn, double x) {return pinvCDF->interp(prn,x);}
};