//#include <cmath>


inline int mymin (int a, int b) {
	return (a < b) ? a : b;
}

inline int mymax (int a, int b) {
	return (a > b) ? a : b;
}

//polynomial interpolation of data based on Numerical Recipes 3rd edition
struct Base_interp {
	protected:
		int n, mm, jsav, cor, dj;
		double *xx, *yy;
	
	public:
		Base_interp( double *x, int numx, double *y, int m) : n(numx), mm(m), jsav(0), cor(0), xx(x), yy(y) {}

		int locate(const double x);

		double virtual rawinterp(int jlo, double x) = 0;

		double interp (double x) {
			int jlo = locate(x);
			return rawinterp(jlo,x);
		}

};

int Base_interp::locate(const double x) {
	int ju,jm,jl;
	if (n < 2 || mm < 2 || mm > n) throw("locate size error");
	bool ascnd=(xx[n-1] >= xx[0]); //is list of x's ascending or descending: it is assumed list x entries are monotonic
	jl=0; //initialize lower limit
	ju=n-1; //initialize upper limit
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1; //integer average of ju and jl (rounded down)
		if (x >= xx[jm] == ascnd) 
			jl = jm;
		else
			ju = jm;
	}
	return mymax(0,mymin(n-mm,jl-((mm-2)>>1)));
}

struct Linear_interp : Base_interp {
	Linear_interp( double * xv, int numx, double * yv) : Base_interp(xv, numx, yv, 2) {}
	double rawinterp(int j, double x) {
		if (xx[j]==xx[j+1]) //table is defective as the points are not monotonically increasing
			return yy[j]; //however we still can return a value for the interpolation
		else
			return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
	}
};

struct Poly_interp : Base_interp {
	private:
		double dy;
	public:
		Poly_interp( double * xv, int numx, double * yv, int m) : Base_interp(xv, numx, yv, m), dy(0.) {}
		double rawinterp(int jl, double x);
};

double Poly_interp::rawinterp(int jl, double x) {
	int i,m,ns=0;
	double y,den,dif,dift,ho,hp,w;
	const double *xa = &xx[jl], *ya = &yy[jl];
	double c[mm], d[mm];

	dif=std::abs(x-xa[0]);
	for (i=0; i<mm; ++i) {
		if ((dift=std::abs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	y=ya[ns--];

	for (m=1; m<mm; ++m) {
		for (i=0; i<mm-m; ++i) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) throw("Poly_interp error");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	}

	return y;
}


//Bilin_interp linearly interpolates values for a function y(x1,x2)
//When called, Bilin_interp needs to be fed the list of x1 points and the list of x2 points that form a regular, 2D grid for which the y values are also given.
//Bilin_interp also needs to know the number of x1 points and the number of x2 points in those respective lists.
//y is assumed to be a 1D array holding the 2D matrix of y(x1,x2) values.  
//y[0] = y(x1 = 0, x2 = 0), y[numx1] = y (x1max, x2 = 0), ..., y[numx1 * numx2] = y(x1max, x2max)
struct Bilin_interp {
	private:
		int numx1, numx2;
		double *x1, *x2, *y;
		Linear_interp x1terp, x2terp;
	public:
		// fill out the x1 and x2 vectors and the matrix of y values.  x1terp and x2terp are dummy objects (in the sense that we are not interested in their y values
		//but only their locate functions to find the right grid points for the interpolation)
		Bilin_interp ( double *xx1, int nx1, double *xx2, int nx2, double *yy) : x1(xx1), numx1(nx1), x2(xx2), numx2(nx2), y(yy), x1terp(x1, numx1, x1), x2terp(x2, numx2, x2) {}

		double interp ( double x1p, double x2p) {
			int i = x1terp.locate(x1p), j = x2terp.locate(x2p);
			double t, u;

			t = (x1p - x1[i]) / (x1[i+1] - x1[i]);
			u = (x2p - x2[j]) / (x2[j+1] - x2[j]);

			return ((1. - t) * (1. - u) * y[(numx2 * j) + i]) + (t * (1. - u) * y[(numx2 * j) + i + 1]) + 
				((1. - t) * u * y[(numx2 * (j + 1)) + i]) + (t * u * y[(numx2 * (j + 1)) + i + 1]);
		}
};

// int main () {
// 	int i, numpts=10;
// 	double xi, yi, dxi=0.01;
// 	double x[numpts], y[numpts];

// 	for (i = 0; i < numpts; ++i) {
// 		x[i] = i + 0.5 * sin (i);
// 		y[i] = i + cos (i * i);
// 	}

// 	Poly_interp myfunc(x, numpts, y, 4);

// 	for (xi = x[0]; xi < x[numpts-1]; xi += dxi){
// 		yi = myfunc.interp(xi);
// 		std::cout << xi << "\t" << yi << std::endl;
// 	}

// 	return 0;
// }