
#include "SOO.h"
double BamSOO(function<pair<double, double>(VectorXd)>, VectorXd, VectorXd, double, int);
struct GP{
	double mu;
	double sig2;
};

double kernel(VectorXd , VectorXd , double );
GP FitGP(const MatrixXd* , const VectorXd* , const VectorXd, double );