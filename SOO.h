

#include "node.h"

VectorXd shift(VectorXd&);

double SOO(function<pair<double,double>(VectorXd)> , VectorXd , VectorXd , double, int);
VectorXd shift(VectorXd &);
void fopt(node*, node*&, int, double&,bool);
VectorXd denorm(VectorXd, VectorXd, VectorXd);
double LOGO(function<double(VectorXd)>, VectorXd, VectorXd, double, int);