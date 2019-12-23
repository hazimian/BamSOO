#include "general.h"
#include "node.h"
#include <chrono>
#include <vector>
#include "BO4.h"

VectorXd shift(VectorXd&);

void fopt(node*, node**, int, double&);

double SOO(function<double(VectorXd, void* dataPtr)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax, void *dataPtr, VectorXd &xbest, const std::shared_ptr<BO4Input> mInp);

double BamSOO(function<double(VectorXd)> , VectorXd , VectorXd, double, int);

double LOGO(function<double(VectorXd, void* dataPtr)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax, void *dataPtr, VectorXd &xbest, const std::shared_ptr<BO4Input> mInp);