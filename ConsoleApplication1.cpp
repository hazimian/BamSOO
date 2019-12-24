// ConsoleApplication1.cpp : Defines the entry point for the console application.
//


#include "IMGPO.h"



pair<double,double> branin(VectorXd x){
	const double a = 1.0;
	const double b = 5.1 / (4.0*pow(PI, 2));
	const double c = 5.0 / PI;
	const double r = 6;
	const double s = 10;
	const double t = 1.0 / (8.0*PI);


	return make_pair(-(a*pow((x[1] - b*pow(x[0], 2) + c*x[0] - r), 2) + s*(1 - t)*cos(x[0]) + s),x(0)-x(1));

}

void PrintNode(node node){
	cout << "---node info:----" << endl;
	cout << "x= " << node.x.transpose() << endl;
	cout << "y= " << node.y.first << endl;
	cout << "d= " << node.d.transpose() << endl;
	cout << "h= " << node.h << endl;
	cout << "j= " << node.j << endl;
	cout << "depth= " << node.get_depth() << endl;
	if (node.left)
		PrintNode(*node.left);
	else
		cout << "left: " << node.left << endl;
	if (node.middle)
		PrintNode(*node.middle);
	else
		cout << "middle: " << node.middle << endl;
	if (node.right)
		PrintNode(*node.right);
	else
		cout << "right: " << node.right << endl;

}
//



double hart3(VectorXd x){
	VectorXd alpha(4);
	alpha << 1.0, 1.2, 3.0, 3.2;
	MatrixXd A(4,3);
	A << 3.0, 10, 30,
		0.1, 10, 35,
		3.0, 10, 30,
		0.1, 10, 35;
	MatrixXd P(4,3);
	P << 3689, 1170, 2673,
		4699, 4387, 7470,
		1091, 8732, 5547,
		381, 5743, 8828;
	P = 1e-4 *P;

	double outer = 0;
	double inner = 0;
	for (auto i : { 0, 1, 2, 3 }){
		inner = 0;
		for (auto j : { 0, 1, 2 })
			inner += A(i, j) * pow(x(j) - P(i, j), 2.0);
		outer = outer + alpha[i] * exp(-inner);
	}
	return outer;
}

double Rosen(VectorXd x){
	int d = x.size();
	double sum = 0;

	for (int i = 0; i < d - 1; i++){
		sum = sum + 100 * pow(x(i+1) - pow(x(i), 2),2) + pow(x(i) - 1,2);
	}

	return -sum;
}

double hart6(VectorXd x){
	VectorXd alpha(4);
	alpha << 1.0, 1.2, 3.0, 3.2;
	MatrixXd A(4, 6);
		A << 10, 3, 17, 3.5, 1.7, 8,
		0.05, 10, 17, 0.1, 8, 14,
		3, 3.5, 1.7, 10, 17, 8,
		17, 8, 0.05, 10, 0.1, 14;
	MatrixXd P(4, 6);
	P << 1312, 1696, 5569, 124, 8283, 5886,
		2329, 4135, 8307, 3736, 1004, 9991,
		2348, 1451, 3522, 2883, 3047, 6650,
		4047, 8828, 8732, 5743, 1091, 381;
	P = 1e-4 *P;
	double outer = 0;
	double inner = 0;
	for (auto i : { 0, 1, 2, 3 }){
		inner = 0;
		for (auto j : { 0, 1,2,3,4,5 })
			inner += A(i, j) * pow(x(j) - P(i, j), 2.0);
		outer = outer + alpha[i] * exp(-inner);
	}
	return outer;

	return (2.58 + outer) / 1.94;
}
double Rast(VectorXd x){
	x(0) = x(0)-.2;
	x(1) = x(1)+.5;
	x(2) = x(2)-.8 ;
	double z = -(30 + pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2) - 10 * (cos(1 * PI*x(0)) + cos(1 * PI*x(1)) + cos(1 * PI*x(2))));
	return z;
}

int node::depth = 0;
int _tmain(int argc, _TCHAR* argv[])
{
	//function attributes

	VectorXd lb(2), ub(2);
	lb << -5, 0;
	ub << 10, 15;
	function<pair<double,double>(VectorXd)> fun = branin;
	const double yopt = -.397887;
  

	
	//VectorXd lb(3), ub(3);
	//lb << 0, 0, 0;
	//ub << 1, 1, 1;
	//function<double(VectorXd)> fun = hart3;
	//const double yopt = 3.86278;

	//VectorXd lb(3), ub(3);
	//lb << -2, -2, -2;
	//ub << 2, 2, 2;
	//function<double(VectorXd)> fun = Rast;
	//const double yopt = 0;

	//VectorXd lb(6), ub(6);
	//lb << 0, 0, 0,0,0,0;
	//ub << 1, 1, 1,1,1,1;
	//function<double(VectorXd)> fun = hart6;
	//const double yopt = 3.32237;


	//VectorXd lb(2), ub(2);
	//lb << -5, -5;
	//ub << 5, 5;
	//function<double(VectorXd)> fun = Rosen;
	//const double yopt = 0;

	const int pmax = 200;
	double ybest;

	ybest=SOO(fun, lb, ub,yopt,pmax);

	

	getchar();
	return 0;
}

