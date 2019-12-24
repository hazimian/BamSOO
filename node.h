#include "general.h"



class node{
public:

  node(VectorXd x, pair<double, double> y, VectorXd d, bool Evaltd);
	node* left;
	node* right;
	node* middle;
	node* parent;

	VectorXd x;
	pair<double,double> y;
	int h;
	int j;
	bool Evaltd;
	VectorXd d;
	static void set_depth(int d){
		depth = d;
	};
	static int get_depth(){
		return depth;
	};
  void insert(char n[], VectorXd x, pair<double, double> y, VectorXd d, bool Evaltd);

private:
	static int depth;
};