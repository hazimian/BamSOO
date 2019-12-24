

#include "node.h"
//define a ternery tree structure

node::node(VectorXd x, pair<double,double> y, VectorXd d, bool Evaltd) : x(x), y(y), d(d), h(0),Evaltd(Evaltd), parent(nullptr), left(nullptr), middle(nullptr), right(nullptr){  }

void node::insert(char n[], VectorXd x, pair<double, double> y, VectorXd d, bool Evaltd){
	//first creat a node
	node* p = new node(x, y, d,Evaltd);
	//now link it to the parent (this)
	if (!strcmp(n, "left")){
		this->left = p;
		this->left->h = this->h + 1;
		this->left->j = 3 * this->j;
		this->left->parent = this;

	}
	if (!strcmp(n, "middle")){
		this->middle = p;
		this->middle->h = this->h + 1;
		this->middle->j = 3 * this->j + 1;
		this->middle->parent = this;

	}
	if (!strcmp(n, "right")){
		this->right = p;
		this->right->h = this->h + 1;
		this->right->j = 3 * this->j + 2;
		this->right->parent = this;

	}
	depth = max(depth, h + 1);

}