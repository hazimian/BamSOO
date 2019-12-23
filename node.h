#include "general.h"
#pragma once

enum class NodeType
{
  LEFT,
  RIGHT,
  MIDDLE
};

class node{
public:

	node(VectorXd x, double y, VectorXd d) : x(x), y(y), d(d), parent(nullptr), left(nullptr), middle(nullptr), right(nullptr), h(0){  }
	node* left;
	node* right;
	node* middle;
	node* parent;
	VectorXd x;
	double y;
	int h;
	int j;
	VectorXd d;
	static void set_depth(int d){
		depth = d;
	};
	static int get_depth(){
		return depth;
	};
	void insert(NodeType nodeType, VectorXd x, double y, VectorXd d);

private:
	static int depth;
};