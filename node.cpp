#include "node.h"
//define a ternery tree structure
int node::depth;

void node::insert(NodeType nodeType, VectorXd x, double y, VectorXd d){
	//first creat a node
	node* p = new node(x, y, d);
	//now link it to the parent (this)
	if (nodeType == NodeType::LEFT){
		this->left = p;
		this->left->h = this->h + 1;
		this->left->j = 3 * this->j;
		this->left->parent = this;

	}
  if (nodeType == NodeType::MIDDLE){
		this->middle = p;
		this->middle->h = this->h + 1;
		this->middle->j = 3 * this->j + 1;
		this->middle->parent = this;
	}
  if (nodeType == NodeType::RIGHT){
		this->right = p;
		this->right->h = this->h + 1;
		this->right->j = 3 * this->j + 2;
		this->right->parent = this;
	}
	depth = max(depth, h + 1);

}