
#include "SOO.h"



VectorXd shift(VectorXd &d){
	VectorXd z = VectorXd::Zero(d.size());
	VectorXd::Index k = 0;
	d.maxCoeff(&k);
	z[k] = 1;
	d[k] = d[k] / 3.0;
	return z;

}

void fopt(node* root, node* &mn, int h, double &m,bool flag){

	if (!root || root->h > h)
		return;
	if (root->h == h){//consider only if has no children
		if (!root->left && !root->middle & !root->right){
			if (root->y.first > m){
				m = root->y.first;
				mn = root;
			}
			return;
		}
		else
			return;
	}
	//otherwise
	if (root->left)
		fopt(root->left, mn, h, m,flag);
	if (root->middle)
		fopt(root->middle, mn, h, m,flag);
	if (root->right)
		fopt(root->right, mn, h, m,flag);



}



VectorXd denorm(VectorXd x, VectorXd lb, VectorXd ub){
	return .5*(ub - lb).cwiseProduct(x) + .5*(lb + ub);
}


double SOO(function<pair<double,double>(VectorXd)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax){
	auto start = chrono::high_resolution_clock::now();
	double ybest = nINF;
  	//create the root node
	VectorXd mid = VectorXd::Zero(lb.rows());//normalized to [-1,1]
	VectorXd delta = 2 * VectorXd::Ones(lb.rows());

  node* root = new node(mid, fun(denorm(mid, lb, ub)), delta, true);
	ybest = max(ybest, root->y.first);


	//now create level h=1
	VectorXd c = shift(root->d);// d passed by reference
	//now insert into h=1
	root->insert("left", root->x - c.cwiseProduct(root->d), fun(denorm(root->x - c.cwiseProduct(root->d), lb, ub)), root->d,true);
	ybest = max(ybest, root->left->y.first);
	root->insert("middle", root->x, root->y, root->d,root->Evaltd);
	root->insert("right", root->x + c.cwiseProduct(root->d), fun(denorm(root->x + c.cwiseProduct(root->d), lb, ub)), root->d,true);
	ybest = max(ybest, root->right->y.first);

	double vmax;
	int count = 3;
	double m;
	node* mn = nullptr;
	int n = 2;
	int depth = 1;
	while (count < pmax){
		vmax = nINF;
		depth = node::get_depth();//update depth
		for (int h = 0; h <= depth; ++h){
			m = nINF;//initialize 
			mn = nullptr;
			fopt(root, mn, h, m,true);//find the optimum to expand
			if (mn && mn->y.first>vmax){
				//now generate children
				n++;
				vmax = mn->y.first;
				c = shift(mn->d);
				
				mn->insert("middle", mn->x, mn->y, mn->d,mn->Evaltd);//just copy parent into middle
				mn->insert("left", mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d),lb,ub)), mn->d,true);//create left node
				count++;
				ybest = max(ybest, mn->left->y.first);
				cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
				mn->insert("right", mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d),lb,ub)), mn->d,true);//create right node
				count++;
				ybest = max(ybest, mn->right->y.first);
				cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;

			}
		}
	}
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	cout <<"elapsed time: "<< elpased.count() << "seconds" << endl;
	return ybest;
}


double LOGO(function<pair<double,double>(VectorXd)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax){
	auto start = chrono::high_resolution_clock::now();
	double ybest = nINF;
	double ybest_ = nINF;

	//create the root node
	//VectorXd mid = 0.5*(lb + ub);
	//VectorXd delta = (ub - lb);
	VectorXd mid = VectorXd::Zero(lb.rows());//normalized to [-1,1]
	VectorXd delta = 2 * VectorXd::Ones(lb.rows());
	node* root = new node(mid, fun(denorm(mid, lb, ub)), delta, true);
	ybest = max(ybest, root->y.first);


	//now create level h=1

	VectorXd c = shift(root->d);// d passed by reference
	//now insert into h=1
	root->insert("left", root->x - c.cwiseProduct(root->d), fun(denorm(root->x - c.cwiseProduct(root->d), lb, ub)), root->d, true);
	ybest = max(ybest, root->left->y.first);
	root->insert("middle", root->x, root->y, root->d, root->Evaltd);
	root->insert("right", root->x + c.cwiseProduct(root->d), fun(denorm(root->x + c.cwiseProduct(root->d), lb, ub)), root->d, true);
	ybest = max(ybest, root->right->y.first);
	int jw = 0;
	double vmax;
	int count = 3;
	double m;
	int n = 2;
	node* mn = nullptr;
	int depth = 1;
	int hplus = 0;
	int W[] = {3,4,5,6,8,30};
	int w = W[0];
	int hupper = 0;
	while (count < pmax){
		vmax = nINF;
		depth = node::get_depth();//update depth
		hplus = hupper;
		for (auto k = 0; k <= max(floor(min(double(w*(sqrt(n)-1)),double(hupper))/w),double(hplus)); ++k){
			m = nINF;//initialize 
			mn = nullptr;
			for (int u = k*w; u<=min(k*w+w-1, depth);u++){
				fopt(root, mn, u, m,true);//find the optimum to expand
			}
			if (mn && mn->y.first>vmax){
				//now generate children
				n++;
				hplus = 0;
				hupper = max(hupper, mn->h + 1);
				vmax = mn->y.first;
				c = shift(mn->d);
				mn->insert("middle", mn->x, mn->y, mn->d, mn->Evaltd);//just copy parent into middle
				mn->insert("left", mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create left node
				count++;
				ybest_ = ybest;
				ybest = max(ybest, mn->left->y.first);
				if (ybest >= ybest_)
					w = W[min(jw, 5)];
				else
					w = W[max(jw-2, 0)];
				cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
				mn->insert("right", mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create right node
				count++;
				ybest_ = ybest;
				ybest = max(ybest, mn->right->y.first);
				if (ybest >= ybest_)
					w = W[min(jw, 5)];
				else
					w = W[max(jw - 2, 0)];
				cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;

			}
		}
	}
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	cout << "elapsed time: " << elpased.count() << "seconds" << endl;
	return ybest;
}


