#include "IMGPO.h"


double IMGPO(function<pair<double,double>(VectorXd)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax){
	auto start = chrono::high_resolution_clock::now();
	double ybest = nINF;
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


	int count = 3;
	double m;
	node* mn = nullptr;
	node* mm = nullptr;
	int depth = 1;
	double B;
	GP GPm;
	MatrixXd *X = new MatrixXd(lb.size(), 3);
	VectorXd *Y = new VectorXd(3);
	(*X) << root->x, root->left->x, root->right->x;
	(*Y) << root->y.first, root->left->y.first, root->right->y.first;
	MatrixXd temp;
	double UCB = 0;
	double vmax;



	int jstar = 0;
	unsigned int zeta = 0;
	int N = 1;
	unsigned int E = 1;

	while (count < pmax){
		vmax = nINF;
		depth = node::get_depth();//update depth
		for (int h = 0; h <= depth; ++h){

			m = nINF;//initialize 
			mn = nullptr;
			fopt(root, mn, h, m,true);//find the optimum to expand
			//check ii
			if (mn && mn->y.first>vmax){

				if (mn->Evaltd){
					vmax = mn->y.first;
				}
				else{
					mn->y = fun(denorm(mn->x, lb, ub));
					mn->Evaltd = true;
					ybest = max(ybest, mn->left->y.first);
					X->conservativeResize(X->rows(), X->cols() + 1);
					Y->conservativeResize(Y->rows() + 1);
					X->col(X->cols() - 1) = mn->x;
					(*Y)(Y->rows() - 1) = mn->y.first;
				}

				//find maximum of smaller intervals (can be further optimized)
				m = nINF;//initialize 
				mm = nullptr;
				for (int hh = h + 1; hh <= depth; hh++){
					fopt(root, mm, hh, m,true);//find the optimum to expand
				}
				//evaluate UCB
				B = sqrt(2 * log2(pow(PI*N, 2) / (12 * .05)));
				GPm = FitGP(X, Y, mn->x, 0);//fit GP
				N++;
				UCB = GPm.mu + B*sqrt(abs(GPm.sig2));
				//check iii
				if (UCB > m){
					//start dividing and insertion
					c = shift(mn->d);
					mn->insert("middle", mn->x, mn->y, mn->d, mn->Evaltd);//just copy parent into middle
					//check UCB
					B = sqrt(2 * log2(pow(PI*N, 2) / (12 * .05)));
					GPm = FitGP(X, Y, mn->x - c.cwiseProduct(mn->d), 0);//fit GP
					N++;
					UCB = GPm.mu + B*sqrt(abs(GPm.sig2));
					//evaluate if greater
					if (UCB > ybest){
						mn->insert("left", mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create left node
						count++;
						ybest = max(ybest, mn->left->y.first);
						cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
					}
					else
						mn->insert("left", mn->x - c.cwiseProduct(mn->d), make_pair(UCB,-1), mn->d, true);//create left node

					//check UCB
					B = sqrt(2 * log2(pow(PI*N, 2) / (12 * .05)));
					GPm = FitGP(X, Y, mn->x + c.cwiseProduct(mn->d), 0);//fit GP
					N++;
					UCB = GPm.mu + B*sqrt(abs(GPm.sig2));
					if (UCB > ybest){
						mn->insert("right", mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create right node
						count++;
						ybest = max(ybest, mn->right->y.first);
						cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
					}
					else
						mn->insert("right", mn->x - c.cwiseProduct(mn->d), make_pair(UCB,-1), mn->d, true);//create left node

				}

			}


		}

	}
	delete X, Y;
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	cout << "elapsed time: " << elpased.count() << "seconds" << endl;
	return ybest;
}