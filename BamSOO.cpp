#include "BamSOO.h"



double kernel(VectorXd x1, VectorXd x2, double p){
	VectorXd V = x1 - x2;
	//return pow(10, p)*exp(-V.transpose()*V);
	double z = sqrt(5 * double(V.transpose()*V) / .25);
	return pow(10, p)*(1 + z + pow(z, 2) / 3)*exp(-z);

}

GP FitGP(const MatrixXd* x, const VectorXd* y, const VectorXd xs, double p){
	MatrixXd Kmat = MatrixXd::Identity(x->outerSize(), x->outerSize());
	MatrixXd iKmat = MatrixXd::Identity(x->outerSize(), x->outerSize());
	MatrixXd Kmat_ = Kmat;
	VectorXd kvec = VectorXd::Zero(x->outerSize());
	vector<double> pvec = { -2,-1.5,-1,-.7,-.5,-.2,0,.5,1};
	//for (int m = 0; m < pvec.size(); m++){
		for (int i = 0; i < x->outerSize(); i++){
			kvec(i) = kernel(x->col(i), xs, p);//creates kvec
			for (int j = 0; j <= i; j++){
				Kmat(i, j) = kernel(x->col(i), x->col(j), p);//creat Kmat
				Kmat(j,i) = Kmat(i, j);
			}
		}

		iKmat = Kmat.llt().solve(MatrixXd::Identity(Kmat.rows(), Kmat.cols()));
		//Kmat_ = Kmat + 1e-8*MatrixXd::Identity(Kmat.innerSize(), Kmat.innerSize());
		//iKmat = (Kmat_).inverse();
		//MatrixXd L(Kmat_.llt().matrixL());
		//LogL = -y->transpose()*iKmat*(*y) - 2 * L.diagonal().array().log().sum();
		//if (isinf(LogL)){
		//	cout << "inf found" << endl;
		//}
		//if (LogL>LogL0 ){
		//	p0 = pvec[m];
		//	iKmat0 = iKmat;
		//	kvec0 = kvec;
		//	LogL0 = LogL;
		//}

	//}
	//cout << p0 << endl;
	double mu = kvec.transpose()*iKmat*(*y);
	double sig2 = kernel(xs, xs, p) - kvec.transpose()*iKmat*kvec;
	return{ mu, sig2 };
}


double BamSOO(function<pair<double,double>(VectorXd)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax){
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
	int N = 1;
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
	while (count < pmax){
		vmax = nINF;
		depth = node::get_depth();//update depth
		for (int h = 1; h <= depth; ++h){
			m = nINF;//initialize 
			mn = nullptr;
			fopt(root, mn, h, m,true);//find the optimum to expand
			if (mn && mn->y.first>vmax){
				//now generate children
				c = shift(mn->d);
				//if selected node has temp value then evaluate
				if (mn->Evaltd)
					mn->insert("middle", mn->x, mn->y, mn->d, mn->Evaltd);//just copy parent into middle
				else
					mn->insert("middle", mn->x, fun(denorm(mn->x, lb, ub)), mn->d, true);

				//do for left element
				GPm = FitGP(X, Y, mn->x - c.cwiseProduct(mn->d), 0);//fit GP
				B = sqrt(2 * log2(pow(PI*N, 2) / (12 * .05)));
				N += 1;
				if (GPm.mu + B*sqrt(abs(GPm.sig2)) > ybest){//if greater than UCB evaluate
					mn->insert("left", mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create left node
					count++;
					ybest = max(ybest, mn->left->y.first);

					X->conservativeResize(X->rows(), X->cols() + 1);
					Y->conservativeResize(Y->rows() + 1);
					X->col(X->cols() - 1) = mn->left->x;
					(*Y)(Y->rows() - 1) = mn->left->y.first;


				}
				else//otherwise use LCB
					mn->insert("left", mn->x - c.cwiseProduct(mn->d), make_pair(GPm.mu - B*sqrt(abs(GPm.sig2)),-1), mn->d, false);//create left node
				//repeat for right element
				B = sqrt(2 * log2(pow(PI*N, 2) / (12 * .05)));
				N += 1;
				GPm = FitGP(X, Y, mn->x + c.cwiseProduct(mn->d), 0);//fit GP
				N += 1;
				if (GPm.mu + B*sqrt(abs(GPm.sig2)) > ybest){//if greater than UCB evaluate
					mn->insert("right", mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d), lb, ub)), mn->d, true);//create right node
					count++;
					ybest = max(ybest, mn->right->y.first);
					X->conservativeResize(X->rows(), X->cols() + 1);
					Y->conservativeResize(Y->rows() + 1);
					X->col(X->cols() - 1) = mn->right->x;
					(*Y)(Y->rows() - 1) = mn->right->y.first;



				}
				else
					mn->insert("right", mn->x + c.cwiseProduct(mn->d), make_pair(GPm.mu - B*sqrt(abs(GPm.sig2)),-1), mn->d, false);//create right node
				vmax = mn->y.first;
			}
		}
		cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
	}
	delete X, Y;
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	cout << "elapsed time: " << elpased.count() << "seconds" << endl;
	return ybest;
}
