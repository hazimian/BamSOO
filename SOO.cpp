#include "node.h"
#include "SOO.h"
#include "Common.h"

struct GP{
	double mu;
	double sig2;
};
double kernel(VectorXd x1, VectorXd x2, double p){
	VectorXd V = x1 - x2;
	return pow(10, p)*exp(-V.transpose()*V);
}

GP FitGP(const MatrixXd* x, const VectorXd* y, const VectorXd xs){
	MatrixXd Kmat = MatrixXd::Identity(x->outerSize(), x->outerSize());
	MatrixXd iKmat = MatrixXd::Identity(x->outerSize(), x->outerSize());
	MatrixXd iKmat0 = iKmat;
	VectorXd kvec = VectorXd::Zero(x->outerSize());
	VectorXd kvec0 = kvec;
	vector<double> pvec = {-1};
	double p0 = 0;
	double LogL = nINF;
	for (int m = 0; m < pvec.size(); m++){
		for (int i = 0; i < x->outerSize(); i++){
			kvec(i) = kernel(x->col(i), xs, pvec[m]);//creates kvec
			for (int j = 0; j < x->outerSize(); j++){
				Kmat(i, j) = kernel(x->col(i), x->col(j), pvec[m]);//creat Kmat
			}
		}
		iKmat = (Kmat + 1e-8*MatrixXd::Identity(Kmat.innerSize(), Kmat.innerSize())).inverse();
		if (y->transpose()*iKmat*(*y) + log(Kmat.determinant())>LogL){
			p0 = pvec[m];
			iKmat0 = iKmat;
			kvec0 = kvec;
		}

	}


	double mu = kvec.transpose()*iKmat*(*y);
	double sig2 = kernel(xs, xs, p0) - kvec.transpose()*iKmat*kvec;
	return{ mu, sig2 };
}

VectorXd shift(VectorXd &d){
	VectorXd z = VectorXd::Zero(d.size());
	VectorXd::Index k = 0;
	d.maxCoeff(&k);
	z[k] = 1;
	d[k] = d[k] / 3.0;
	return z;

}

void fopt(node* root, node** mn, int h, double &m){

	if (!root || root->h > h)
		return;
	if (root->h == h){//consider only if has no children
		if (!root->left && !root->middle & !root->right){
			if (root->y > m){
				m = root->y;
				*mn = root;
			}
			return;
		}
		else
			return;
	}
	//otherwise
	if (root->left)
		fopt(root->left, mn, h, m);
	if (root->middle)
		fopt(root->middle, mn, h, m);
	if (root->right)
		fopt(root->right, mn, h, m);



}

VectorXd denorm(VectorXd x, VectorXd lb, VectorXd ub){
	return .5*(ub - lb).cwiseProduct(x) + .5*(lb + ub);
}


double SOO(function<double(VectorXd, void* dataPtr)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax, void *dataPtr, VectorXd &xbest, const std::shared_ptr<BO4Input> mInp){
	auto start = chrono::high_resolution_clock::now();
	double ybest = nINF;
  
	//create the root node
	//VectorXd mid = 0.5*(lb + ub);
	//VectorXd delta = (ub - lb);
	VectorXd mid = VectorXd::Zero(lb.rows());//normalized to [-1,1]
	VectorXd delta = 2 * VectorXd::Ones(lb.rows());
	node* root = new node(mid, fun(denorm(mid, lb, ub), dataPtr), delta);
	
  if (root->y > ybest)
  {
    ybest = root->y;
    xbest = root->x;
  }

	//now create level h=1

	VectorXd c = shift(root->d);// d passed by reference
	//now insert into h=1
  root->insert(NodeType::LEFT, root->x - c.cwiseProduct(root->d), fun(denorm(root->x - c.cwiseProduct(root->d), lb, ub), dataPtr), root->d);
	
  if (root->left->y > ybest)
  {
    ybest = root->left->y;
    xbest = root->left->x;
  }

  root->insert(NodeType::MIDDLE, root->x, root->y, root->d);
  root->insert(NodeType::RIGHT, root->x + c.cwiseProduct(root->d), fun(denorm(root->x + c.cwiseProduct(root->d), lb, ub), dataPtr), root->d);
	
  if (root->right->y > ybest)
  {
    ybest = root->right->y;
    xbest = root->right->x;
  }
  
  TextInfo txtInfo;
  stringstream strStream;
	int count = 3;
	double m;
	node* mn = nullptr;
	int N = 3;
	int depth = 1;
	while (count < pmax){
		depth = node::get_depth();//update depth
		for (int h = 1; h <= depth; ++h){
			m = nINF;//initialize 
			mn = nullptr;
			fopt(root, &mn, h, m);//find the optimum to expand
			if (mn){
				//now generate children
				c = shift(mn->d);
				N += 1;
        mn->insert(NodeType::MIDDLE, mn->x, mn->y, mn->d);//just copy parent into middle
        mn->insert(NodeType::LEFT, mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d), lb, ub), dataPtr), mn->d);//create left node
				count++;
				
        if (mn->left->y > ybest)
        {
          ybest = mn->left->y;
          xbest = mn->left->x;
        }
        //cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
        /*strStream.str("");
        strStream.clear();
        strStream << count << ",,,,," << (abs(ybest - yopt)) << "\r\n";
        txtInfo.txt = strStream.str();
        txtInfo.txtType = "CONSOLE";
        mInp->callSetText(txtInfo);*/

        mn->insert(NodeType::RIGHT, mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d), lb, ub), dataPtr), mn->d);//create right node
				count++;
        if (mn->right->y > ybest)
        {
          ybest = mn->right->y;
          xbest = mn->right->x;
        }
        //cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
        /*strStream.str("");
        strStream.clear();
        strStream << count << ",,,,," << (abs(ybest - yopt)) << "\r\n";
        txtInfo.txt = strStream.str();
        txtInfo.txtType = "CONSOLE";
        mInp->callSetText(txtInfo);*/
			}
		}
    strStream.str("");
    strStream.clear();
    strStream << count << ",,,,," << (abs(ybest - yopt)) << "\r\n";
    txtInfo.txt = strStream.str();
    txtInfo.txtType = "CONSOLE";
    mInp->callSetText(txtInfo);
	}
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	//cout <<"elapsed time: "<< elpased.count() << "seconds" << endl;
  xbest = denorm(xbest, lb, ub);
	return ybest;
}






double BamSOO(function<double(VectorXd)> fun, VectorXd lb, VectorXd ub, double yopt, int pmax){
	auto start = chrono::high_resolution_clock::now();
	double ybest = nINF;
	//create the root node
	//VectorXd mid = 0.5*(lb + ub);
	//VectorXd delta = (ub - lb);
	VectorXd mid = VectorXd::Zero(lb.rows());//normalized to [-1,1]
	VectorXd delta = 2*VectorXd::Ones(lb.rows());
	node* root = new node(mid, fun(mid), delta);
	ybest = min(ybest, root->y);


	//now create level h=1
	VectorXd c = shift(root->d);// d passed by reference
	//now insert into h=1
	root->insert(NodeType::LEFT, root->x - c.cwiseProduct(root->d), fun(denorm(root->x - c.cwiseProduct(root->d), lb, ub)), root->d);
	ybest = min(ybest, root->left->y);
  root->insert(NodeType::MIDDLE, root->x, root->y, root->d);
  root->insert(NodeType::RIGHT, root->x + c.cwiseProduct(root->d), fun(denorm(root->x + c.cwiseProduct(root->d), lb, ub)), root->d);
	ybest = min(ybest, root->right->y);


	int count = 3;
	double m;
	node* mn = nullptr;
	int N = 1;
	int depth = 1;
	double B;
	GP GPm;
	MatrixXd *X=new MatrixXd(lb.size(),3);
	VectorXd *Y=new VectorXd(3);
	(*X) << root->x, root->left->x, root->right->x;
	(*Y) << root->y, root->left->y, root->right->y;
	MatrixXd temp;
	double UCB = 0;
	while (count < pmax){
		depth = node::get_depth();//update depth
		for (int h = 1; h <= depth; ++h){
			m = nINF;//initialize 
			mn = nullptr;
			fopt(root, &mn, h, m);//find the optimum to expand
			if (mn){
				//now generate children
				c = shift(mn->d);

        mn->insert(NodeType::MIDDLE, mn->x, mn->y, mn->d);//just copy parent into middle
				B = sqrt(2 * log2(pow(PI*N, 2) / (6 * .05)));
				//do for left element
				GPm = FitGP(X, Y, mn->x - c.cwiseProduct(mn->d));//fit GP
				if (GPm.mu + B*sqrt(abs(GPm.sig2)) > ybest){//if greater than UCB evaluate
          mn->insert(NodeType::LEFT, mn->x - c.cwiseProduct(mn->d), fun(denorm(mn->x - c.cwiseProduct(mn->d), lb, ub)), mn->d);//create left node
					count++;
					ybest = max(ybest, mn->left->y);

					X->conservativeResize(X->rows(), X->cols() + 1);
					Y->conservativeResize(Y->rows() + 1);
					X->col(X->cols()-1)=mn->left->x;
					(*Y)(Y->rows() - 1) = mn->left->y;
					//cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
					N += 1;
				}
				else//otherwise use LCB
          mn->insert(NodeType::LEFT, mn->x - c.cwiseProduct(mn->d), GPm.mu - B*sqrt(abs(GPm.sig2)), mn->d);//create left node
				//repeat for right element
				GPm = FitGP(X, Y, mn->x + c.cwiseProduct(mn->d));//fit GP
				if (GPm.mu + B*sqrt(abs(GPm.sig2)) > ybest){//if greater than UCB evaluate
          mn->insert(NodeType::RIGHT, mn->x + c.cwiseProduct(mn->d), fun(denorm(mn->x + c.cwiseProduct(mn->d), lb, ub)), mn->d);//create right node
					count++;
					ybest = max(ybest, mn->right->y);
					X->conservativeResize(X->rows(), X->cols() + 1);
					Y->conservativeResize(Y->rows() + 1);
					X->col(X->cols() - 1) = mn->right->x;
					(*Y)(Y->rows() - 1) = mn->right->y;
					//cout << count << ",,,,," << log10(abs(ybest - yopt)) << endl;
					N += 1;
				}
				else
          mn->insert(NodeType::RIGHT, mn->x + c.cwiseProduct(mn->d), GPm.mu - B*sqrt(abs(GPm.sig2)), mn->d);//create right node
			}
		}
	}
	delete X, Y;
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elpased = finish - start;
	//cout << "elapsed time: " << elpased.count() << "seconds" << endl;
	return ybest;
}