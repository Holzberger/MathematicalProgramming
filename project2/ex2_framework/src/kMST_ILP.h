#ifndef __KMST_ILP__H__
#define __KMST_ILP__H__

#include "Tools.h"
#include "Instance.h"
#include "CutCallback.h"
#include <ilcplex/ilocplex.h>

#include <stdlib.h>

using namespace std;

ILOSTLBEGIN

class kMST_ILP
{

private:

	Instance &instance;
	string model_type;
	int k;

	IloEnv env;
	IloModel model;
	IloCplex cplex;

	IloBoolVarArray x; // edge or arc selection variables
	IloBoolVarArray z; // node selection variables
	IloIntVarArray f; // flow variables
	IloIntVarArray d; // node distance variables

	IloNumArray values; // to store variable values

	double epsInt, epsOpt;

	void initCPLEX();
	void modelCommon();
	void modelSCF();
	void modelMCF();
	void modelMTZ();
	
	// count cuts
	u_int u_cuts; // user cuts
	u_int c_cuts; // candidate cuts
	

public:

	kMST_ILP( Instance &_instance, string _model_type, int _k );
	~kMST_ILP();
	void solve();
	void write_output();

};

#endif //__KMST_ILP__H__
