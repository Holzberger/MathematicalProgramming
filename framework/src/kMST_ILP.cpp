#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance &_instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), epsInt( 0.0 ), epsOpt( 0.0 )
{
	if( k == 0 ) k = instance.n_nodes;
}

void kMST_ILP::solve()
{
	try {

		// initialize CPLEX solver
		initCPLEX();

		// add common constraints
		modelCommon();
		// add model-specific constraints
		if( model_type == "scf" ) modelSCF();
		else if( model_type == "mcf" ) modelMCF();
		else if( model_type == "mtz" ) modelMTZ();
		else if( model_type == "dcc" ) cout << "DC-CUT model: no additional constraints\n";
		else if( model_type == "cec" ) cout << "CE-CUT model: no additional constraints\n";
		else {
			cerr << "No existing model chosen\n";
			exit( -1 );
		}

		// build model
		cplex = IloCplex( model );
		cplex.exportModel( "model.lp" );

		// set parameters
		cplex.setParam( IloCplex::Param::Threads, 1 ); // only use a single thread
		cplex.setParam( IloCplex::Param::TimeLimit, 3600 ); // set time limit to 1 hour
		cplex.setParam( IloCplex::Param::WorkMem, 8192 ); // set memory limit to 8 GB

		epsInt = cplex.getParam( IloCplex::Param::MIP::Tolerances::Integrality );
		epsOpt = cplex.getParam( IloCplex::Param::Simplex::Tolerances::Optimality );

		// set cut-callback for cycle-elimination cuts ("cec") or directed connection cuts ("dcc")
		// both for integer (mandatory!) and fractional (optional) solutions
		CutCallback cb( model_type, epsOpt, instance, x, z );
		if( model_type == "dcc" || model_type == "cec" ) {
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation;
			cplex.use( &cb, contextmask );
		}

		// solve model
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished.\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";
		IloNumArray vals(env);
		IloNumArray valsf(env);
		cplex.getValues(vals, x);
		cplex.getValues(valsf, f);
		
		for(u_int i=0; i<instance.n_edges*2; i+=2){
			if (abs(vals[i]) > 1e-4)
				env.out() << i/2 <<"["<<valsf[i]<< "("<<instance.edges[i/2].v1 <<","<< instance.edges[i/2].v2<<"),";
			if (abs(vals[i+1]) > 1e-4)
				env.out() << i/2 <<"["<<valsf[i+1]<< "("<< instance.edges[i/2].v2 <<","<< instance.edges[i/2].v1<<"),";
		}
		env.out() << "\n";
		cplex.getValues(vals, f);
		env.out() << "Values = " << vals << endl;

	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
		exit( -1 );
	}
}

// ----- private methods -----------------------------------------------

void kMST_ILP::initCPLEX()
{
	cout << "initialize CPLEX ... ";
	try {
		env = IloEnv();
		model = IloModel( env );
	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
	}
	cout << "done.\n";
}

void kMST_ILP::modelCommon()
{
	try {

		// +++++++++++++++++++++++++++++++++++++++++++++++
		// TODO create variables, build common constraints
		// +++++++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelCommon: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelCommon: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelSCF()
{
	try {
		
		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build single commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++
		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes);
		f = IloIntVarArray(env, instance.n_edges*2, 0, k);
		
		IloExpr objective( env );
		IloExpr sum_x( env );
		IloExpr sum_x1( env );
		IloExpr sum_f( env );
		IloExpr sum_f1( env );
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				objective += (x[i] + x[i+1])*instance.edges[i/2].weight;
		}
		model.add(IloMinimize(env, objective));
		
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if( instance.edges[i/2].v1==0 ){
				sum_x  += x[i]; // out of root x0j
				sum_x1 += x[i+1]; // to root xj0
				sum_f  += f[i];
				sum_f1 += f[i+1];
			}
			if( instance.edges[i/2].v2==0 ){
				sum_x  += x[i+1]; // out of root x0j
				sum_x1 += x[i]; // to root xj0
				sum_f  += f[i+1]; 
				sum_f1 += f[i];
			}
		}
		model.add(sum_x == 1);
		model.add(sum_x1 == 0);
		model.add(sum_f == k);
		model.add(sum_f1 == 0);
		sum_x.clear();
		sum_x1.clear();
		sum_f.clear();
		sum_f1.clear();
		
		// select arc whenever there is a flow
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			model.add(f[i] <= k*x[i]);
			model.add(f[i+1] <= k*x[i+1]);
		}
		
		IloExpr sum_z( env );
		
		// select exactly k+1 nodes
		for( u_int i = 0; i <  instance.n_nodes; i++ ) {
			sum_z += z[i];
		}
		model.add(sum_z == k+1);
		sum_z.clear();
		
		// flow continuity
		for( u_int j = 1; j < instance.n_nodes; j++ ){
			for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				if(instance.edges[i/2].v1 ==j)
					sum_x += -f[i] + f[i+1];
				if(instance.edges[i/2].v2 ==j)
					sum_x +=  f[i] - f[i+1];
			}
			model.add(sum_x == z[j]);
			sum_x.clear();
		}

		
		sum_f.end();
		sum_x1.end();
		sum_x.end();
		objective.end();
	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelSCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelSCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMCF()
{
	try {

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build multi commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMTZ()
{
	try {
		for( u_int i = 0; i <  instance.n_edges; i++ ) {
		cout << instance.edges[i].v1 << " ";
		cout << instance.edges[i].v2 << " ";
		cout << instance.edges[i].weight << "\n";
	}

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build Miller-Tucker-Zemlin model
		// ++++++++++++++++++++++++++++++++++++++++++
		u_int root_id = 0;
		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes);
		f = IloIntVarArray(env);
		f.add(IloIntVar(env, 0, 0)); // root has 0 value
		f.add(IloIntVarArray(env, instance.n_nodes-1, 1, k-1));
		
		
		IloExpr objective( env );
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				objective += (x[i] + x[i+1])*instance.edges[i/2].weight;
		}
		model.add(IloMinimize(env, objective));
		
		IloExpr sum_x( env );
		IloExpr sum_z( env );
		
		// root has exactly one outgoing edge
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if(instance.edges[i/2].v1==0)
				sum_x += x[i];
			if(instance.edges[i/2].v2==0)
				sum_x += x[i+1];
		}
		model.add(sum_x == 1);
		sum_x.clear();
		
		// select exactly k nodes
		for( u_int i = 0; i <  instance.n_nodes; i++ ) {
			sum_z += z[i];
		}
		model.add(sum_z == k);
		sum_z.clear();
		
		// MTZ constraint
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			model.add( f[instance.edges[i/2].v1] + x[i] <= f[instance.edges[i/2].v2] + k*(1-x[i]) );
			//if(instance.edges[i/2].v1!=0)
				model.add(f[instance.edges[i/2].v2] + x[i+1] <= f[instance.edges[i/2].v1] + k*(1-x[i+1]) );
		}
		
		// exactly one outgoing edge at any node
		for( u_int i = 0; i <  instance.n_nodes; i++ ) {
			for( u_int j = 0 ; j < instance.n_edges*2; j+=2 ){
				if(instance.edges[j/2].v2 == i)
					sum_x += x[j];
				if(instance.edges[j/2].v1 == i)
					sum_x += x[j+1];
			}
			model.add(sum_x == z[i]);
			sum_x.clear();
		}
		
		// there is only an outgoing edge possible if outgoing node is selected
		for( u_int i = 1; i <  instance.n_nodes; i++ ) {
			for( u_int j = 0 ; j < instance.n_edges*2; j+=2 ){
				if(instance.edges[j/2].v1 == i)
					sum_x += x[j];
				if(instance.edges[j/2].v2 == i)
					sum_x += x[j+1];
			}
			model.add(sum_x <= (k-1)*z[i]);
			sum_x.clear();
		}
		


		sum_z.end();
		sum_x.end();
		objective.end();
	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMTZ: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMTZ: unknown exception.\n";
		exit( -1 );
	}
}

kMST_ILP::~kMST_ILP()
{
	// free CPLEX resources
//	x.end();
//	z.end();
//	if( model_type == "scf" || model_type == "mcf" ) f.end();
//	else if( model_type == "mtz" ) d.end();
	cplex.end();
	model.end();
	env.end();
}
