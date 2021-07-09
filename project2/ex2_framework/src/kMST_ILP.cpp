#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance &_instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), epsInt( 0.0 ), epsOpt( 0.0 )
{
	if( k == 0 ) k = instance.n_nodes;
}

void kMST_ILP::write_output()
{	
	ofstream myfile;
  	myfile.open (model_type+".txt", ios::out | ios::app);
  	myfile << model_type<<"\n";
  	myfile << k<<"\n";
  	myfile << instance.n_nodes<<"\n";
  	myfile << cplex.getStatus()<<"\n";
  	myfile << cplex.getNnodes()<<"\n";
  	myfile << cplex.getObjValue()<<"\n";
  	myfile << Tools::CPUtime()<<"\n";
  	myfile << cplex.getMIPRelativeGap()<<"\n";
  	myfile << cplex.getIncumbentNode()<<"\n";
  	myfile.close();
	
}

void kMST_ILP::solve()
{
	try {

		// initialize CPLEX solver
		initCPLEX();

		// add common constraints
		//modelCommon();
		// add model-specific constraints
		if( model_type == "scf" ) modelSCF();
		else if( model_type == "mcf" ) modelMCF();
		else if( model_type == "mtz" ) modelMTZ();
		else if( model_type == "dcc" ) modelCommon();//cout << "DC-CUT model: no additional constraints\n";
		else if( model_type == "cec" ) modelCommon();//cout << "CE-CUT model: no additional constraints\n";
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
		cplex.setParam( IloCplex::Param::WorkMem, 6144 ); // set memory limit to 8 GB 8192

		epsInt = cplex.getParam( IloCplex::Param::MIP::Tolerances::Integrality );
		epsOpt = cplex.getParam( IloCplex::Param::Simplex::Tolerances::Optimality );

		cout<<"integer tolerance "<<epsInt<<"\n";
		cout<<"lp tolerance "<<epsOpt<<"\n";
		

		// set cut-callback for cycle-elimination cuts ("cec") or directed connection cuts ("dcc")
		// both for integer (mandatory!) and fractional (optional) solutions
		CutCallback cb( model_type,epsInt, epsOpt, instance, x, z );
		if( model_type == "dcc" || model_type == "cec" ) {
			cout<<"setting cut callback \n";
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation;
			cplex.use( &cb, contextmask );
		}

		// solve model cplex.getMIPRelativeGap();
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished.\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";
		
		//IloNumArray vals_x(env, instance.n_edges*2);
		//cplex.getValues(vals_x, x);
		
		//for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
		//	cout<<vals_x[i]<<", ";
		//}
		//cout<<"\n";
		//for( u_int i = 1; i <  instance.n_edges*2; i+=2 ) {
		//	cout<<vals_x[i]<<", ";
		//}
		//cout<<"\n";
		//for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
		//	if(vals_x[i]>epsInt)
		//		cout<<"("<<instance.edges[i/2].v1<<","<<instance.edges[i/2].v2<<","<<vals_x[i]<<") ";
		//	if(vals_x[i+1]>epsInt)
		//		cout<<"("<<instance.edges[i/2].v2<<","<<instance.edges[i/2].v1<<","<<vals_x[i+1]<<") ";
				
		//}
		//cout<<"\n";
		
		//IloNumArray vals_z(env, instance.n_nodes);
		//cplex.getValues(vals_z, z);
		//for( u_int i = 0; i <  instance.n_nodes; i++ ) {
		//if(vals_z[i]>epsInt)
		//	cout<<i<<", ";
		//}
		


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
		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes);
		
		IloExpr objective( env );
		IloExpr sum_x( env );
		IloExpr sum_z( env );

		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				objective += (x[i] + x[i+1])*instance.edges[i/2].weight;
		}
		model.add(IloMinimize(env, objective));

		// root has exactly one outgoing edge
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if(instance.edges[i/2].v1==0)
				sum_x += x[i];
			if(instance.edges[i/2].v2==0)
				sum_x += x[i+1];
		}
		model.add(sum_x == 1);
		sum_x.clear();
		
		// no bidirectional arcs
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			model.add(x[i] + x[i+1] <= 1);
		}

		
		// no edge back into root
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if(instance.edges[i/2].v2==0)
				sum_x += x[i];
			if(instance.edges[i/2].v1==0)
				sum_x += x[i+1];
		}
		model.add(sum_x == 0);
		sum_x.clear();
		
		// dont select root
		model.add(z[0] == 0);
		sum_z.clear();
		
		// select exactly k nodes
		for( u_int i = 1; i <  instance.n_nodes; i++ ) {
			sum_z += z[i];
		}
		model.add(sum_z == k);
		sum_z.clear();
		
		// exactly one ingoing edge at any node
		for( u_int i = 1; i <  instance.n_nodes; i++ ) {
			for( u_int j = 0 ; j < instance.n_edges*2; j+=2 ){
				if(instance.edges[j/2].v2 == i)
					sum_x += x[j];
				if(instance.edges[j/2].v1 == i)
					sum_x += x[j+1];
			}
			model.add(sum_x == z[i]);
			sum_x.clear();
		}
		
		// exactly one ingoing edge at any node
		for( u_int i = 1; i <  instance.n_nodes; i++ ) {
			for( u_int j = 0 ; j < instance.n_edges*2; j+=2 ){
				if(instance.edges[j/2].v1 == i)
					model.add(x[j] <= z[i]);
				if(instance.edges[j/2].v2 == i)
					model.add(x[j+1] <= z[i]);
			}
		}
		
		// k-1 arcs in total, since tree
		for( u_int i = 0; i <  instance.n_edges*2; i++ ) {
			sum_x += x[i];
		}
		model.add(sum_x == k);
		sum_x.clear();


		objective.end();
		sum_x.end();
		sum_z.end();
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
		
		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes);
		f = IloIntVarArray(env, instance.n_edges*2, 0, k);
		
		IloExpr objective( env );
		IloExpr sum_x( env );
		IloExpr sum_x1( env );
		IloExpr sum_f( env );
		IloExpr sum_f1( env );
		IloExpr sum_z( env );
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				objective += (x[i] + x[i+1])*instance.edges[i/2].weight;
		}
		model.add(IloMinimize(env, objective));
		
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if( instance.edges[i/2].v1==0 ){
				sum_x  += x[i]; 
				sum_x1 += x[i+1]; 
				sum_f  += f[i];
				sum_f1 += f[i+1];
			}
			if( instance.edges[i/2].v2==0 ){
				sum_x  += x[i+1]; 
				sum_x1 += x[i]; 
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
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			model.add(f[i] <= k*x[i]);
			model.add(f[i+1] <= k*x[i+1]);
		}
		
		for( u_int i = 0; i <  instance.n_nodes; i++ ) {
			sum_z += z[i];
		}
		model.add(sum_z == k+1);
		sum_z.clear();
		
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
		sum_f1.end();
		sum_x1.end();
		sum_x.end();
		sum_z.end();
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

		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes-1);
		f = IloIntVarArray(env, instance.n_edges*2*(instance.n_nodes-1), 0, 1);
		
		IloExpr objective( env );
		IloExpr sum_x( env );
		IloExpr sum_x1( env );
		IloExpr sum_f( env );
		IloExpr sum_f1( env );
		IloExpr sum_z( env );
		
		u_int stride = 2*instance.n_edges;
		
		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
				objective += (x[i] + x[i+1])*instance.edges[i/2].weight;
		}
		model.add(IloMinimize(env, objective));

		for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			if( instance.edges[i/2].v1==0 ){
				sum_x  += x[i]; 
				sum_x1 += x[i+1]; 
			}
			if( instance.edges[i/2].v2==0 ){
				sum_x  += x[i+1]; 
				sum_x1 += x[i]; 
			}
		}
		model.add(sum_x == 1);
		model.add(sum_x1 == 0);
		sum_x.clear();
		sum_x1.clear();
		
		for( u_int i = 0; i < instance.n_nodes-1; i++ ) {
			for( u_int j = 0; j<2*instance.n_edges; j+=2 ) {
				if( instance.edges[j/2].v1==0 ){
					sum_f  += f[j +i*stride]; // out of root x0j
					sum_f1 += f[j+1 +i*stride]; // to root xj0
				}
				if( instance.edges[j/2].v2==0 ){
					sum_f  += f[j+1 +i*stride]; // out of root x0j
					sum_f1 += f[j +i*stride]; // to root xj0
				}
			}
			model.add(sum_f == z[i]);
			model.add(sum_f1 == 0);
			sum_f.clear();
			sum_f1.clear();
		}
		
		for( u_int i = 0; i < instance.n_nodes-1; i++ ) {
			for( u_int j = 0; j<2*instance.n_edges; j+=2 ) {
				model.add( f[j +i*stride] <= x[j]);
				model.add( f[j+1 +i*stride] <= x[j+1]);
			}
		}
		
		for( u_int i = 0; i <  instance.n_nodes-1; i++ ) {
			sum_z += z[i];
		}
		model.add(sum_z == k);
		sum_z.clear();
		
		for( u_int n = 0; n < instance.n_nodes-1; n++ ) {
			for( u_int j = 1; j < instance.n_nodes; j++ ){
				for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
					if(instance.edges[i/2].v1 ==j)
						sum_x += -f[i +n*stride] + f[i+1 +n*stride];
					if(instance.edges[i/2].v2 ==j)
						sum_x +=  f[i +n*stride] - f[i+1 +n*stride];
				}
				if(j==(n+1)){
					model.add(sum_x == z[n]);}
				else{
					model.add(sum_x == 0);}
				sum_x.clear();
			}
		}

		objective.end();
		sum_x.end();
    sum_x1.end();
		sum_f.end();
		sum_f1.end();
		sum_z.end();

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

		x = IloBoolVarArray(env, instance.n_edges*2);
		z = IloBoolVarArray(env, instance.n_nodes);
		f = IloIntVarArray(env);
		f.add(IloIntVar(env, 0, 0)); 
		f.add(IloIntVarArray(env, instance.n_nodes-1, 1, k));
		
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
		
		// exactly one ingoing edge at any node
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
	if( model_type == "scf" || model_type == "mcf" || model_type == "mtz" ) f.end();
	x.end();
	z.end();
//	if( model_type == "scf" || model_type == "mcf" ) f.end();
//	else if( model_type == "mtz" ) d.end();
	cplex.end();
	model.end();
	env.end();
}
