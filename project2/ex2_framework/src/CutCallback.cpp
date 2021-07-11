#include "CutCallback.h"

CutCallback::CutCallback( string _cut_type, double _epsI, double _epsO, 
					      Instance &_instance, IloBoolVarArray &_x, IloBoolVarArray &_z,
						  u_int &_u_cuts, u_int &_c_cuts) :
	cut_type( _cut_type ), epsI( _epsI ), epsO( _epsO), instance( _instance ), context( NULL ), x( _x ), z( _z ), u_cuts(_u_cuts), c_cuts(_c_cuts)
{
	arc_weights.resize( 2 * instance.n_edges );
	

	for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
			arcs.push_back({instance.edges[i/2].v1, instance.edges[i/2].v2});
			arcs.push_back({instance.edges[i/2].v2, instance.edges[i/2].v1});
	}
	
	capacities = vector<double>(instance.n_edges*2, 0);
	
	n_ucuts = 0;
	multiplier_ucuts = 1;
	
	u_cuts = 0;
	c_cuts = 0;
	
	

}

CutCallback::~CutCallback()
{
}

void CutCallback::invoke( const IloCplex::Callback::Context &_context )
{
	context = &_context;
	env = context->getEnv();
	xsol = IloNumArray( env, 2 * instance.n_edges );
	zsol = IloNumArray( env, instance.n_nodes );

	// integer solution -> mandatory check for feasibility
	if( context->inCandidate() ) {
		if( !context->isCandidatePoint() ) throw IloCplex::Exception( -1, "Unbounded solution" );
		context->getCandidatePoint( x, xsol );
		context->getCandidatePoint( z, zsol );
		if( cut_type == "dcc" ) connectionCuts();
		else if( cut_type == "cec" ) cycleEliminationCuts();
	}
	// fractional solution -> optional to improve dual bounds
	else if( context->inRelaxation() ) {
		context->getRelaxationPoint( x, xsol );
		context->getRelaxationPoint( z, zsol );
		if( cut_type == "dcc" ) connectionCuts();
		else if( cut_type == "cec" ) cycleEliminationCuts();
	}
	else {
		xsol.end();
		zsol.end();
		throw IloCplex::Exception( -1, "Unexpected contextID" );
	}

	xsol.end();
	zsol.end();
}

/*
 * separation of directed connection cut inequalities
 */
void CutCallback::connectionCuts()
{
	try {
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// TODO find violated directed connection cut inequalities
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		IloExpr sum_x( env );
		
		// update capacities
		for( u_int i = 0; i <  instance.n_edges*2; i++ )
			capacities[i] = xsol[i];
		
		Maxflow mc(instance.n_nodes, instance.n_edges*2, arcs );
		// current target node for maxflow-algorithm
		u_int t = 1;
		mc.update( 0, t, &capacities[0] );
		// save cutsets here
		vector<int> my_cut(instance.n_nodes, 0);
		double cap_mincut = 0;
		// arcs in the cut
		u_int arc_counter = 0;
		
		while(t<instance.n_nodes){
			if(zsol[t]> (1-epsI) ){ // check if target node is in solution
				cap_mincut = -1;
				cap_mincut = mc.min_cut( 10, &my_cut[0] );
				for( u_int i = 0; i <  instance.n_edges*2; i+=2 ) {
					if(my_cut[instance.edges[i/2].v1] < 2 && my_cut[instance.edges[i/2].v2] == 2){
						sum_x += x[i];
						arc_counter++;
					}
					if(my_cut[instance.edges[i/2].v2] < 2 && my_cut[instance.edges[i/2].v1] == 2){
						sum_x += x[i+1];
						arc_counter++;
					}
				}
				if( (cap_mincut  - zsol[t]) < -((double)arc_counter+1)*epsO){ // true if violated dcc is found
					break;
				}else{ // no violated dcc
					sum_x.clear();
				}
			}
			arc_counter = 0;
			// next target
			t++;
			mc.update( 0, t);
		}
			// add found violated cut to model
			switch( context->getId() ) {
				case IloCplex::Callback::Context::Id::Candidate:
					if(arc_counter > 0){ // if there is a violated dcc with more than 0 arcs  
						context->rejectCandidate( IloRange( env,  -((double)arc_counter+1)*epsO, sum_x - z[t], IloInfinity  ) );
						// decrease max number of user-cuts till next integer solution
						multiplier_ucuts *= 0.75;
						n_ucuts = 0;
						
						c_cuts++;
					}
					break;
				case IloCplex::Callback::Context::Id::Relaxation:
					if( (arc_counter > 0) && (n_ucuts < (double)instance.n_nodes*multiplier_ucuts) ){ // if there is a violated dcc with more than 0 arcs  
						context->addUserCut( IloRange( env,  -((double)arc_counter+1)*epsO, sum_x - z[t], IloInfinity  ), 
											 IloCplex::UseCutForce, IloFalse  );
						n_ucuts++;
						
						u_cuts++;
					}
					break;
				default:
					throw IloCplex::Exception( -1, "Unexpected contextID" );
			}
			sum_x.clear();
			sum_x.end();
 
	}
	catch( IloException &e ) {
		cerr << "CutCallback: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "CutCallback: unknown exception.\n";
		exit( -1 );
	}
	
	
}

/*
 * separation of cycle elimination cut inequalities
 */
void CutCallback::cycleEliminationCuts()
{
	try {
		

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// TODO find violated cycle elimination cut inequalities
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		
		
		
		IloExpr sum_x( env );
		
		

		switch( context->getId() ) {
			case IloCplex::Callback::Context::Id::Candidate:
				for( u_int j = 0 ; j < instance.n_edges*2; j+=2 ){
					arc_weights[j/2] = 1 - xsol[j];
					arc_weights[instance.n_edges+j/2] = 1 - xsol[j+1];
				}
				violated_path = 0;
				for( u_int j = 0; j < instance.n_edges; j++ ){
					sp = shortestPath(instance.edges[j].v2, instance.edges[j].v1 );
					if( ( sp.weight+arc_weights[j]) < (1-epsO*((double)sp.path.size()+1)) ){
						violated_path = 1;
						// sum_x += x[j*2];
						sum_x += x[j*2]+x[j*2+1];
						break;
					}
					sp = shortestPath(instance.edges[j].v1, instance.edges[j].v2 );
					if( ( sp.weight+arc_weights[j+instance.n_edges]) < (1-epsO*((double)sp.path.size()+1)) ){
						violated_path = 1;
						// sum_x += x[j*2+1];
						sum_x += x[j*2]+x[j*2+1];
						break;
					}
				}
				
				
				if( violated_path==1 ){
					for(auto it = sp.path.begin(); it !=sp.path.end(); it++)
						if(*it<instance.n_edges){
							sum_x += x[(*it)*2];
							sum_x += x[(*it)*2+1];//
						}else{	
							sum_x += x[(*it-instance.n_edges)*2];//
							sum_x += x[(*it-instance.n_edges)*2+1];
						}
				}
				
				if(violated_path==1){
				context->rejectCandidate( IloRange( env, -IloInfinity , sum_x, (double)sp.path.size()  ) );
				c_cuts++;
				}
				sum_x.clear();
				break;
			case IloCplex::Callback::Context::Id::Relaxation:
				break;
			default:
				throw IloCplex::Exception( -1, "Unexpected contextID" );
		}
		sum_x.end();

	}
	catch( IloException &e ) {
		cerr << "CutCallback: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "CutCallback: unknown exception.\n";
		exit( -1 );
	}
	
}

/*
 * Dijkstra's algorithm to find a shortest path
 * Note: slow implementation with vectors in time O(n^2), instead of heaps etc.
 */
CutCallback::SPResultT CutCallback::shortestPath( u_int source, u_int target )
{
	vector<SPNodeT> nodes( instance.n_nodes );
	vector<bool> finished( instance.n_nodes, false ); // indicates finished nodes

	// initialization
	for( u_int v = 0; v < instance.n_nodes; v++ ) {
		nodes[v].pred = -1;
		nodes[v].pred_arc_id = -1;
		if( v == source ) nodes[v].weight = 0;
		else nodes[v].weight = numeric_limits<double>::max();
	}

	while( true ) {

		// find unfinished node with minimum weight to examine next
		// (should usually be done with heap or similar data structures)
		double wmin = numeric_limits<double>::max();
		u_int v = 0;
		for( u_int u = 0; u < instance.n_nodes; u++ ) {
			if( !finished[u] && nodes[u].weight < wmin ) {
				wmin = nodes[u].weight;
				v = u;
			}
		}

		// if all reachable nodes are finished or target node is reached -> stop
		if( wmin == numeric_limits<double>::max() || v == target ) break;

		// this node is finished now
		finished[v] = true;

		// update all adjacent nodes on outgoing arcs
		for( u_int e : instance.incidentEdges[v] ) {
			u_int a; // according arc id
			u_int u; // adjacent node
			if( instance.edges[e].v1 == v ) {
				a = e;
				u = instance.edges[e].v2;
			}
			else {
				a = e + instance.n_edges;
				u = instance.edges[e].v1;
			}
			// only examine adjacent node if unfinished
			if( !finished[u] ) {
				// check if weight at node u can be decreased
				if( nodes[u].weight > nodes[v].weight + arc_weights[a] ) {
					nodes[u].weight = nodes[v].weight + arc_weights[a];
					nodes[u].pred = v;
					nodes[u].pred_arc_id = a;
				}
			}
		}
	}

	SPResultT sp;
	sp.weight = 0;
	int v = target;
	while( v != (int) source && v != -1 ) {
		int a = nodes[v].pred_arc_id;
		if( a < 0 ) break;
		sp.weight += arc_weights[a];
		sp.path.push_back( a );
		v = nodes[v].pred;
	}
	return sp;
}
