#ifndef CUTCALLBACK_H_
#define CUTCALLBACK_H_

#include "Instance.h"
#include "Maxflow.h"
#include <ilcplex/ilocplex.h>

using namespace std;

class CutCallback: public IloCplex::Callback::Function
{

private:

	IloEnv env;
	string cut_type;
	double epsI; 
	double epsO;
	Instance &instance;
	const IloCplex::Callback::Context *context;
	IloBoolVarArray &x;
	IloBoolVarArray &z;
	IloNumArray xsol, zsol;
	
	
	// count cuts
	u_int &u_cuts; // user cuts
	u_int &c_cuts; // candidate cuts
	
	// save arcs, <node1 to node2> by same order as in x here
	list<pair<u_int, u_int>> arcs; 
	
	// capacities of current arcs, same order as arcs
	vector<double> capacities;
	
	// number of usercuts
	u_int n_ucuts;
	double multiplier_ucuts;
	
	
	
	// 1 if violated path identified 0 else
	u_int violated_path;

	// separate directed connection cuts
	void connectionCuts();

	// separate cycle elimination cuts
	void cycleEliminationCuts();


	// SHORTEST PATHS

	struct SPNodeT
	{
		int pred; // -1 when uninitialized
		int pred_arc_id; // -1 when uninitialized
		double weight; // path weight to node
	};

	// result of shortest path computation
	struct SPResultT
	{
		list<u_int> path;
		double weight;
	}sp;
	
	// shortest path structure
	//SPResultT sp;

	// arc weights for shortest path computation
	// (should be set according to current LP solution)
	vector<double> arc_weights;

	// computes a shortest path from source to target according to arc_weights
	// in a directed graph:
	//    number of edges: m -> number of arcs: 2*m
	//    edge (v1,v2) with id <i> -> arc (v1,v2) with id <i>,
	//                                arc (v2,v1) with id <i+m>
	// (returns list of arc ids of a shortest path and the according weight)
	SPResultT shortestPath( u_int source, u_int target );

public:

	CutCallback( string _cut_type, double _epsI, double _epsO, Instance &_instance, IloBoolVarArray &_x, IloBoolVarArray &_z,
				 u_int &_u_cuts, u_int &_c_cuts );
	virtual ~CutCallback();

	void invoke( const IloCplex::Callback::Context &_context );

};

#endif /* CUTCALLBACK_H_ */
