//============================================================================
//                                  I B E X
// File        : ibex_EvalMax.cpp
// Author      : Dominique Monnet, Jordan Ninin, Joris Tillet
// License     : See the LICENSE file
// Created     : Oct 1, 2021
//============================================================================


#include <fstream>
#include "ibex_EvalMax.h"
#include "ibex_Timer.h"
#include "ibex_SystemFactory.h"
#include "ibex_CtcIdentity.h"
#include "ibex_NoBisectableVariableException.h"
//using namespace std;
namespace ibex {

//Default parameters for light optim min max solver
const double EvalMax::default_timeout = 20;
const int EvalMax::default_prob_heap = 10 ; //10% to pop second heap in light_solver
const double EvalMax::default_goal_abs_prec = 1e-4;
const double EvalMax::default_goal_rel_prec = 1e-6;
const int EvalMax::default_iter =10000;
const bool EvalMax::default_visit_all = false;
const double EvalMax::default_prec_y= 1.e-14;



EvalMax::EvalMax(IntervalVector& y_box_init, Function &fxy)  :
		id(next_id()),
		trace(false),
		timeout(default_timeout),
		list_elem_max(0),
		nb_iter(default_iter),
		prec_y(default_prec_y),
		monitor(false),
		local_search_iter(2),
		visit_all(default_visit_all),
		goal_abs_prec(default_goal_abs_prec),
		goal_rel_prec(default_goal_rel_prec),
		xy_sys(NULL),
		ctc_xy(NULL),
		minus_goal_y_at_x(NULL),
		local_solver(NULL),
		bsc(new LargestFirst()),
		found_point(false),
		time(0),
		y_box_init(y_box_init),
		crit_heap(default_prob_heap),
		cost1(new CellCostMaxPFub_MinMax(*this)),
		cost2(new CellCostPFlb_MinMax(*this)),
		save_heap_ub(NEG_INFINITY)
{
	SystemFactory fac;
	fac.add_var(fxy.args());
	fac.add_goal(fxy);

	xy_sys = &rec(new System(fac));
	ctc_xy = &rec(new CtcIdentity(xy_sys->box.size()));

	if ( (y_box_init.size()<xy_sys->box.size()) ) {

		// goal function reformulation as min instead of max for local solver
		Array<const ExprNode> args(xy_sys->goal->nb_arg());
		Array<const ExprSymbol> var;
		for(int i = 0;i<xy_sys->goal->nb_arg();i++) {
			const ExprSymbol& a = ExprSymbol::new_(xy_sys->goal->arg(i).dim);
			var.add(a);
			args.set_ref(i,a);
		}
		minus_goal_y_at_x = new Function(var,-(*xy_sys->goal)(args));

		// TODO replace the with a constrained local search (when it will be done)
		local_solver = new UnconstrainedLocalSearch(*minus_goal_y_at_x,IntervalVector(1));

	} else {
		ibex_error("[ibex_EvalMax] -- xy_sys must have a goal function.");
	}
}


EvalMax::EvalMax(IntervalVector& y_box_init, System &xy_sys, Ctc &ctc_xy) :
		id(next_id()),
		trace(false),
		timeout(default_timeout),
		list_elem_max(0),
		nb_iter(default_iter),
		prec_y(default_prec_y),
		monitor(false),
		local_search_iter(2),
		visit_all(default_visit_all),
		goal_abs_prec(default_goal_abs_prec),
		goal_rel_prec(default_goal_rel_prec),
		xy_sys(&xy_sys),
		ctc_xy(&ctc_xy),
		minus_goal_y_at_x(NULL),
		local_solver(NULL),
		bsc(new LargestFirst()),
		found_point(false),
		time(0),
		y_box_init(y_box_init),
		crit_heap(default_prob_heap),
		cost1(new CellCostMaxPFub_MinMax(*this)),
		cost2(new CellCostPFlb_MinMax(*this)),
		save_heap_ub(NEG_INFINITY)
{
	if ((xy_sys.goal !=NULL)|| (y_box_init.size()<xy_sys.box.size()) ) {

		// goal function reformulation as min instead of max for local solver
		Array<const ExprNode> args(xy_sys.goal->nb_arg());
		Array<const ExprSymbol> var;
		for(int i = 0;i<xy_sys.goal->nb_arg();i++) {
			const ExprSymbol& a = ExprSymbol::new_(xy_sys.goal->arg(i).dim);
			var.add(a);
			args.set_ref(i,a);
		}
		minus_goal_y_at_x = new Function(var,-(*xy_sys.goal)(args));

		// TODO replace the with a constrained local search (when it will be done)
		local_solver = new UnconstrainedLocalSearch(*minus_goal_y_at_x,IntervalVector(1));

	} else {
		ibex_error("[ibex_EvalMax] -- xy_sys must have a goal function.");
	}
}

EvalMax::~EvalMax() {
	delete bsc;
    delete local_solver;
    delete minus_goal_y_at_x;
    delete cost1;
    delete cost2;
	delete_save_heap();
}

void EvalMax::add_property(const IntervalVector& init_box, BoxProperties& map) {

	if (!map[BxpMinMax::get_id(*this)]) {
		BxpMinMax* data_x= new BxpMinMax(*this);

		Cell * y_cell = new Cell(y_box_init);

		// add data required by the bisector
		bsc->add_property(y_cell->box,y_cell->prop);
		// add data required by the contractor
		// attention on met dans y_prop des info sur xy_box(init_box|y_box_init) : ctc_xy.add_property(xy_box,y_cell->prop);
		ctc_xy->add_property(xy_sys->box,y_cell->prop);
		// add data required by the buffer
		data_x->y_heap.add_property(y_cell->box,y_cell->prop);
		BxpMinMaxSub * data_y = dynamic_cast<BxpMinMaxSub *>(y_cell->prop[BxpMinMaxSub::get_id(*this)]);
		if (!data_y)
			ibex_error("[ibex_EvalMax] -- null ptr in eval: y_cell element has no BxpMinMaxSub.");
		else
			data_y->feasible = (xy_sys->nb_ctr==0);
		// add the entire domain of y (i.e. y_box_init) in the y_heap
		data_x->y_heap.push(y_cell);

		map.add(data_x);
	}
}

Interval EvalMax::eval(const IntervalVector &x, double loup) {
	BoxProperties prop(x);
	return this->eval(x, prop, loup);
}

//Interval EvalMax::eval(Cell &x, double loup) {
//	return this->eval(x.box, x.prop, loup);
//}

Interval EvalMax::eval(const IntervalVector &x_box, BoxProperties &x_prop, double loup) {
	Timer timer;
	timer.start();

	bool res =optimize(x_box, x_prop, loup);

	time = timer.get_time();
	timer.stop();

	if (res) {
		BxpMinMax *data_x = dynamic_cast<BxpMinMax* >(x_prop[BxpMinMax::get_id(*this)]);
		return data_x->fmax;
	} else {
		return Interval::empty_set();
	}

}

bool EvalMax::optimize(const IntervalVector &x_box, BoxProperties &x_prop, double loup) {

	add_property(x_box, x_prop);

	delete_save_heap();

	found_point  = false;

	BxpMinMax *data_x = dynamic_cast<BxpMinMax* >(x_prop[BxpMinMax::get_id(*this)]);

//	std::cout <<"    DEB "<<data_x->fmax <<std::endl;
//	std::cout<<std::endl<<"*************************"<<std::endl;
//	std::cout<<"get y_heap of size: "<<y_heap.size()<<std::endl;
//	std::cout<<"initial fmax: "<<data_x->fmax<<std::endl;
//	std::cout<<"box: "<<x_box<<std::endl;

	CellHeapMinMax& y_heap = data_x->y_heap; // current cell

	// Define the TimeOut of to compute the bounds of x_box
	Timer timer;
	timer.start();

	// ********** contract x_box with ctc_xy***************
	/* ca n' a rien a faire ici. Ce genre de contraction de x_box doit être dans un contracteur
	IntervalVector xy_box(x_box.size()+y_box_init.size());
	xy_box.put(0,x_box);
	xy_box.put(x_box.size(),y_box_init);
	ctc_xy.contract(xy_box);
	if (xy_box.is_empty()) {
		return false;
	} else { // contract the result on x
		x_box &= xy_box.subvector(0, x_box.size()-1);
	}
	*/
	// ********************************************************

	// *************************************************************************************
	//monitoring variables, used to track upper bound, lower bound, number of elem in y_heap and heap_save at each iteration
	std::vector<double> ub, lb, nbel, nbel_save;

	save_heap_ub = NEG_INFINITY;
//    std::cout<<"initial fmax: "<<data_x->fmax<<std::endl;


	// *********** loop ********************
	try {
		int current_iter = 0;
		while(!stop_crit_reached(current_iter, *data_x) ) {

			found_point  = false;

			Cell* y_cell = y_heap.pop(); // we extract an element with critprob probability to take it according to the first crit
			current_iter++;

//			     std::cout << *y_cell << std::endl;
//			     std::cout<<"current_iter: "<<current_iter<<std::endl;
			if (((list_elem_max != 0 && ((y_heap.size() + heap_save.size())>list_elem_max)) || (y_cell->box.max_diam())<=prec_y) ) {
				// continue to evaluate cells of y_heap without increasing size of list, need to do it else nothing happend if list already reached max size
//				std::cout<< " TOO small : "<< y_cell->box  << std::endl;
				bool res = handle_cell(x_box, data_x, y_cell, loup);
				if (!res) return false; // x_cell has been deleted

			} else {
				try {
					std::pair<Cell*,Cell*> subcells_pair = bsc->bisect(*y_cell);// bisect tmp_cell into 2 subcells
					delete y_cell;

//					std::cout<< " G : "<< subcells_pair.first->box <<"  \t  D : "<< subcells_pair.second->box << std::endl;

					// std::cout<<"handle first cell in light optim"<<std::endl;
					bool res = handle_cell(x_box, data_x, subcells_pair.first, loup);
					// std::cout<<"first cell handled"<<std::endl;
					if (!res) { // x_cell has been deleted
						delete subcells_pair.second;
						return false;
					}

					// std::cout<<"handle second cell in light optim"<<std::endl;
					res = handle_cell( x_box, data_x, subcells_pair.second,loup);
					// std::cout<<"second cell handled"<<std::endl;
					if (!res) { // x_cell has been deleted
						return false;
					}
				}
				catch (NoBisectableVariableException& ) {
//					std::cout<< " NO cut : "<< y_cell->box  << std::endl;

					bool res = handle_cell(x_box, data_x, y_cell,loup);
					if (!res) return false; // x_cell has been deleted
				}

			}

			// Reduction de la taille de y_heap
			if (found_point) {
				//elimine cost1 > loup
				//elimine  ( -data->maxfxy.ub() ) > (-(data_x->fmax.lb()))
				//elimine  tous les y, tel que : data_y->maxfxy.ub() < (data_x->fmax.lb())
				//elimine tous les y, dont on est sur qu'ils ne contiennent pas le max (peu importe le x).
				y_heap.contract(-(data_x->fmax.lb()));
			}

			///////////////////////////////
			if (monitor) {
				if (!y_heap.empty()) {
					lb.push_back(data_x->fmax.lb());
					BxpMinMaxSub * top1_minmax = dynamic_cast<BxpMinMaxSub *>(y_heap.top1()->prop[BxpMinMaxSub::get_id(*this)]);
					if (!top1_minmax)
						ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMaxSub.");
					if (save_heap_ub < top1_minmax->maxfxy.ub())
						ub.push_back(top1_minmax->maxfxy.ub());
					else
						ub.push_back(save_heap_ub);
					nbel.push_back(y_heap.size());
					nbel_save.push_back(heap_save.size());
				}
			}
			///////////////////////////

			// std::cout<<"nb iter: "<<current_iter<<std::endl;
			timer.check(timeout);

		}

	}
	catch (TimeOutException& ) { }

	//  ***************   force visit of all elements of the heap *************************
	if (visit_all) {
//		std::cout<<"  visit_all  y_heap detail: "<<std::endl;
		while (!y_heap.empty()) {
			Cell * y_cell = y_heap.pop();

//			BxpMinMaxSub * data_y = dynamic_cast<BxpMinMaxSub *>(y_cell->prop[BxpMinMaxSub::get_id(*this)]);
//			std::cout<<"   \t \t      "<<y_cell->box<< "  fmax= " << data_y->maxfxy << std::endl;

			bool res = handle_cell(x_box, data_x, y_cell, loup, true);
			if (!res) return false;

			///////////////////////////////
			if (monitor) {
				if (!y_heap.empty()) {
					lb.push_back(data_x->fmax.lb());
					BxpMinMaxSub * top1_minmax = dynamic_cast<BxpMinMaxSub *>(y_heap.top1()->prop[BxpMinMaxSub::get_id(*this)]);
					if (!top1_minmax)
						ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMaxSub.");
					if (save_heap_ub < top1_minmax->maxfxy.ub())
						ub.push_back(top1_minmax->maxfxy.ub());
					else
						ub.push_back(save_heap_ub);
					nbel.push_back(y_heap.size());
					nbel_save.push_back(heap_save.size());
				}
			}
			//////////////////////////////
		}

	}

	// ******************  merge save_heap in y_heap  *************************************
	fill_y_heap(y_heap);

	if (y_heap.empty()) {// cela veut dire que y_box_init ne verifie pas les contraintes
//		std::cout <<"   WARNING    OUT 3 "<< std::endl;
		return false;
	}


	// ********** Update the lower and upper bound of "max f(x,y_heap)"  **********************
	BxpMinMaxSub * top1_minmax = dynamic_cast<BxpMinMaxSub *>(y_heap.top1()->prop[BxpMinMaxSub::get_id(*this)]);
	if (!top1_minmax) ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMaxSub.");
	// get the upper bound of max f(x,y_heap)
	double new_fmax_ub = top1_minmax->maxfxy.ub();

	BxpMinMaxSub * top2_minmax = dynamic_cast<BxpMinMaxSub *>(y_heap.top2()->prop[BxpMinMaxSub::get_id(*this)]);
	if (!top2_minmax) ibex_error("[ibex_EvalMax] -- null ptr in eval: top2 element has no BxpMinMaxSub.");
	// get the lower bound of max f(x,y_heap)
	double new_fmax_lb = top2_minmax->maxfxy.lb();


	if (new_fmax_ub< new_fmax_lb) {
		ibex_error("[ibex_EvalMax] --  error, please report this bug.");
	}

	data_x->fmax &= Interval(new_fmax_lb, new_fmax_ub);

	// ***********************Reduction de la taille de y_heap*********************************
	//elimine  tous les y, tel que : data_y->maxfxy.ub() < (data_x->fmax.lb())
	//elimine tous les y, dont on est sur qu'il ne contienne pas le max (peu importe le x).
	y_heap.contract(-(data_x->fmax.lb()));

	if (y_heap.empty()) {// cela veut dire que y_box_init ne verifie pas les contraintes
//		std::cout <<"   WARNING    OUT 4 "<< std::endl;
		return false;
	}

	//************************************************************************************

	if (data_x->fmax.is_empty() || data_x->fmax.lb() > loup) {  // on ne peut pas trouver mieux que "loup"
		//     std::cout <<"       OUT 4 "<<std::endl;
		return false;
	}

	///////////////////////////////
	if (monitor) {
		lb.push_back(data_x->fmax.lb());
		ub.push_back(data_x->fmax.ub());
		nbel.push_back(y_heap.size());
		nbel_save.push_back(heap_save.size());
		export_monitor(&ub, &lb, &nbel, &nbel_save, x_box);
	}
	///////////////////////////////

	//                std::cout <<"    FIN "<<data_x->fmax <<std::endl;
	return true;
}

//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

bool EvalMax::handle_cell(const IntervalVector& x_box, BxpMinMax* data_x , Cell* y_cell, double loup, bool no_stack) {

	IntervalVector xy_box(x_box.size()+y_cell->box.size());
	xy_box.put(0, x_box);
	xy_box.put(x_box.size(), y_cell->box);

	// recuperer les data
	BxpMinMaxSub * data_y = dynamic_cast<BxpMinMaxSub *>(y_cell->prop[BxpMinMaxSub::get_id(*this)]);

	// ************************ Check constraints **************************************
	if (!(data_y->feasible) ){
		if (handle_constraint(xy_box, y_cell->box, y_cell->prop)) {
			delete y_cell;
			return true;
		}
	} else {
		handle_ctrfree(xy_box, y_cell->box);
	}

	//************** mid point test and  local optim test point  *****************
	// To improve the lower bound of fmax
	std::pair<double , IntervalVector >  res = get_best_lb(x_box,y_cell->box, data_y->feasible);
	double best_lb = res.first;

	if (best_lb > data_y->maxfxy.lb()) {
		data_y->maxfxy &= Interval(best_lb,POS_INFINITY);
	}
	if (data_x->fmax.lb() < best_lb) {
		//  data_x->best_sol = mid_y_box.mid();
		found_point = true;
		data_x->fmax &= Interval(best_lb,POS_INFINITY);; // yes we found a feasible solution for all x
		data_x->best_sol = res.second;
	}

	if ( data_x->fmax.is_empty() || (best_lb > data_x->fmax.ub())) {
		ibex_error("[ibex_EvalMax] -- error, lb >ub from local search. Please report this bug.");
	}
	if (loup < best_lb) {
		// there exists y such as constraint is respected and f(x,y)>best max, the max on x will be worst than the best known solution
		delete y_cell;
		return false; // no need to go further, x_box does not contain the solution
	}


	//************ part below add a contraction w.r.t f(x,y)<best_max, this part may not be efficient on every problem ******************************
	if (data_y->feasible && (data_x->fmax.ub()>loup) ) {
		IntervalVector xy_box_mem(xy_box);
		// std::cout<<"contraction in light solver, init box: "<<mid_y_box<<std::endl;
		xy_sys->goal->backward(Interval(NEG_INFINITY,loup),xy_box);

		if(xy_box.is_empty()) {
			delete y_cell;
			// std::cout<<"delete x cell because empty contraction"<<std::endl;
			return false;
		}

		// y contracted => E y, for all x f(x,y)>loup, => x deleted
		for (int i=x_box.size(); i<xy_box.size(); i++) {
			if (xy_box[i] != xy_box_mem[i]) {
				// std::cout<<"delete x cell because contraction on y"<<std::endl;
				delete y_cell;
				return false;
			}
		}
		// no contraction on y but contraction on  x: keep contraction on x
		// NON on ne contracte pas x_box, ici on fait que du eval. On mets cette contraction dans un contracteur ailleur
		//for (int k=0; k<x_box.size(); k++) {
		//	x_box[k] &= xy_box[k];
		//}
	}
	//********************************************

	// Update the lower and upper bound on y
	// objective function evaluation
	Interval tmp_eval_goal = xy_sys->goal->eval(xy_box);
	if (tmp_eval_goal.ub()>data_y->maxfxy.ub()) {
		ibex_error("[ibex_EvalMax] -- get worst upper bound, should not happen due to monotonicity of ifunc. ");
	}

//	std::cout<<"    dealing with box: "<<xy_box<<std::endl;
//	std::cout<<"    previous value of maxfxy:  "<<data_y->maxfxy<<std::endl;
//    std::cout<<"    eval goal  res: "<<tmp_eval_goal<<std::endl;

	data_y->maxfxy &= tmp_eval_goal;

    if (data_y->maxfxy.is_empty() || data_x->fmax.lb() > data_y->maxfxy.ub()) {  // y_box cannot contains max f(x,y)
		delete y_cell;
		return true;
	}

	// check if it is possible to find a better solution than those already found on x
	if ((data_y->maxfxy.lb() > loup) && (data_y->feasible)) {
		// box verified condition and eval is above best max, x box does not contain the solution
		// I think this case was already check with the mid-point in get_best_lb.
		delete y_cell;
		//std::cout <<"           OUT 7 "<<std::endl;
		return false; // no need to go further, x_box does not contains the solution
	}

	//*************************************************
	// store y_cell
	if (y_cell->box.max_diam() <= prec_y) {
//		            std::cout<<"y_cell pushed in heap_save, box: "<<y_cell->box<<" pf: "<<data_y->maxfxy<<" pu: "<<data_y->feasible<<std::endl;
		save_heap_ub = save_heap_ub<data_y->maxfxy.ub()? data_y->maxfxy.ub() : save_heap_ub;
		heap_save.push_back(y_cell);
	}
	else {
		if (!no_stack) {
			data_x->y_heap.push(y_cell);
		}
		else
			heap_save.push_back(y_cell);
	}
	return true;
}

//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

bool EvalMax::handle_constraint( IntervalVector& xy_box, IntervalVector& y_box, BoxProperties& y_prop ) {
	//    std::cout<<"handle constraint on xy, xy init: "<<xy_box<<std::endl;

	BxpMinMaxSub * data_y = dynamic_cast<BxpMinMaxSub *>(y_prop[BxpMinMaxSub::get_id(*this)]);
	if (!(data_y->feasible) )  {
		switch (check_constraints(xy_box)){
		case 2: { // all the constraints are satisfied
			data_y->feasible=true;
			break;
		}
		case 0: { // One constraint is false
			// constraint on x and y not respected, move on.
			return true;
		}
		default: // nothing to do
			break;
		}
	}
	if (!(data_y->feasible) )  {// there is a constraint on x and y
		ContractContext y_context(y_prop);
		ctc_xy->contract(xy_box, y_context); // TODO regarder ce que l on peut faire de mieux avec ce context
		if (xy_box.is_empty()) { // constraint on x and y not respected, move on. Delete y
			return true;
		} else {
			// on peut propager la contraction sur le y, mais pas sur le x
			for (int k=0; k<y_box.size(); k++) {
				y_box[k] &= xy_box[xy_box.size()-y_box.size()+k];
			}
		}
	} else {
		handle_ctrfree(xy_box,y_box);
	}
	return false;
}



void EvalMax::handle_ctrfree( IntervalVector& xy_box, IntervalVector& y_box ) {
	//    std::cout<<"init box: "<<*xy_box<<std::endl;
	IntervalVector grad(xy_box.size());
	xy_sys->goal->gradient(xy_box,grad);
	int k =0;
	for (int i=xy_box.size()-y_box.size(); i<xy_box.size(); i++) {
		if (grad[i].lb() > 0) {
			xy_box[i] = xy_box[i].ub();
			y_box[k] = xy_box[i];
		}
		if (grad[i].ub() < 0) {
			//                        (xy_box)[i] = Interval((xy_box)[i].lb(),(xy_box)[i].lb()+1.e-15);
			xy_box[i] = xy_box[i].lb();
			y_box[k] = xy_box[i];
		}
		k++;
	}
	//    std::cout<<"final box: "<<*xy_box<<std::endl;
	//    std::cout<<"free cst contraction done, contracted box: "<<*xy_box<<std::endl;
}

//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________


//mid point test and  local optim test : to find a better lower bound of fmax
std::pair<double, IntervalVector> EvalMax::get_best_lb(const IntervalVector& x_box, const IntervalVector& y_box, bool y_feasible)  {
	double best_lb= NEG_INFINITY;
	IntervalVector try_y_point(x_box.size()+y_box.size());
	try_y_point.put(0,x_box);

	IntervalVector best_point(x_box.size()+y_box.size());
	best_point.set_empty();

	IntervalVector xy_box(x_box.size()+y_box.size());
	xy_box.put(0,x_box);
	xy_box.put(x_box.size(),y_box);
	Vector final_point = xy_box.mid();

	// initialise the box constraint of the local minimizer
	local_solver->set_box(xy_box);

	for (int k=0; k<local_search_iter; k++)  {

		if (k==0) {
			// first try the middle
			// get the box (x,mid(y))
			for (int i=0; i<y_box.size(); i++) {
				try_y_point[x_box.size()+i] = y_box[i].mid();
			}

		} else {
			// try with a local minimizer
			Vector xy_rand = xy_box.random(); // pick a random x

			//std::cout<<"             box: "<< xy_box<<" random start: "<<xy_rand<<std::endl;
			local_solver->minimize(xy_rand,final_point,1e-3); // maximizes
			//std::cout<<"             local solution: "<<final_point<<std::endl;

			// update current y box in the try_y_point for evaluation at a point y given by local optim
			for (int k =x_box.size(); k<xy_box.size(); k++) {
				try_y_point[k] = final_point[k];
			}
		}

		//check if the constraint is respected
		bool ctr_ok = y_feasible;
		if (!ctr_ok) { // constraint on xy exist and is not proved to be satisfied
			ctr_ok = (check_constraints(try_y_point) ==2);
		}
		if (ctr_ok) {
			// found y such as xy constraint is respected

			// x y constraint respected for all x_box and a point y, try_y_point is a candidate for evaluation
			Interval midres = xy_sys->goal->eval(try_y_point);
			if (midres.lb() > best_lb) {
				best_lb= midres.lb();
				best_point= try_y_point;
			}
		}
	}
	return std::pair<double, IntervalVector>(best_lb,best_point);

}

//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

int EvalMax::check_constraints(const IntervalVector& xy_box) {
	int out = 2;
	if (xy_sys->nb_ctr>0) {
		bool toto= true;
		IntervalVector res = xy_sys->f_ctrs.eval_vector(xy_box);
		bool val;
		for (int c=0; c<xy_sys->f_ctrs.image_dim(); c++) {
			switch (xy_sys->ops[c]) {
			case LT:  val=res[c].lb()>0; break;
			case LEQ: val=res[c].lb()>0; break;
			case EQ:  val=(!res[c].contains(0)); break;
			case GEQ: val=res[c].ub()<0; break;
			case GT:  val=res[c].ub()<0; break;
			}
			if (val)
				return 0;
			else if (toto) {
				switch (xy_sys->ops[c]) {
				case LT:  val=res[c].ub()>=0; break;
				case LEQ: val=res[c].ub()>0; break;
				case EQ:  val=res[c]!=Interval::zero(); break;
				case GEQ: val=res[c].lb()<0; break;
				case GT:  val=res[c].lb()<=0; break;
				}
				toto = false;
				out = 1;
			}
		}
	}
	return out;
}


//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

void EvalMax::delete_save_heap() {
	while (!heap_save.empty()) {
		delete heap_save.back();
		heap_save.pop_back();
	}
}




bool EvalMax::stop_crit_reached(int current_iter, const BxpMinMax& data_x) const {
	if (nb_iter !=0 && current_iter >= nb_iter)  {// nb_iter ==  0 implies minimum precision required (may be mid point x case)
//		std::cout<<"STOP: nb_iter max reached"<<std::endl;
		return true;
	}
	if(list_elem_max != 0 && data_x.y_heap.size() + heap_save.size()>list_elem_max) {
//		std::cout<<"STOP: list reached max elem"<<std::endl;
		return true;
	}
	if (data_x.y_heap.empty()) {
//		std::cout<<"STOP: empty buffer"<<std::endl;
		return true;
	}

	if (data_x.fmax.diam()<goal_abs_prec) { //desired precision on global maximum enclosure reached
//		std::cout<<"STOP: absolute precision reached "<<std::endl;
		return true;
	}
	if (data_x.fmax.diam()<goal_rel_prec*data_x.fmax.mid()) { //desired precision on global maximum enclosure reached
//		std::cout<<"STOP: relative precision reached "<<std::endl;
		return true;
	}
	return false;
}




void EvalMax::fill_y_heap(DoubleHeap<Cell>& y_heap) {
	Cell* tmp_cell;
	//        int nb_pt(0);
	while (!heap_save.empty()) {
		// push all boxes of heap_save in y_heap
		tmp_cell = heap_save.back();
		//std::cout << tmp_cell->box << std::endl;

		y_heap.push(tmp_cell);
		heap_save.pop_back();
	}
}



void export_monitor(std::vector<double>* ub, std::vector<double>* lb, std::vector<double>* nbel, std::vector<double>* nbel_save, const IntervalVector& box) {
	//    std::ofstream out;
	std::ofstream out("log.txt", std::ios::app);
	//    std::string output_name = "log.txt";
	//    out.open("log.txt");
	out <<"--------------------------------"<<std::endl;
    out <<"x box: "<<box<<std::endl;
	out <<"--------------------------------"<<std::endl;
	out <<"  lb \t  ub \t heap \t save  "<<std::endl;
	for (size_t i = 0; i<ub->size(); i++) {
		out<<lb->at(i)<<"\t"<<ub->at(i)<<"\t"<<nbel->at(i)<<"\t"<<nbel_save->at(i)<<std::endl;
	}
	out.close();
}

} // end namespace ibex
