//
// Created by joris on 17/11/2021.
//

#ifndef __IBEX_HCONTROL_IBEX_EVALMAX_H
#define __IBEX_HCONTROL_IBEX_EVALMAX_H

#include <vector>
#include "ibex_Interval.h"
#include "ibex_Function.h"
#include "ibex_Ctc.h"
#include "ibex_Bsc.h"
#include "ibex_Memory.h"
#include "ibex_BoxProperties.h"
#include "ibex_Bxp.h"
#include "ibex_Cell.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_UnconstrainedLocalSearch.h"
#include "ibex_LargestFirst.h"
#include "ibex_ContractContext.h"
#include "ibex_CellFuncMinMax.h"
//#include "ibex_CellMinMaxHeap.h"
#include "ibex_BxpMinMax.h"

//#include "ibex-affine/ibex_AffineEval.h"

namespace ibex {
//class BxpMinMax;
//class BxpMinMaxSub;
class CellHeapMinMax;
//class CellCostMaxPFub_MinMax ;
//class CellCostPFlb_MinMax;


class EvalMax : protected Memory {

public:

	/* Constructor without constraint on y */
	EvalMax(IntervalVector& y_box_init, Function& fxy);

	/* Constructor with constraint */
	EvalMax(IntervalVector& y_box_init, System& xy_sys, Ctc& ctc_xy);

	/* Destructor */
	~EvalMax();

	/* returns an enclosure of the maximum of the objective function: max f(x,y)
	 * modifies y_heap inherited from father box
	 * This function works as a classic B&B algorithm, with an initial stack of box instead of a initial box.
	 * Inputs: -y_heap: heap of cells containing y boxes sorted w.r.t the ub of eval, the top cell has the greatest ub eval.
	 *         -x_box: current x box, needed to evaluate f(x,y)
	 *         -objective_function: f
	 *         -nb_iter: number of times boxes of y_heap are cut
	 *         -best_max: minimum found for max(f), best current solution
	 *         -fmax: enclosure of max(f(x_box),y_heap), this result is inherited from early computation with a box containing x_box
	 *         -min_prec: minimum size of boxes in y_heap
	 *         -is_midp: true if optimize run with x midpoint eval, false else
	 * */
	Interval eval(const IntervalVector& X, double loup = POS_INFINITY);
	Interval eval(const IntervalVector& X, BoxProperties& prop, double loup = POS_INFINITY);

	/**
	 * Allows to add the properties data required
	 * by this MinMax optimizer to the root cell
	 */
	void add_property(const IntervalVector& init_box, BoxProperties& map);

	int get_size() const;
	const IntervalVector& get_best_point_eval() const;
	double get_goal_abs_prec() const;
	void set_goal_abs_prec(double goalAbsPrec) ;
	double get_goal_rel_prec() const ;
	void set_goal_rel_prec(double goalRelPrec) ;
	unsigned int get_list_elem_max() const ;
	void set_list_elem_max(unsigned int listElemMax);
	int get_local_search_iter() const ;
	void set_local_search_iter(int localSearchIter) ;
	bool is_monitor() const ;
	void set_monitor(bool monitor);
	int get_nb_iter() const ;
	void set_nb_iter(int nbIter);
	double get_prec_y() const ;
	void set_prec_y(double precY);
	double get_time() const ;
	double get_timeout() const;
	void set_timeout(double timeout);
	int get_trace() const ;
	void set_trace(int trace);
	bool is_visit_all() const ;
	void set_visit_all(bool visitAll);
	const IntervalVector& get_y_box_init() const;

	/**
	 * \brief Identifying number.
	 */
	const long id;


private:

	int trace;
	double timeout;
	unsigned int  list_elem_max;
//	double ext_crit_prob;
	int nb_iter;
	double prec_y;
	bool monitor;
	int local_search_iter;
	bool visit_all;

	double goal_abs_prec; // absolute precision on goal evaluation, stop maximization when reached
	double goal_rel_prec; // absolute precision on goal evaluation, stop maximization when reached

	System* xy_sys; // contains constraints on x and y

	Ctc* ctc_xy; //contractor for constraints on xy

	Function* minus_goal_y_at_x; // goal function f becomes -f to solve a minimization problem over y at a fixed x
	UnconstrainedLocalSearch *local_solver;

	Bsc* bsc; // bisector

	std::vector<Cell*> heap_save;

	bool found_point;

	double time;

	IntervalVector y_box_init;

	// To create the y_heap of every BxpMinMax
	friend class CellHeapMinMax;
	int crit_heap;
	CellCostMaxPFub_MinMax* cost1;
	CellCostPFlb_MinMax* cost2;

	double save_heap_ub;


	bool optimize(const IntervalVector& X, BoxProperties& prop, double loup);

	bool handle_cell(const IntervalVector& x, BxpMinMax* data_x , Cell* y_cell, double loup, bool no_stack = false);
	// no_stack=true =>  visit all leaves in y_heap.

	bool handle_constraint( IntervalVector &xy_box, IntervalVector &y_box, BoxProperties& y_prop);
	void handle_ctrfree( IntervalVector& xy_box,  IntervalVector &y_box);

	/**
	 * Delete the elements in the save heap
	 */
	void delete_save_heap();

	/**
	 * run local search algorithm for a particular x and maximizes over y to provide y_max a local maximum.
	 * Objective function is then evaluate at (x_box,max_y) to try to provide a better lower bound.
	 * update x_data->fmax  and y_data->pf
	 */
	//double local_search_process(const IntervalVector& x_box, const IntervalVector& xy_box, double loup);
	std::pair<double, IntervalVector> get_best_lb(const IntervalVector& x_box, const IntervalVector& y_box, bool y_feasible);

	/**
	 * return true if the stop criterion is reached
	 */
	bool stop_crit_reached(int current_iter, const BxpMinMax& data_x) const;

	/**
	 * add elements of Heap_save into y_heap
	 */
	void fill_y_heap(DoubleHeap<Cell>& y_heap);

	/**
	 * return 0 if box is non feasible w.r.t constraints on xy, 1 if not known, 2 if box is entirely feasible
	 */
	int check_constraints(const IntervalVector& xy_box);

    //Default parameters for light optim min max solver
	static const double default_timeout;
	static const double default_goal_abs_prec;
	static const double default_goal_rel_prec;
    static const int default_iter;
    static const int default_prob_heap;
    static const bool default_visit_all;
    static const double default_prec_y;
	//        bool check_already_in(Cell * const y_cell, DoubleHeap<Cell> * y_heap);

};

void export_monitor(std::vector<double>* ub, std::vector<double>* lb, std::vector<double>* nbel, std::vector<double>* nbel_save, const IntervalVector& box);


inline double EvalMax::get_goal_abs_prec() const {
	return goal_abs_prec;
}

inline void EvalMax::set_goal_abs_prec(double goalAbsPrec) {
	goal_abs_prec = goalAbsPrec;
}

inline double EvalMax::get_goal_rel_prec() const {
	return goal_rel_prec;
}

inline void EvalMax::set_goal_rel_prec(double goalRelPrec) {
	goal_rel_prec = goalRelPrec;
}

inline unsigned int EvalMax::get_list_elem_max() const {
	return list_elem_max;
}

inline void EvalMax::set_list_elem_max(unsigned int listElemMax) {
	list_elem_max = listElemMax;
}

inline int EvalMax::get_local_search_iter() const {
	return local_search_iter;
}

inline void EvalMax::set_local_search_iter(int localSearchIter) {
	local_search_iter = localSearchIter;
}

inline bool EvalMax::is_monitor() const {
	return monitor;
}

inline void EvalMax::set_monitor(bool monitor) {
	this->monitor = monitor;
}

inline int EvalMax::get_nb_iter() const {
	return nb_iter;
}

inline void EvalMax::set_nb_iter(int nbIter) {
	nb_iter = nbIter;
}

inline double EvalMax::get_prec_y() const {
	return prec_y;
}

inline void EvalMax::set_prec_y(double precY) {
	if (precY>default_prec_y) {
		prec_y = precY;
	} else {
		prec_y = default_prec_y;
	}
}

inline double EvalMax::get_time() const {
	return time;
}

inline double EvalMax::get_timeout() const {
	return timeout;
}

inline void EvalMax::set_timeout(double timeout) {
	this->timeout = timeout;
}

inline int EvalMax::get_trace() const {
	return trace;
}

inline void EvalMax::set_trace(int trace) {
	this->trace = trace;
}

inline bool EvalMax::is_visit_all() const {
	return visit_all;
}

inline void EvalMax::set_visit_all(bool visitAll) {
	visit_all = visitAll;
}

inline const IntervalVector& EvalMax::get_y_box_init() const {
	return y_box_init;
}

inline int EvalMax::get_size() const {
	return xy_sys->box.size();
}


} // end namespace ibex

#endif //IBEX_HCONTROL_IBEX_EVALMAX_H
