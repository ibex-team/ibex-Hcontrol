//
// Created by joris on 17/11/2021.
//

#ifndef IBEX_HCONTROL_IBEX_EVALMAX_H
#define IBEX_HCONTROL_IBEX_EVALMAX_H

#include <vector>
#include "ibex_Function.h"
#include "ibex_Ctc.h"
#include "ibex_Bsc.h"
#include "ibex_BoxProperties.h"
#include "ibex_Bxp.h"
#include "ibex_Cell.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_UnconstrainedLocalSearch.h"
#include "ibex_LargestFirst.h"
#include "ibex_BxpMinMax.h"
#include "ibex-affine/ibex_AffineEval.h"

namespace ibex {
    class BxpMinMax;
    class BxpMinMaxSub;

    class EvalMax {

    public:

        /* Constructor without constraint on y */
//        EvalMax(Function& f, int nx, int ny); TODO

        /* Constructor with constraint */
        EvalMax(ExtendedSystem& xy_sys, int nx, int ny, Ctc& ctc_xy);

        /* Constructor*/
//    EvalMax(NormalizedSystem& y_sys,Ctc& ctc_xy,UnconstrainedLocalSearch* local_solver,bool csp_actif = false);

        /* Constructor*/
//    EvalMax(NormalizedSystem& y_sys);

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
        Interval eval(IntervalVector& X, double loup = 1e15);
        Interval eval(IntervalVector& X, BoxProperties& prop, double loup = 1e15);
//        bool eval(Cell* x_cell, double loup);

        int trace;
        double timeout;
        int list_elem_max;
        double ext_crit_prob;
        int nb_iter;
        double prec_y;
        bool monitor;
        int local_search_iter;
        bool visit_all;
        ExtendedSystem& xy_sys; // contains constraints on x and y
        double goal_abs_prec; // absolute precision on goal evaluation, stop maximization when reached

    private:
        friend class OptimMinMax;
        Affine2Eval* affine_goal;
        Ctc& ctc_xy; //contractor for constraints on xy
        UnconstrainedLocalSearch *local_solver;
        //double abs_min_prec; // absolute minimum prec bisection on y
        Bsc* bsc; // bisector
        std::vector<Cell*> heap_save;
        bool found_point;
        double time;
        bool csp_actif;
        IntervalVector best_point_eval;


        double save_heap_ub;

        /* contract xy_box and xy_box_ctc w.r.t max_ctc contractor
         * */
//    void contract_best_max_cst( Ctc* max_ctc,IntervalVector* xy_box,IntervalVector* xy_box_ctc,y_heap_elem* elem);



        bool handle_cell(Cell* x_cell, Cell* y_cell, double loup, bool no_stack = false);
        // no stack needed for visit all leaves.

        /**
         * Delete the elements in the save heap
         */
        void delete_save_heap();

        IntervalVector xy_box_hull(const IntervalVector& x_box);

        /**
         * run local search algorithm for a particular x and maximizes over y to provide y_max a local maximum. Objective function is then evaluate at (xbox,max_y) to try to provide a better lower bound.
         * Inputs: x_box: current x box, xy_box: box after contraction w.r.t contraction, loup: current lower upper.
         */
        double local_search_process(const IntervalVector& x_box, const IntervalVector& xy_box, double loup);


        /**
         * return true if the stop criterion is reached
         */
        bool stop_crit_reached(int current_iter,DoubleHeap<Cell> * y_heap,const Interval& fmax) const;

//        void set_y_sol(Vector& start_point);

        /**
         * add elements of Heap_save into y_heap
         */
        void fill_y_heap(DoubleHeap<Cell>& y_heap);

        /**
         * set y part of xy_box with y_box
         */
        static IntervalVector init_xy_box(const IntervalVector& x_box, const IntervalVector& y_box);

        bool handle_constraint(BxpMinMaxSub *data_y, IntervalVector &xy_box, IntervalVector &y_box);
        bool handle_cstfree(IntervalVector& xy_box, Cell * y_cell);

        /**
         * return a feasible point in y_box w.r.t constraints on xy
         */
        IntervalVector get_feasible_point(Cell* x_cell, Cell * y_cell);

        static Interval eval_all(Function* f, const IntervalVector& box);

        /**
         * return 0 if box is non feasible w.r.t constraints on xy, 1 if not known, 2 if box is entirely feasible
         */
        int check_constraints(const IntervalVector& xy_box);

        /**
         * returns a box composed of x_box(not modified) and the middle of y_box, needed for midpoint evaluation
         * Inputs: -xy_box: whole box
         *         -y_box: y box to get the middle
         */
        static IntervalVector get_mid_y(const IntervalVector& x_box, const IntervalVector& y_box);


        /* Default timeout */
        static const double default_timeout;
        static const double default_goal_abs_prec;


//        bool check_already_in(Cell * const y_cell, DoubleHeap<Cell> * y_heap);

    };

    void export_monitor(std::vector<double>* ub, std::vector<double>* lb, std::vector<double>* nbel, std::vector<double>* nbel_save, const IntervalVector& box);
} // end namespace ibex

#endif //IBEX_HCONTROL_IBEX_EVALMAX_H
