//============================================================================
//                                  I B E X
// File        : ibex_OptimMinMax.cpp
// Author      : Dominique Monnet, Jordan Ninin
// License     : See the LICENSE file
// Created     : Oct 1, 2016
//============================================================================

#include "ibex_OptimMinMax.h"
#include "ibex_LargestFirst.h"
#include "ibex_Timer.h"
#include <stdio.h>
#include "ibex_DataMinMax.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_SystemFactory.h"

using namespace std;

namespace ibex {

//********* Default parameters for light optim minmax solver ************
const int OptimMinMax::default_iter = 10;
const int OptimMinMax::default_list_rate = 0;
const double OptimMinMax::default_min_prec_coef = 10;
const int OptimMinMax::default_list_elem_absolute_max = 500;
const int OptimMinMax::default_prob_heap = 10; //10% to pop second heap in light_solver
const int OptimMinMax::default_local_iter = 0;
const bool OptimMinMax::default_visit_all = false;

//********* Default parameters for light local solver ************

const double OptimMinMax::default_nb_sols = 5;
const double OptimMinMax::default_min_acpt_diam=1e-2;
const int OptimMinMax::default_nb_sivia_iter=0;
const int OptimMinMax::default_nb_optim_iter=0;
const double OptimMinMax::default_y_sol_radius = 0.15;
const double OptimMinMax::default_reg_acpt_error = 1e-1;

const int OptimMinMax::default_nb_point = 1;
const int OptimMinMax::default_eval_period = 1;
const double OptimMinMax::default_perf_thresh = 0.3;

// Csp default parameters for light solver
const int OptimMinMax::default_iter_csp = 10;
const int OptimMinMax::default_list_rate_csp = 0;
const double OptimMinMax::default_min_prec_coef_csp = 100;
const int OptimMinMax::default_list_elem_absolute_max_csp = 100;
const int OptimMinMax::default_prob_heap_csp = 0; //0% to pop second heap in light_solver
const int OptimMinMax::default_local_iter_csp = 0; //10% to pop second heap in light_solver
const bool OptimMinMax::default_visit_all_csp = false;



OptimMinMax::OptimMinMax(NormalizedSystem& x_sys,NormalizedSystem& xy_sys, Ctc& x_ctc,Ctc& xy_ctc,double prec_x,double prec_y,double goal_rel_prec):
                                                                Optim(x_sys.nb_var, new CellDoubleHeap(*new CellCostFmaxlb_opt(), *new CellCostFmaxub_opt()),
                                                                                prec_x, goal_rel_prec, 0, 1), // attention meme precision en relatif et en absolue
                                                                                propag(true),
                                                                                x_box_init(x_sys.box), y_box_init(xy_sys.box.subvector(x_sys.nb_var, xy_sys.nb_var-1)), trace_freq(10000),
                                                                                lsolve(xy_sys,xy_ctc,NULL),
                                                                                loc_solve(xy_sys,NULL,NULL,x_ctc,x_sys.box.size(),xy_sys.box.size()-x_sys.box.size(),false),
                                                                                x_ctc(x_ctc),x_sys(x_sys),
                                                                                bsc(new LargestFirst()),
                                                                                prec_y(prec_y),
                                                                                iter(default_iter),
                                                                                list_rate(default_list_rate),
                                                                                min_prec_coef(default_min_prec_coef),
                                                                                critpr(default_prob_heap),
                                                                                list_elem_absolute_max(default_list_elem_absolute_max),
                                                                                local_iter(0),
                                                                                fa_lsolve(xy_sys,xy_ctc,NULL,true), // useless if no fa cst but need to construct it...
                                                                                fa_loc_solve(xy_sys,NULL,NULL,x_ctc,x_sys.box.size(),xy_sys.box.size()-x_sys.box.size(),true),
                                                                                fa_y_cst(false),
                                                                                y_box_init_fa(IntervalVector(1)),
                                                                                monitor(false),
                                                                                only_csp(false),
                                                                                monitor_csp(false),
                                                                                visit_all(default_visit_all),
                                                                                nb_point(default_nb_point),
                                                                                eval_period(default_eval_period),
                                                                                nb_sols(default_nb_sols),
                                                                                min_acpt_diam(default_min_acpt_diam),
                                                                                nb_sivia_iter(default_nb_sivia_iter),
                                                                                nb_optim_iter(default_nb_optim_iter),
                                                                                y_sol_radius(default_y_sol_radius),
                                                                                reg_acpt_error(default_reg_acpt_error),
                                                                                perf_thresh(default_perf_thresh),
                                                                                min_goal(x_sys.goal != NULL)

{
    if(!min_goal && xy_sys.goal !=NULL) {

    // goal function reformulation as min instead of max for local solver
        Array<const ExprNode> args(xy_sys.goal->nb_arg());
        Array<const ExprSymbol> var;
        for(int i = 0;i<xy_sys.goal->nb_arg();i++) {
            const ExprSymbol& a = ExprSymbol::new_(xy_sys.goal->arg(i).dim);
            var.add(a);
            args.set_ref(i,a);
        }
        minus_goal_y_at_x = new Function(var,-(*xy_sys.goal)(args));

//        minus_goal_csp = minus_goal;// does not matter, fa solver not used

        local_search = new UnconstrainedLocalSearch(*minus_goal_y_at_x,IntervalVector(1));
        lsolve.local_solver = local_search;
        loc_solve.local_solver_over_y = local_search;
        loc_solve.local_solver_over_x = new UnconstrainedLocalSearch(*xy_sys.goal,IntervalVector(1));

        //create affine eval of goal function
        Affine2Eval* aff_eval = new Affine2Eval(*(xy_sys.goal));
        lsolve.affine_goal = aff_eval;
    //    fa_lsolve(xy_sys,xy_ctc,UnconstrainedLocalSearch(minus_goal_csp,IntervalVector(1)),true); // useless if no fa cst but need to construct it...

        lsolve.goal_abs_prec = goal_rel_prec/100; // set goal prec of maximization problem lower than minimization
        fa_lsolve.goal_abs_prec = 1e-2;
    }





};

OptimMinMax::OptimMinMax(NormalizedSystem& x_sys,NormalizedSystem& xy_sys,NormalizedSystem& max_fa_y_cst_sys, Ctc& x_ctc,Ctc& xy_ctc,Ctc& y_fa_ctc,
                double prec_x,double prec_y,double goal_rel_prec,double fa_cst_prec):
                                                Optim(x_sys.nb_var, new CellDoubleHeap(*new CellCostFmaxlb_opt(), *new CellCostFmaxub_opt()),
                                                                prec_x, goal_rel_prec, goal_rel_prec, 1), // attention meme precision en relatif et en absolue
                                                                propag(true),
                                                                x_box_init(x_sys.box), y_box_init(xy_sys.box.subvector(x_sys.nb_var, xy_sys.nb_var-1)), trace_freq(10000),
                                                                x_ctc(x_ctc),x_sys(x_sys),
                                                                lsolve(xy_sys,xy_ctc,NULL),
                                                                loc_solve(xy_sys,NULL,NULL,x_ctc,x_sys.box.size(),xy_sys.box.size()-x_sys.box.size(),false),
                                                                bsc(new LargestFirst()),
                                                                prec_y(prec_y),
                                                                iter(default_iter),
                                                                list_rate(default_list_rate),
                                                                min_prec_coef(default_min_prec_coef),
                                                                critpr(default_prob_heap),
                                                                list_elem_absolute_max(default_list_elem_absolute_max),
                                                                local_iter(default_local_iter),
                                                                prec_fa_y(fa_cst_prec),
                                                                fa_lsolve(max_fa_y_cst_sys,y_fa_ctc,NULL,true), // construct light solver for for all y cst with fake contractor, contractor useless since no constraint on xy only the objective function matter
                                                                fa_loc_solve(xy_sys,NULL,NULL,x_ctc,x_sys.box.size(),xy_sys.box.size()-x_sys.box.size(),true),
                                                                fa_y_cst(true),
                                                                y_box_init_fa(max_fa_y_cst_sys.box.subvector(x_sys.nb_var, max_fa_y_cst_sys.nb_var-1)),
                                                                iter_csp(default_iter_csp),
                                                                list_rate_csp(default_list_rate_csp),
                                                                min_prec_coef_csp(default_min_prec_coef_csp),
                                                                critpr_csp(default_prob_heap_csp),
                                                                list_elem_absolute_max_csp(default_list_elem_absolute_max_csp),
                                                                local_iter_csp(default_local_iter_csp),
                                                                monitor(false),
                                                                only_csp(false),
                                                                monitor_csp(false),
                                                                visit_all(default_visit_all),
                                                                nb_point(default_nb_point),
                                                                eval_period(default_eval_period),
                                                                nb_sols(default_nb_sols),
                                                                min_acpt_diam(default_min_acpt_diam),
                                                                nb_sivia_iter(default_nb_sivia_iter),
                                                                nb_optim_iter(default_nb_optim_iter),
                                                                y_sol_radius(default_y_sol_radius),
                                                                reg_acpt_error(reg_acpt_error),
                                                                min_goal(x_sys.goal != NULL)

{
    if(!min_goal && xy_sys.goal !=NULL) {
        Array<const ExprNode> args(xy_sys.goal->nb_arg());
        Array<const ExprSymbol> var;
        for(int i = 0;i<xy_sys.goal->nb_arg();i++) {
            const ExprSymbol& a = ExprSymbol::new_(xy_sys.goal->arg(i).dim);
            var.add(a);
            args.set_ref(i,a);
        }
        minus_goal_y_at_x = new Function(var,-(*xy_sys.goal)(args));



        local_search = new UnconstrainedLocalSearch(*minus_goal_y_at_x,IntervalVector(1));
        lsolve.local_solver = local_search;
        loc_solve.local_solver_over_y = local_search;
        loc_solve.local_solver_over_x = new UnconstrainedLocalSearch(*xy_sys.goal,IntervalVector(1));


        Affine2Eval* aff_eval = new Affine2Eval(*(xy_sys.goal));
        lsolve.affine_goal = aff_eval;

        lsolve.goal_abs_prec = goal_rel_prec/100; // set goal prec of maximization problem lower than minimization

        }

        Array<const ExprNode> args_csp(max_fa_y_cst_sys.goal->nb_arg());
        Array<const ExprSymbol> var_csp;
        for(int i = 0;i<max_fa_y_cst_sys.goal->nb_arg();i++) {
            const ExprSymbol& a = ExprSymbol::new_(max_fa_y_cst_sys.goal->arg(i).dim);
            var_csp.add(a);
            args_csp.set_ref(i,a);
        }

        minus_goal_csp_y_at_x = new Function(var_csp,-(*max_fa_y_cst_sys.goal)(args_csp));

        local_search_csp = new UnconstrainedLocalSearch(*minus_goal_csp_y_at_x,IntervalVector(1));
        fa_lsolve.local_solver = local_search_csp;
        fa_loc_solve.local_solver_over_y = local_search_csp;
        fa_loc_solve.local_solver_over_x = new UnconstrainedLocalSearch(*(max_fa_y_cst_sys.goal),IntervalVector(1));


        Affine2Eval* aff_eval_csp = new Affine2Eval(*(max_fa_y_cst_sys.goal));
        fa_lsolve.affine_goal = aff_eval_csp;

        fa_lsolve.goal_abs_prec = 1e-2;

}


OptimMinMax::~OptimMinMax() {
    cout<<"Call destructor for optiminmax object: "<<this<<endl;
    cout<<"try to flush buffer"<<endl;
    if(buffer!=NULL){
	buffer->flush();
        cout<<"buffer flushed"<<endl;
	delete &(buffer->cost1());
	delete &(buffer->cost2());
        cout<<"buffer cost func deleted"<<endl;
	delete buffer;
    }
    cout<<"buffer deleted"<<endl;
        delete bsc;
    cout<<"bsc deleted"<<endl;
//    delete minus_goal_y_at_x;
//    delete minus_goal_csp_y_at_x;
//    cout<<"minus goal deleted"<<endl;
}


bool OptimMinMax::check_optimizer() {
    if(x_sys.goal != NULL && lsolve.xy_sys.goal != NULL) { // two objective functions-> error
        cout<<" Error: Two ojective functions found, choose either x depending function to minimize or x and y depending function for min max optim."<<endl;
        return false;
    }
    return true;
}


Optimizer::Status OptimMinMax::optimize() {
    bool problem_ok = check_optimizer();
    if(!problem_ok)
        return INFEASIBLE;
    return optimize(x_sys.box,POS_INFINITY);
}

void OptimMinMax::init_lsolve() {
    lsolve.local_search_iter = local_iter;
    lsolve.visit_all = visit_all;
}

void OptimMinMax::init_fa_lsolve() {
    fa_lsolve.local_search_iter = local_iter_csp;
    fa_lsolve.visit_all = visit_all_csp;
}

void OptimMinMax::init_loc_solve() {
    loc_solve.min_acpt_diam = min_acpt_diam;
    loc_solve.nb_optim_iter = nb_optim_iter;
    loc_solve.nb_sols = nb_sols;
    loc_solve.y_sol_radius = y_sol_radius;
    loc_solve.reg_acpt_error = reg_acpt_error;
    loc_solve.nb_sivia_iter = nb_sivia_iter;
}

void OptimMinMax::local_optimize(const IntervalVector& x_box_ini, double obj_init_bound) {


    cout<<" nb arg of fa cst: "<<fa_lsolve.xy_sys.goal->nb_arg()<<endl;
    cout<<" nb var of fa cst: "<<fa_lsolve.xy_sys.goal->nb_var()<<endl;
    Array<const ExprNode> args(fa_lsolve.xy_sys.goal->nb_arg());
    Array<const ExprSymbol> var;
    for(int i = 0;i<1;i++) {
        const ExprSymbol& a = ExprSymbol::new_(fa_lsolve.xy_sys.goal->arg(i).dim);
        var.add(a);
        args.set_ref(i,a);
    }

    cout<<" variable created"<<endl;

    Vector best_sol(x_box_init.size());


    loc_solve.nb_optim_iter = 0;
    loc_solve.nb_sivia_iter = 0;
    vector<Vector> max_points;
    vector<Matrix> max_sols;
    double loup = obj_init_bound;
    Vector x_ini = x_box_init.random();
    cout<<" xini random: "<<x_ini<<endl;

    Cell * root = new Cell(x_ini);
    lsolve.add_backtrackable(*root,y_box_init,critpr);
    buffer->cost1().add_backtrackable(*root);
    buffer->cost2().add_backtrackable(*root);

    lsolve.nb_iter = choose_nbiter(true,false,root);
    lsolve.prec_y = prec_y;
    lsolve.list_elem_max = 0;
    lsolve.visit_all = false;
    lsolve.optimize(root,loup); // eval maxf(midp,heap_copy), go to minimum prec on y to get a thin enclosure
    DataMinMax * data_x = &(root->get<DataMinMaxOpti>());
    cout<<" maximum reached for y = "<<(*(data_x->best_sol))<<endl;
    cout<<" fmax: "<<data_x->fmax<< " eval at best max point: "<<lsolve.xy_sys.goal->eval(*(data_x->best_sol))<<"  which is "<<*data_x->best_sol<<endl;
    max_points.push_back(data_x->best_sol->subvector(x_box_init.size(),x_box_init.size()+y_box_init.size()-1).mid());

    delete root;

    // also put some random points....
    cout<<"y_box init: "<<y_box_init<<endl;
    for(int i =0;i<2;i++) {
        max_points.push_back(2*y_box_init.random()-y_box_init.lb());
        cout<<" add random point "<<max_points.back()<<endl;
    }


//    if(fa_y_cst) {
//        Cell* fay_cell = new Cell(x_ini);
//        fa_lsolve.add_backtrackable(*fay_cell,y_box_init_fa,critpr);
//        cout<<"local fa ctc solve at x = "<<x_ini<<endl;
//        fa_lsolve.nb_iter = choose_nbiter(true,true);
//        fa_lsolve.prec_y = prec_y;
//        fa_lsolve.list_elem_max = 0;
//        fa_lsolve.visit_all = false;
//        fa_lsolve.optimize(fay_cell,POS_INFINITY);
//        DataMinMax * data_x_cst = &(root->get<DataMinMaxCsp>());
//        cout<<"  fa cst max: "<<data_x_cst->fmax<<endl;
//        if(data_x_cst->fmax.lb()>0)
//            best_sol = x_ini;
//        IntervalVector y_max_fa_cst = data_x_cst->best_sol->subvector(x_box_init.size(),x_box_init.size()+y_box_init.size()-1);
//        cout<<" try to set args1..."<<endl;
//        args[1] = y_max_fa_cst;
//        cout<<"done, try to create function..."<<endl;
//        Function* f_cst_at_y = new Function(var,(*fa_lsolve.xy_sys.goal)(args,y_max_fa_cst));
//        cout<<" try to eval new cst at y function: ";
//        cout<<f_cst_at_y->eval(x_ini)<<endl;
//        NumConstraint cst_at_y(*f_cst_at_y,LEQ);
//        x_sys.ctrs.add(cst_at_y); // add g(x,y)<0 in x_sys constraint, containing Cx typed constrains

//        delete fay_cell;

//    }


    //	max_points.push_back(Vector(1,3));
//	max_points.push_back(Vector(1,-3))
    bool max_not_found = true;


    for(int i = 0; i<100;i++) {
            loc_sols = pair<std::vector<Vector>,std::vector<Matrix> >(max_points,max_sols);
            cout<<" optimize with y max ";
            for(int j=0;j<max_points.size();j++) {
                    cout<<max_points.at(j)<<",  ";
            }
            cout<<endl;
            optimize(x_box_ini,obj_init_bound);
            cout<<" optimization done, reset loup "<<endl;
            loup = POS_INFINITY;
            Cell* x_cell = new Cell(loup_point);
            lsolve.add_backtrackable(*x_cell,y_box_init,critpr);
            cout<<"local solve at x = "<<loup_point<<endl;
            lsolve.nb_iter = choose_nbiter(true,false,x_cell);
            lsolve.prec_y = prec_y;
            lsolve.list_elem_max = 0;
            lsolve.visit_all = false;
            lsolve.optimize(x_cell,loup);
            data_x = &(x_cell->get<DataMinMaxOpti>());
            cout<<" get fmax "<<data_x->fmax<<endl;
            Vector y_max = data_x->best_sol->subvector(x_box_init.size(),x_box_init.size()+y_box_init.size()-1).mid();
            cout<<" max_y sol : "<< y_max<<" gives: "<<lsolve.xy_sys.goal->eval(*(data_x->best_sol))<<endl;
            max_points.push_back(y_max);
//		max_not_found &= loc_solve.check_twins(y_max,max_points);
    }
    return;
	
}



Optimizer::Status OptimMinMax::optimize(const IntervalVector& x_box_ini1, double obj_init_bound) {
    cout<<"start optimization"<<endl;
//	bool monitor = false;
	if (trace) lsolve.trace =trace;
	// we estimate that one iteration is at most 1% of the total time
//	lsolve.timeout=timeout/100;
//        fa_lsolve.timeout = timeout/100;

	loup=obj_init_bound;

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	x_box_init = x_box_ini1;

	//****** initialization of the first Cell********
	Cell * root = new Cell(x_box_init);
	bsc->add_backtrackable(*root);
    lsolve.add_backtrackable(*root,y_box_init,critpr);
    buffer->cost1().add_backtrackable(*root);
    buffer->cost2().add_backtrackable(*root);

    if(fa_y_cst) {
        fa_lsolve.add_backtrackable(*root,y_box_init_fa,critpr_csp);
    }
    buffer->critpr = heap_prob;

    //************* set light optim minmax solver param *********
    init_lsolve();

    //************* set csp light optim minmax solver param *********
    init_fa_lsolve();

    //************* set  light local solver param *********
    init_loc_solve();

	//****** x_heap initialization ********
	nb_cells=0;
        buffer->flush();

	// *** initialisation Algo ***
	loup_changed=false;
	double ymax;
    initial_loup=obj_init_bound;
	loup_point=x_box_init.mid();
    time=0;
    Timer::reset_time();
    Timer::start();

    if(!handle_cell(root))
        return INFEASIBLE;
//    spawn(root);

    update_uplo();

	//monitoring variables, used to track upper bound, lower bound, number of elem in y_heap and heap_save at each iteration
	std::vector<double> ub,lb,nbel,nbyel;
	long long int nbel_count(0);
        bool handle_res1,handle_res2;

        try {
                while (!buffer->empty() && (loup-uplo)>goal_rel_prec) {
			//			if (trace >= 2) cout << " buffer " << buffer << endl;
			if (trace >= 2) buffer->print(cout);


			loup_changed=false;

//            cout<<" get cell "<<endl;

            Cell *c = buffer->pop();
                        if(monitor) {
                            DataMinMax * data_x = &(c->get<DataMinMaxOpti>());
                            nbel_count -= data_x->y_heap->size();
                        }



//            if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;


            //====== only for purtpose of non inheritance comparison
            if(!propag){
                DataMinMax * data_x = &(c->get<DataMinMaxOpti>());

                data_x->y_heap->flush();
                Cell * y_cell = new Cell(y_box_init);
                y_cell->add<OptimData>();
                data_x->y_heap->push(y_cell);
                if(fa_y_cst) {
                    DataMinMax * data_x_csp = &(c->get<DataMinMaxCsp>());
                    data_x_csp->y_heap->flush();
                    Cell * y_cell_csp = new Cell(y_box_init_fa);
                    y_cell_csp->add<OptimData>();
                    data_x_csp->y_heap->push(y_cell_csp);
                }
            }


//            cout<<"******************************"<<endl;
//            cout<<"initial cell, box: "<<c->box<<endl<<"         fmax: "<<data_x->fmax<<endl;
//            cout<<"        y_heap size: "<<data_x->y_heap->size()<<endl;
            try {
//                cout<<"bisect cell "<<endl;
                pair<Cell*,Cell*> new_cells=bsc->bisect_cell(*c);

//                DataMinMax * data_x1 = &((c)->get<DataMinMaxCsp>());
//                Cell *y_cell = data_x1->y_heap->pop1();
//                fa_lsolve.check_already_in(y_cell,data_x1->y_heap);

//                show_yheap(c);
//                cout<<"delete cell "<<endl;
                delete c; // deletes the cell.
//                cout<<"deleted"<<endl;

//                DataMinMax * data_x1 = &((new_cells.first)->get<DataMinMaxOpti>());
//                Interval fmax_ini = data_x1->fmax;
//                cout<<" ===== init first cell box: "<<new_cells.first->box<<" of fmax: "<<data_x1->fmax<<endl;
//                cout<<"handle first cell"<<endl;
                handle_res1 = handle_cell(new_cells.first);

//                if(handle_res1) show_yheap(new_cells.first);
//                if(handle_res1) {
//                    if(!data_x1->fmax.is_subset(fmax_ini)) {
//                        cout<<" CRITICAL ERROR: Bound computation error"<<endl;
//                        return INFEASIBLE;
//                    }
//                    cout<<" ===== first cell box after handel cell: "<<new_cells.first->box<<" of fmax: "<<data_x1->fmax<<endl;

//                }

//                else
//                    cout<<" ==== first cell deleted"<<endl;

                if(handle_res1 && monitor) {
//					cout<<"try to get backtrack 1"<<endl;
                    DataMinMax * data_x1 = &((new_cells.first)->get<DataMinMaxOpti>());
//                    cout<<"get backtrack 1"<<endl;
                    nbel_count += data_x1->y_heap->size();
                }


//                DataMinMax * data_x2 = &((new_cells.second)->get<DataMinMaxOpti>());
//                cout<<" ===== init second cell box: "<<new_cells.second->box<<" of fmax: "<<data_x2->fmax<<endl;
                handle_res2 = handle_cell(new_cells.second);
//                if(handle_res2) show_yheap(new_cells.second);
//                if(handle_res2) {
//                    if(!data_x2->fmax.is_subset(fmax_ini)) {
//                        cout<<" CRITICAL ERROR: Bound computation error"<<endl;
//                        return INFEASIBLE;
//                    }
//                    cout<<" ===== second cell box after handel cell: "<<new_cells.second->box<<" of fmax: "<<data_x2->fmax<<endl;
//                }
//                    else
//                      cout<<" ==== second cell deleted"<<endl;

                if(handle_res2 && monitor) {
//                    cout<<"try to get backtrack 2"<<endl;
                    DataMinMax * data_x2 = &((new_cells.second)->get<DataMinMaxOpti>());
//                    cout<<"get backtrack 2"<<endl;
                    nbel_count += data_x2->y_heap->size();
                }

//                if(handle_res2) {
//                    bool spwn = spawn(new_cells.second);
////                    cout<<" spawn cell res nd: "<<spwn<<endl;

//                }
//                cout<<"second cell ok"<<endl;


                if (uplo_of_epsboxes == NEG_INFINITY) {
                        cout << " possible infinite minimum " << endl;
                        break;
                }
                if (loup_changed) {
                        // In case of a new upper bound (loup_changed == true), all the boxes
                        // with a lower bound greater than (loup - goal_prec) are removed and deleted.
                        // Note: if contraction was before bisection, we could have the problem
                        // that the current cell is removed by contractHeap. See comments in
                        // older version of the code (before revision 284).

                        ymax=compute_ymax();

                        buffer->contract(ymax);
                        //cout << " now buffer is contracted and min=" << buffer.minimum() << endl;


                        if (ymax <= NEG_INFINITY) {
                                if (trace) cout << " infinite value for the minimum " << endl;
                                break;
                        }
                        if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
                }
                if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
                update_uplo();


                if(monitor)
                {

                        lb.push_back(uplo);
                        ub.push_back(loup);
                        nbel.push_back(buffer->size());
                        nbyel.push_back(nbel_count);
                }


                Timer::check(timeout);

            }
            catch (NoBisectableVariableException& ) {
                    bool res=handle_cell(c);
                    if (res) update_uplo_of_epsboxes(c->get<DataMinMaxOpti>().fmax.lb());
                    //if (trace>=1) cout << "epsilon-box found: uplo cannot exceed " << uplo_of_epsboxes << endl;
                    update_uplo(); // the heap has changed -> recalculate the uplo

            }
        }
        if(monitor)
            export_monitor(&ub,&lb,&nbel,&nbyel);
    }
    catch (TimeOutException& ) {
            Timer::stop();
            time = Timer::get_time();
            return TIME_OUT;
    }

    Timer::stop();
    time = Timer::get_time();

    if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && goal_abs_prec==0 && goal_rel_prec==0)))
            return INFEASIBLE;
    else if (loup==initial_loup)
            return NO_FEASIBLE_FOUND;
    else if (uplo_of_epsboxes == NEG_INFINITY)
            return UNBOUNDED_OBJ;
    else
            return SUCCESS;
}





bool  OptimMinMax::handle_cell(Cell * x_cell) {
	bool local_eval = (!loc_sols.first.empty() || !loc_sols.second.empty());

        DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());

        //clear points of fsbl_point_list inherited from mother box that do not belong to this box.
//        cout<<" clear point "<<endl;
//       data_x->clear_notin_point(x_cell->box,true);

        ofstream out;
        if(monitor_csp) {
            out.open("paver.txt",std::ios_base::app);
            if(!out.is_open())
                cout<<"ERROR: cannot open paver.txt"<<endl;
        }
	// ATTENTION CECI EST FAUX !!!!!
	//double ymax;
	//if (loup==POS_INFINITY) ymax=POS_INFINITY;
	//else ymax= compute_ymax()+1.e-15;
	// update the upper bound with the loup
	//data_x->fmax &= Interval(NEG_INFINITY,loup);
	//std::cout <<" debut "<<data_x->fmax <<std::endl;
	//if (data_x->fmax.is_empty()) {
	//	delete x_cell;
	//	return false;
	//}
	// ATTENTION CECI EST FAUX !!!!!

	//***************** contraction w.r.t constraint on x ****************
//        cout<<"checking ctr.... "<<endl;
//        cout<<"current box: "<<x_cell->box <<" fmax: "<<data_x->fmax<<endl;
        IntervalVector tmpbox(x_cell->box);
        int res_cst = check_constraints(x_cell,false);
//        cout<<" constraint res: "<<res_cst<<endl;
        if (res_cst == 0) {
            if(monitor_csp) {
                for(int i=0; i<tmpbox.size();i++) {
                    out<<tmpbox[i].lb()<<" "<<tmpbox[i].ub()<<" ";
                }
            out<<res_cst<<endl;
            }
            out.close();
            return false;
        }
        else if (res_cst==2 && only_csp) {
            if(monitor_csp) {
                for(int i=0; i<x_cell->box.size();i++)
                    out<<x_cell->box[i].lb()<<" "<<x_cell->box[i].ub()<<" ";
                out<<res_cst<<endl;
            }

            data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
            delete x_cell; // need to delete x_cell because not deleted in check cst.
            out.close();
            return false;
        }
        else if (data_x->pu != 1)     {
                x_ctc.contract(x_cell->box);
                if(x_cell->box.is_empty()) {
                        //vol_rejected += x_cell->box.volume();
//                        data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                        delete x_cell;
                        //cout<<"loup : "<<loup<<" get for point: x = "<<best_sol<<" y = "<<max_y<<" uplo: "<<uplo<< " volume rejected: "<<vol_rejected/init_vol*100<<endl;
                        return false;
                }
                /*  TODO faire des options de TRACE
            else if(x_cell->box !=  box_mem) {
                vol_rejected += box_mem.volume()-x_cell->box.volume();
                cout<<"loup : "<<loup<<" get for point: x = "<<best_sol<<" y = "<<max_y<<" uplo: "<<uplo<< " volume rejected: "<<vol_rejected/init_vol*100<<endl;
            }
                 */
        }
        //clear points of fsbl_point_list inherited from mother box that do not belong to this box.

//        cout<<"done"<<endl;
//        cout<<"midpoint evaluation... "<<endl;

        if(!only_csp) {
//            cout<<"size of sol struct: "<<data_x->fsbl_pt_list.size()<<endl;

            //************* point evaluation ****************
            eval_counter++;
            if(eval_counter%eval_period == 0) {
                for(int i=0;i<nb_point;i++)
                {
                    Cell *x_copy = new Cell(*x_cell); // copy of the Cell and the y_heap
                    bool found_feas_pt = get_feasible_point(x_copy);
    //                cout<<"constraints fa y checked"<<endl;
                    if(found_feas_pt)    { // we found a feasible point

                        //                cout<<endl<<"**************"<<endl;
    //                    cout<<"try feasiblem point "<<x_copy->box<<", init fmax: "<<x_copy->get<DataMinMaxOpti>().fmax<<endl;
                        lsolve.nb_iter = choose_nbiter(true,false,x_cell);
                        lsolve.prec_y = prec_y;
                        lsolve.list_elem_max = 0; // no limit on heap size
                        lsolve.visit_all = false; // no need to visit all leaves in midpoint
    //                    cout<<"run point eval at "<<x_copy->box<<endl;
                        bool res1;
                        if(!local_eval) {
                            if(min_goal) {
                                res1 = eval_goal(x_copy,loup);
                            }
                            else {
                                lsolve.nb_iter = choose_nbiter(true,false,x_cell);
                                lsolve.prec_y = prec_y;
                                lsolve.list_elem_max = 0; // no limit on heap size
                                lsolve.visit_all = false; // no need to visit all leaves in midpoint
                                cout<<"run point eval at "<<x_copy->box<<endl;
                                res1 = lsolve.optimize(x_copy,loup); // eval maxf(midp,heap_copy), go to minimum prec on y to get a thin enclosure
                                DataMinMaxOpti * data_x_copy = &(x_copy->get<DataMinMaxOpti>());
                                cout<<" mid point result: "<<data_x_copy->fmax<<endl;

                            }
                        }
                        else {
    //						cout<<"  eval midpoint"<<endl;
                                Interval eval = loc_solve.eval_backward_max_solutions(loc_sols,x_copy->box,loup);
    //						cout<<" res eval midpoint : "<<eval<<endl;
                                if(eval.is_empty())
                                        res1 = false;
                                else {
                                        res1 = true;
                                        DataMinMaxOpti * data_x_copy = &(x_copy->get<DataMinMaxOpti>());
                                        data_x_copy->fmax = eval;
                                }
    //						cout<<" res1 = "<<res1<<endl;
                        }
                        lsolve.visit_all = visit_all; // reset visit all to initial value
                        if(!res1) {
                                delete x_copy;
                        }
                        else {
                            Interval ev = x_copy->get<DataMinMaxOpti>().fmax;
                            IntervalVector ysol = x_copy->get<DataMinMaxOpti>().y_heap->top()->box;
    //                        cout<<" point eval at: "<<x_copy->box<<endl<<"         res: "<<ev<<endl<<"         at: "<<ysol<<endl<<"         at: "<<lsolve.best_point_eval<<endl;
                            Vector sol = lsolve.best_point_eval.mid();

    //                        feasible_point * sol_struct = new feasible_point(sol,x_copy->get<DataMinMaxOpti>().fmax);
    //                        data_x->fsbl_pt_list.push_back(sol_struct);
    //                        cout<<"size of sol struct: "<<data_x->fsbl_pt_list.size()<<endl;

                            //                    cout<<"res midp eval: "<<x_copy->get<DataMinMaxOpti>().fmax<<endl;
                            double new_loup = x_copy->get<DataMinMaxOpti>().fmax.ub();

                            //data_x->fmax &= Interval(NEG_INFINITY,new_loup);

                            if(new_loup<loup) { // update best current solution
                                //cout<<"worst case over Y: "<<x_copy->get<DataMinMax>().y_heap->top1()->box<<endl;
                                loup = new_loup;
                                loup_changed = true;
                                loup_point = (x_copy->box.mid());

                                if (trace) cout << "[mid]"  << " loup update " << loup  << " loup point  " << loup_point << endl;
                                //cout<<"loup : "<<loup<<" get for point: x = "<<best_sol<<" y = "<<max_y<<" uplo: "<<uplo<< " volume rejected: "<<vol_rejected/init_vol*100<<endl;
                            }
                            delete x_copy; // delete copy of the heap, no more use and it was not delete in light optim since res1 = 1
                        }
                    }
                }
            }
//            cout<<"end midpoint eval"<<endl;
//            cout<<" evaluation of objectif... "<<endl;

            //************ evaluation of f(x,y_heap) *****************
            lsolve.prec_y = compute_min_prec(x_cell->box,false);
            lsolve.nb_iter = choose_nbiter(false,false,x_cell);
            lsolve.list_elem_max = compute_heap_max_size(data_x->y_heap->size(),false);
    //        CellCostmaxFmaxub cost;
    //        cout<<"cost after eval: "<<cost.cost(*x_cell)<<endl;
    //        DataMinMax *data_opt_x = &(x_cell->get<DataMinMax>());
    //        cout<<"fmax from DataMinMax cast: "<<data_opt_x->fmax<<endl;
            // compute
//                    cout<<"run light optim with box: "<<x_cell->box<<endl;
//                    cout<<"            init fmax: "<<data_x->fmax<<endl;
//                    cout<<"            init x_box: "<<x_cell->box<<endl;

//            cout<<"current x box: "<<x_cell->box<<endl;
//            bool res(true);
            bool res;
            if(!local_eval) { // regular algo
                if(min_goal) { // not a minmax problem, only a min problem-> suffices
                    res = eval_goal(x_cell,loup);

                }

                else {
//                    cout<<" befor lightsolver, fmax: "<<data_x->fmax<<endl;
//                    cout<<" nb element in yheap: "<<data_x->y_heap->size()<<endl;
                    res =lsolve.optimize(x_cell,loup);
//                    cout<<" after lightsolver, fmax: "<<data_x->fmax<<endl;
//                    cout<<" max get for point "<<*data_x->best_sol<<" which gives "<<lsolve.xy_sys.goal->eval(*data_x->best_sol)<<endl;
                    double how_far = (data_x->fmax.lb()-uplo)/uplo;
        //            cout<<"  uplo: "<<uplo<<  "how far befor calling lightloc: "<<how_far<<endl;

        //            cout<<" entering locsolve, fmax: "<<data_x->fmax<<endl;
                    how_far = 2;
                    if (res && (nb_optim_iter!=0||nb_sivia_iter!=0) && how_far>1 )
                            res = loc_solve.compute_supy_lb(x_cell,uplo,loup,minus_goal_y_at_x);

                }
            }
            else { // mitsos derived algo, evaluate function at a finite number of point or regression y_x^* = Ax contained in loc_solve.
                Interval eval = loc_solve.eval_backward_max_solutions(loc_sols,x_cell->box,loup);
                if(eval.is_empty())
                        res = false;
                else {
                        res = true;
                        data_x->fmax &= Interval(eval.lb(),POS_INFINITY);
                }
            }

//            cout<<" done "<<endl;


//            cout<<"res : "<<res<<endl;
//            cout<<" respects constraint? "<<data_csp_x->fmax<<endl;
//			cout<<"fmax after lightoptimMinMax: "<<data_x->fmax<<endl;
//            if(res) {
//                OptimData* ydata = &(data_x->y_heap->top1()->get<OptimData>());
//                cout<< "top 1 pf: "<<ydata->pf<<endl;
//            }
//            cout<<"             x_box after lightoptimMinMax: "<<x_cell->box<<endl;
//            cout<<"out of locsolve, fmax:"<<data_x->fmax<<endl;
//            cout<<"             x_box after lightoptimMinMax: "<<x_cell->box<<endl;
//            if (res)
//                res &= loc_solve.compute_supy_lb(x_cell,uplo,loup,minus_goal_y_at_x);
//            cout<<"lightloc res: "<<res<<endl;
//                cout<<"        res : "<<res<<endl;
//            cout<<"            fmax after lightLocalOptim: "<<data_x->fmax<<endl;
//            cout<<"            xbox after lightLocalOptim: "<<x_cell->box<<endl;
//            cout<<"optimization on box finished"<<endl;
//            if(res) cout<<"          light optim result: "<<data_x->fmax<< endl;
            //std::cout<<" apres res="<<res<<" bound= "<<data_x->fmax <<std::endl;

    //        cout<<"done"<<endl;

            if(!res) { // certified that x box does not contains the solution
                    // TODO Trace
                    //vol_rejected += x_cell->box.volume();
                    //x_cell->y_heap.flush();
                    //                delete x_cell;
                    //cout<<"loup : "<<loup<<" get for point: x = "<<best_sol<<" y = "<<max_y<<" uplo: "<<uplo<< " volume rejected: "<<vol_rejected/init_vol*100<<endl;
//                cout<<" delete x_cell"<<endl;
//                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
//                cout<<" fsbl_list cleared"<<endl;
                delete x_cell;
//                cout<<"x_cell deleted"<<endl;
                return false;
            }
//            cout<<"res checked"<<endl;
    //        else // need to delete x_copy
    //            delete x_copy;
    //        cout<<"done"<<endl;
    //        cout<<"check min prec..."<<endl;
        }

        //***** if x_cell is too small ******************
        if(x_cell->box.max_diam()<prec) {
            cout<<"Min prec, box: "<<x_cell->box<<endl;
            int ctr_ok = 2;
            tmpbox = x_cell->box;
//            if(fa_y_cst) { // go min precision on csp problem, else suppose the ctr on x is respect because x_cell not deleted when first check constraints at the begining of handle_cell()
                ctr_ok= check_constraints(x_cell, true);
//            }
            if(only_csp && monitor_csp) {
                for(int i=0; i<x_cell->box.size();i++) {
                    out<<tmpbox[i].lb()<<" "<<tmpbox[i].ub()<<" ";
                }
                out<<ctr_ok<<endl;
                if(ctr_ok!=0) {
                    data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                    delete x_cell;
                }
                out.close();
                return false;
            }
            if (ctr_ok !=0 && !only_csp) {
                lsolve.nb_iter = choose_nbiter(true,false,x_cell);   // need to be great enough so the minimum precision on y is reached
                lsolve.prec_y = prec_y;
                lsolve.list_elem_max = 0; // no limit on heap size
//                cout<<"fmax ini: "<<data_x->fmax<<endl;
                //cout<<"for box: "<<x_cell->box<<endl;
                bool res;
                if(min_goal)
                    res = eval_goal(x_cell,loup);
                else
                    res =lsolve.optimize(x_cell,loup); // eval maxf(midp,heap_copy), go to minimum prec on y to get a thin enclosure
//                cout<<"fmax min prec: "<<data_x->fmax<<endl;
                if(!res) {
                    data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                    delete x_cell;
                }
                else{
                        update_uplo_of_epsboxes(data_x->fmax.lb());
                        data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                        delete x_cell;
                }
            }
                //cout<<"minprec reached! "<<" box: "<<x_cell->box<<" eval full prec: "<<x_cell->fmax <<endl;
                //cout<<"loup : "<<loup<<" get for point: x = "<<best_sol<<" y = "<<max_y<<" uplo: "<<uplo<< " volume rejected: "<<vol_rejected/init_vol*100<<endl;
                //min_prec_reached = true;
                return false;
        }
//        cout<<"done"<<endl;
//        cout<<"update cell data..."<<endl;


        /*	if (data_x->fmax.is_empty()) {
                delete x_cell;
                return false;
        }*/

        // update optim data of the cell
        buffer->cost1().set_optim_data(*x_cell,x_sys);
        buffer->cost2().set_optim_data(*x_cell,x_sys);
        buffer->push(x_cell);
        nb_cells++;
//        std::cout<<"      final value "<<data_x->fmax <<std::endl;
//        cout<<"done, exit handle_cell()"<<endl;
	return true;
}

bool OptimMinMax::eval_goal(Cell* x_cell,double loup) {
//    cout<<"in eval goal"<<endl;
    x_sys.goal->backward(Interval(NEG_INFINITY,loup),x_cell->box);
//    cout<<" backward contraction: "<<x_cell->box<<endl;
    if(x_cell->box.is_empty())
        return false;
    DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());
    data_x->fmax = x_sys.goal->eval(x_cell->box);
//    cout<<"box: "<<x_cell->box<<"   evaluation: "<<data_x->fmax<<endl;
    if(data_x->fmax.lb()>loup)
        return false;
    return true;

}



double OptimMinMax::compute_min_prec( const IntervalVector& x_box,bool csp) {

    if(!csp) {
        if(min_prec_coef == 0)
                return prec_y;

//        double ratio = 0;
//        for(int r=0;r<x_box_init.size();r++)
//                ratio += (x_box_init[r]).diam()/(x_box[r]).diam();
//        ratio/=x_box.size();
//        double new_prec = y_box_init.max_diam()*ratio/min_prec_coef;
//        double prec_ratio = prec/prec_y;
//        double new_prec = x_box.max_diam()*prec_ratio*min_prec_coef;
        double new_prec = x_box.max_diam()*min_prec_coef;
        return new_prec>prec_y?new_prec:prec_y;
    }
    else {
        if(min_prec_coef_csp == 0)
                return prec_y;

//        double ratio = 0;
//        for(int r=0;r<x_box_init.size();r++)
//                ratio += (x_box_init[r]).diam()/(x_box[r]).diam();
//        ratio/=x_box.size();
//        double new_prec = y_box_init_fa.max_diam()*ratio/min_prec_coef_csp;

        double new_prec = x_box.max_diam()*min_prec_coef;
        return new_prec>prec_fa_y?new_prec:prec_fa_y;
    }
}

int OptimMinMax::choose_nbiter(bool midpoint_eval,bool csp,Cell* x_cell) {

    if(!midpoint_eval) {
        if(!csp) {
            DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());
            data_x->nb_bisect+=iter;
            if(propag)
                return iter;
            else
                return data_x->nb_bisect;
        }
        else {
            DataMinMaxCsp * data_x = &(x_cell->get<DataMinMaxCsp>());
            data_x->nb_bisect+=iter_csp;
//            cout<<"nb bisect: "<<data_x->nb_bisect;
            if(propag)
                return iter_csp;
            else
                return data_x->nb_bisect;

        }
    }
        else
            return 0; // go to min prec
}

int OptimMinMax::compute_heap_max_size(int y_heap_size,bool csp) {
    if(!csp) {
        if(list_rate == 0 || y_heap_size>= list_elem_absolute_max) { // no constraint on list growth rate or max size already reached
                return list_elem_absolute_max;
        }
        return  y_heap_size+list_rate;
    }
    else {
        if(list_rate_csp == 0 || y_heap_size>= list_elem_absolute_max_csp) { // no constraint on list growth rate or max size already reached
                return list_elem_absolute_max_csp;
        }
        return  y_heap_size+list_rate_csp;
    }
}

bool OptimMinMax::get_feasible_point(Cell * elem) {
        cout<<endl<<"     get_feasible_point"<<endl;
        elem->box = elem->box.random(); //
//        elem->box = elem->box.mid(); //get the box (x,mid(y))
        int res = check_constraints(elem,true);
        cout<<"done, res = "<<res<<endl;
        if(res == 2)
            return true;
        return false;
}

int OptimMinMax::check_regular_ctr(const IntervalVector& box) {
        int res =2;
        for(int i=0;i<x_sys.nb_ctr;i++) {
                Interval int_res = x_sys.ctrs[i].f.eval(box);
                if(int_res.lb()>=0)
                        return 0;
                else if(int_res.ub()>=0)
                        res = 1;
        }
        return res;
}

int OptimMinMax::check_fa_ctr(Cell* x_cell,bool midp) {
    DataMinMax * data_csp = &(x_cell->get<DataMinMaxCsp>());
//    cout<<"data_csp->pu: "<<data_csp->pu<<endl;
    if(data_csp->pu != 1) {
        fa_lsolve.nb_iter = choose_nbiter(midp,true,x_cell);
        if(midp) {
            fa_lsolve.visit_all = false; // no need to visit the heap in midp
            fa_lsolve.list_elem_max = 0;
            fa_lsolve.prec_y = prec_fa_y;
        }
        else {
            fa_lsolve.list_elem_max = compute_heap_max_size(data_csp->y_heap->size(),true);
            fa_lsolve.prec_y = compute_min_prec(x_cell->box,true);
        }
        bool ok= fa_lsolve.optimize(x_cell,0);
        if(midp)
            cout<<" res of csp for all maximization eval: "<<data_csp->fmax<<endl;
        fa_lsolve.visit_all = visit_all_csp; // reset visit all to initial value
        if(!ok) {
//            cout<<"     fa ctr not satisfied for box: "<<x_cell->box<<endl;
            return 0;
        }
        else {
            if(data_csp->y_heap->empty()) {
                data_csp->pu = 1;
                cout<<" sic validated from empty list"<<endl;
                return 2;
            }
            else if(data_csp->y_heap->top1()->get<OptimData>().pf.ub()<0) {
                data_csp->pu = 1;
                data_csp->y_heap->flush();
                return 2;
            }

            else return 1;
        }
    }
    return 2;
}

int OptimMinMax::check_constraints(Cell * x_cell,bool midp) {
//        cout<<"     check constraints for box: "<<x_cell->box<<endl;
        int res_rctr = 2;
        int res_factr = 2;
        DataMinMaxOpti * data_opt = &(x_cell->get<DataMinMaxOpti>());
//        cout<<"     check regular ctr..."<<endl;
        if(data_opt->pu != 1)
            res_rctr = check_regular_ctr(x_cell->box);
        if(res_rctr == 2)
            data_opt->pu = 1;
        else if(res_rctr == 0) {
            if(!midp)// do not delete solution list if not box to discard (if mid point evaluation
                data_opt->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
            delete x_cell;
            return 0;
        }
//        cout<<"     done, res_rctr = "<<res_rctr<<endl;
//        cout<<"     check fa ctr..."<<endl;
        if(fa_y_cst)
            res_factr = check_fa_ctr(x_cell,midp);
        if(res_factr==0)
            delete x_cell;
//        cout<<"     done, res_factr = "<<res_factr<<endl;
//        cout<<"res_factr: "<<res_factr<<endl;
        if(res_rctr == 2 &&  res_factr == 2) // all ctr satisfied
            return 2;
        else if(res_rctr == 0 || res_factr == 0)
            return 0;
        else
            return 1;
}

//bool OptimMinMax::spawn(Cell * x_cell) {
//    if(loc_solve.perf_ctc_coef<perf_thresh) {
//        while(!loc_solve.csp_queue.empty()) {
//            Cell* spawn_cell = new Cell(*x_cell);
//            spawn_cell->box = loc_solve.csp_queue.front().first;
//            DataMinMaxOpti * data_opt = &(spawn_cell->get<DataMinMaxOpti>());
//            data_opt->fmax &= Interval(loc_solve.csp_queue.front().second.lb(),POS_INFINITY);
//            bool res = true;
////            bool res =lsolve.optimize(spawn_cell,loup);
//            if(!res)
//                delete spawn_cell;
//            else {
//                buffer->cost2().set_optim_data(*spawn_cell,x_sys);
//                buffer->push(spawn_cell);
//            }
//            loc_solve.csp_queue.pop();
//        }
//        delete x_cell;
//        return true;
//    }
//    buffer->cost2().set_optim_data(*x_cell,x_sys);
//    buffer->push(x_cell);
//    return false;
//}


void export_monitor(std::vector<double> * ub,std::vector<double> * lb,std::vector<double> * nbxel,std::vector<double> * nbyel) {
        std::ofstream out;
        //    std::string output_name = "log.txt";
        out.open("log.txt");
        //    out<<"x box: "<<box<<std::endl;
        for(int i = 0;i<ub->size();i++) {
                out<<ub->at(i)<<" "<<lb->at(i)<<" "<<nbxel->at(i)<<" "<<nbyel->at(i)<<std::endl;
        }
        out.close();

}

void OptimMinMax::show_yheap(Cell * x_cell) {
    DataMinMaxCsp * data_opt = &(x_cell->get<DataMinMaxCsp>());
    std::vector<Cell *> y_cells;
    Cell * y_cell;
    OptimData  *data_y;

    cout<<"   y_heap of box  "<<x_cell->box<< "contains "<<data_opt->y_heap->size()<<" elements: "<<endl;
    while(!data_opt->y_heap->empty()) {
        y_cell = data_opt->y_heap->pop1();
        data_y = &(y_cell->get<OptimData>());
//        cout<<"     "<<y_cell->box<<"  eval: "<<data_y->pf<<endl;
        cout<<"  pointer: "<<y_cell<<endl;
        IntervalVector box(x_cell->box.size()+y_cell->box.size());
        box.put(0,x_cell->box);
        box.put(x_cell->box.size(),y_cell->box);
//        cout<<"  evaluation here at "<<box<<" : "<<this->fa_lsolve.xy_sys.goal->eval(box)<<endl;
        y_cells.push_back(y_cell);
    }
    while(!y_cells.empty()){
        y_cell = y_cells.back();
        fa_lsolve.check_already_in(y_cell,data_opt->y_heap);
        data_opt->y_heap->push(y_cell);
        y_cells.pop_back();
    }
}
}

