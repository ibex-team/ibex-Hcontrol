//
// Created by joris on 17/11/2021.
//

#include "ibex_EvalMax.h"
#include "ibex_Timer.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;
namespace ibex {

    const double EvalMax::default_timeout = 20;
    const double EvalMax::default_goal_abs_prec = 1e-2;

//EvalMax::EvalMax(Function &f, int nx, int ny) {}
    EvalMax::EvalMax(NormalizedSystem &xy_sys, int nx, int ny, Ctc &ctc_xy) :
            trace(false), timeout(default_timeout),
            list_elem_max(0), nb_iter(0), prec_y(0),
            monitor(false), local_search_iter(0), xy_sys(xy_sys), goal_abs_prec(default_goal_abs_prec), ctc_xy(ctc_xy),
            bsc(new LargestFirst()), found_point(false), time(0), best_point_eval(xy_sys.box) {
//        TODO check
    }

    EvalMax::~EvalMax() = default;

    Interval EvalMax::eval(IntervalVector &X) {
//        TODO
    }

    Interval EvalMax::eval(IntervalVector &X, BoxProperties &prop, double loup) {
//        delete_save_heap(); TODO

        found_point  = false;
        BxpMinMax *data_x;
        Cell *x_cell = new Cell(X);

        if (csp_actif)
            data_x = dynamic_cast<BxpMinMax *>(prop[BxpMinMaxCsp::id]);
        else
            data_x = dynamic_cast<BxpMinMax *>(prop[BxpMinMaxOpti::id]);

        //std::cout <<"    DEB "<<data_x->fmax <<std::endl;
//                cout<<endl<<"*************************"<<endl;

        DoubleHeap<Cell> *y_heap = data_x->y_heap; // current cell
//        cout<<"get y_heap of size: "<<y_heap->size()<<endl;
//        cout<<"initial fmax: "<<data_x->fmax<<endl;
//        cout<<"box: "<<x_cell->box<<endl;

//        // Define the TimeOut of to compute the bounds of x_cell
        Timer timer;
        time = timer.get_time();
//
//        // ********** contract x_box with ctc_xy***************
//        //        if (!csp_actif) // no constraint if dealing with fa cst, better not use the useless contractor
//        //        {
        IntervalVector xy_box = xy_box_hull(x_cell->box); // TODO
        ctc_xy.contract(xy_box);

        if (xy_box.is_empty()) {
//                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
//                delete x_cell;
            return false;
        } else {
            // contract the result on x
            x_cell->box &= xy_box.subvector(0,x_cell->box.size()-1);
        }
        if (local_search_iter>0)
        {
            double lower_bound_ls = local_search_process(x_cell->box, xy_box, loup); // TODO
            if (lower_bound_ls > data_x->fmax.ub()) {
                ibex_error("ibex_LightOptimMinMax: error, lb >ub from local search. Please report this bug.");
            }
//            cout<<"local optim return new lb: "<< lower_bound_ls<<endl;
            if (lower_bound_ls > loup) {
//                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
//                delete x_cell;
                return false;
            }
            else {
                data_x->fmax &= Interval(lower_bound_ls, POS_INFINITY); // update lower bound on fmax
            }
        }


        //monitoring variables, used to track upper bound, lower bound, number of elem in y_heap and heap_save at each iteration
        std::vector<double> ub, lb, nbel, nbel_save;

        save_heap_ub = NEG_INFINITY;
//        cout<<"initial fmax: "<<data_x->fmax<<endl;



        // *********** loop ********************
        try {
            int current_iter = 1;
            while(!stop_crit_reached(current_iter, y_heap,data_x->fmax)) { // TODO
                if (trace >= 3) cout << *y_heap << endl;


                found_point  = false;

                Cell * y_cell = y_heap->pop(); // we extract an element with critprob probability to take it according to the first crit
                current_iter++;
//                                                std::cout<<"current_iter: "<<current_iter<<std::endl;
                if ((list_elem_max != 0 && ((y_heap->size() + heap_save.size())>list_elem_max)) || (y_cell->box.size())<prec_y) { // continue to evaluate cells of y_heap without increasing size of list, need to do it else nothing happend if list already reached max size
                    bool res = handle_cell(x_cell, y_cell, loup); // TODO
                    if (!res) { // x_cell has been deleted
                        return false;
                    }
                }
                else
                {
                    try {
                        pair<Cell*,Cell*> subcells_pair = bsc->bisect(*y_cell);// bisect tmp_cell into 2 subcells
                        delete y_cell;
                        //                                                                std::cout<<"handle first cell in light optim"<<std::endl;
                        bool res = handle_cell(x_cell, subcells_pair.first, loup);
                        //                                                                std::cout<<"first cell handled"<<std::endl;
                        if (!res) { // x_cell has been deleted
                            delete subcells_pair.second;
                            //                                                                                std::cout <<"       OUT 1 "<<std::endl;
                            return false;
                        }
                        //                                                                std::cout<<"handle second cell in light optim"<<std::endl;
                        res = handle_cell( x_cell, subcells_pair.second,loup);
                        //                                                                std::cout<<"second cell handled"<<std::endl;
                        if (!res) { // x_cell has been deleted
                            //                                                                                std::cout <<"       OUT 2 "<<std::endl;
                            return false;
                        }

//                                if (found_point && !csp_actif) { // a feasible solution has been found
//                                    y_heap->contract(-(data_x->fmax.lb())); // to check

//                                }
                        // ** contract not yet
                        //if (found_point) { // a feasible solution has been found
                        //	y_heap->contract(-(data_x->fmax.lb())); // to check
                        //}
                    }
                    catch (NoBisectableVariableException& ) {
                        bool res = handle_cell(x_cell, y_cell,loup);
                        cout << " no bisectable caught" << endl;

//                                if (res) heap_save.push_back(y_cell);
                        if (!res) return false;

                    }
                }
                if (monitor)
                {
                    if (!y_heap->empty()) {
                        lb.push_back(data_x->fmax.lb());
                        auto top1_minmax = dynamic_cast<BxpMinMax *>(y_heap->top1()->prop[BxpMinMax::id]);
                        if (!top1_minmax)
                            ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMax.");
                        if (save_heap_ub < top1_minmax->pf.ub())
                            ub.push_back(top1_minmax->pf.ub());
                        else
                            ub.push_back(save_heap_ub);
                        nbel.push_back(y_heap->size());
                        nbel_save.push_back(heap_save.size());
                    }
                }
//                                                std::cout<<"nb iter: "<<current_iter<<std::endl;
                timer.check(time+timeout);

            }

        }
        catch (TimeOutException& ) { }

//        cout<<"visit all, y_heap size: "<<y_heap->size()<<endl;
        //force visit of all elements of the heap
//        cout<<"visit all = "<<visit_all<<endl;
        if (visit_all) {
//            cout<<"             y_heap detail: "<<endl;
            while (!y_heap->empty()) {
                Cell * y_cell = y_heap->pop();
//                cout<<"                 "<<y_cell->box<<endl;
                bool res = handle_cell(x_cell, y_cell, loup, true);
                if (!res) return false;
                if (monitor)
                {
                    if (!y_heap->empty()) {
                        lb.push_back(data_x->fmax.lb());
                        auto top1_minmax = dynamic_cast<BxpMinMax *>(y_heap->top1()->prop[BxpMinMax::id]);
                        if (!top1_minmax)
                            ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMax.");
                        if (save_heap_ub < top1_minmax->pf.ub())
                            ub.push_back(top1_minmax->pf.ub());
                        else
                            ub.push_back(save_heap_ub);
                        nbel.push_back(y_heap->size());
                        nbel_save.push_back(heap_save.size());
                    }
                }
            }
        }
//        cout<<"all eval done"<<endl;
//        cout<<"out of loop"<<endl;
//                std::cout<<"light optim: out of loop"<<std::endl;
        /*
??    if (is_midp && !midp_hit) // midpoint x eval case: if no y found such as xy constraint respected, cannot keep the result
??        return Interval::EMPTY_SET;
         */
        // insert the y_box with a too small diameter (those which were stored in heap_save)
        //         std::cout<<"try to fill y_heap"<<std::endl;

//        if (visit_all == true) {
//            for(int i = 0; i<heap_save.size();i++) {
//                cout<<"         lightoptimMinMax: box: "<<heap_save.at(i)->box<<endl<<"                eval: "<<heap_save.at(i)->get<OptimData>().pf<<endl;
////                   <<"  mid eval: "<<xy_sys.goal->eval(heap_save.at(i)->box.mid())<<endl;
//            }
//        }
        fill_y_heap(*y_heap); // TODO
//        cout<<"heap filled"<<endl;

        //if (found_point) { // a feasible solution has been found
//        if (!csp_actif)
//            y_heap->contract(-(data_x->fmax.lb())); // to check
//            cout<<"contraction done"<<endl;
//        else
//            y_heap->contract(-loup); // in csp solve, remove all boxes y that satisfy the constraint for all
        //}
        //        std::cout<<"y_heap filled"<<std::endl;

        // ** contract y_heap now

        //        std::cout<<"found point pass"<<std::endl;


        if (y_heap->empty()) {
            cout <<"       OUT 3 "<< endl;
//            data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
//            delete x_cell;
//            if (csp_actif) // sic case: empty set y means that the set defined by the constraints g(x,y)<0 is empty, and therefore that the sic is respected over the empty set
//                return true;
//            else
//                return false; // minmax case: empty set means ????? suppose that x can be discarded????
            return csp_actif;
        }
//        cout<<"non empty heap"<<endl;

        // Update the lower and upper bound of of "max f(x,y_heap)"
        auto top1_minmax = dynamic_cast<BxpMinMax *>(y_heap->top1()->prop[BxpMinMax::id]);
        if (!top1_minmax)
            ibex_error("[ibex_EvalMax] -- null ptr in eval: top1 element has no BxpMinMax.");
        double new_fmax_ub = top1_minmax->pf.ub(); // get the upper bound of max f(x,y_heap)
        auto top2_minmax = dynamic_cast<BxpMinMax *>(y_heap->top2()->prop[BxpMinMax::id]);
        if (!top2_minmax)
            ibex_error("[ibex_EvalMax] -- null ptr in eval: top2 element has no BxpMinMax.");
        double new_fmax_lb = top2_minmax->pf.lb(); // get the lower bound of max f(x,y_heap)

        //std::cout<<"new_fmax_ub: "<<new_fmax_ub<<std::endl<<"new_fmax_lb: "<<new_fmax_lb<<std::endl<<"fmax_lb (from found point): "<<data_x->fmax.lb()<<std::endl;

        if (new_fmax_ub< new_fmax_lb) {
            ibex_error("ibex_LightOptimMinMax: error, please report this bug.");
        }

//                std::cout<<"fmax ini: "<<data_x->fmax<<std::endl;
        data_x->fmax &= Interval(new_fmax_lb, new_fmax_ub);
//                std::cout<<"new_fmax_lb: "<<new_fmax_lb<<" new_fmax_ub: "<<new_fmax_ub<<std::endl;
//                std::cout<<"       fmax final: "<<data_x->fmax<<std::endl;
        //        if (data_x->fmax.lb()>100000){
        //            std::cout<<"Issue: fmax_lb > 100,000"<<std::endl;

        //        }
                cout << "loup: " << loup << endl;

//        if (csp_actif) { // if SIC, all boxes lower than  can be deleted since of no interest (verify the constraint)
//            y_heap->contract(loup);
//            if (y_heap->empty())
//                cout<<" y_heap empty after loup contraction"<<endl;
//        }

        if (data_x->fmax.is_empty() || data_x->fmax.lb() > loup) {
//                                std::cout <<"       OUT 4 "<<std::endl;
//                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
//                delete x_cell;
            return false;
        }


        if (monitor)
        {
            lb.push_back(data_x->fmax.lb());
            ub.push_back(data_x->fmax.ub());
            nbel.push_back(y_heap->size());
            nbel_save.push_back(heap_save.size());
            export_monitor(&ub,&lb,&nbel,&nbel_save,x_cell->box); // TODO
        }
        best_point_eval = xy_box.mid();
        for(int i=0;i<xy_sys.box.size()-x_cell->box.size();i++) {
            best_point_eval[x_cell->box.size()+i] = y_heap->top1()->box[i].mid();
        }

//                std::cout <<"    FIN "<<data_x->fmax <<std::endl;
        return true;
    }


} // end namespace ibex