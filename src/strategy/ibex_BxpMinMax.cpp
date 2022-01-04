//
// Created by joris on 17/11/2021.
//

#include "ibex_BxpMinMax.h"

namespace ibex {


feasible_point::feasible_point(const Vector& box,const Interval& eval) : point(box), eval(eval) { }
feasible_point::feasible_point(const feasible_point& pt) : point(pt.point), eval(pt.eval) { }


Map<long,long,false>& BxpMinMax::ids() {
	static Map<long,long,false> _ids;
	return _ids;
}

Map<long,long,false>& BxpMinMaxSub::ids() {
	static Map<long,long,false> _ids;
	return _ids;
}


BxpMinMax::BxpMinMax(EvalMax& evalmax, int crit_heap) : Bxp(BxpMinMax::get_id(evalmax)),
		best_sol(NULL),
		y_heap(evalmax,crit_heap),
		nb_bisect(0),
		pu(0),
		evalmax(evalmax) {

}

BxpMinMax::~BxpMinMax() {
	y_heap.flush();
	delete best_sol;
}


void BxpMinMax::clear_fsbl_list() {
	fsbl_pt_list.clear();
}

void BxpMinMax::clear_notin_point(const IntervalVector &x_box, bool strong_del) {
	int size = fsbl_pt_list.size();
	std::vector<feasible_point> save_vect;
	//    std::cout<<"======================= "<<std::endl;
	//    std::cout<<" check in box: "<<x_box<<std::endl;
	for(int i=0;i<size;i++) {
		feasible_point pt = fsbl_pt_list.back();
		fsbl_pt_list.pop_back();
		//        std::cout<<" feas pt "<<pt->point<<std::endl;
		if((x_box.contains(pt.point.subvector(0,x_box.size()-1)))) {
			//            if(strong_del)
			//                delete pt;
			save_vect.push_back(pt);
			//            std::cout<<"      is deleted"<<std::endl;
		}
		//        else {
		//            std::cout<<"      is NOT deleted"<<std::endl;

		//        }
	}
	fsbl_pt_list = save_vect;
	//    std::cout<<"***********************"<<std::endl;
}

BxpMinMax::BxpMinMax(const BxpMinMax &e) : Bxp(get_id(e.evalmax)),
		fmax(e.fmax),
		best_sol(NULL),
		y_heap(e.y_heap,true),
		nb_bisect(e.nb_bisect+1),
		pu(e.pu),
		evalmax(e.evalmax)
 {
	if (!e.best_sol) {
		best_sol = new IntervalVector(*e.best_sol);
	}

}

long BxpMinMax::get_id(const EvalMax& evalmax) {
	try {
		return ids()[evalmax.id];
	} catch(Map<long,long,false>::NotFound&) {
		long new_id=next_id();
		ids().insert_new(evalmax.id, new_id);
		return new_id;
	}
}



BxpMinMaxSub::BxpMinMaxSub(const EvalMax& evalmax) : Bxp(BxpMinMaxSub::get_id(evalmax)),maxfxy(Interval::all_reals()), feasible(0), evalmax(evalmax) {

}

BxpMinMaxSub::~BxpMinMaxSub() = default;


BxpMinMaxSub::BxpMinMaxSub(const BxpMinMaxSub &e) : Bxp(get_id(e.evalmax)), maxfxy(e.maxfxy),
		feasible(e.feasible), evalmax(e.evalmax) {

}

long BxpMinMaxSub::get_id(const EvalMax& evalmax) {
	try {
		return ids()[evalmax.id];
	} catch(Map<long,long,false>::NotFound&) {
		long new_id=next_id();
		ids().insert_new(evalmax.id, new_id);
		return new_id;
	}

}

//long BxpMinMax::get_id() {
	//    auto *pBxpOpti = dynamic_cast<BxpMinMaxOpti*>(this);
	//    if (pBxpOpti)
		//        return BxpMinMaxOpti::id;
	//    else
		//        return BxpMinMaxCsp::id;
	//}


//    long BxpMinMaxOpti::get_id(const NormalizedSystem& sys) {
//        try {
//            return ids()[sys.original_sys_id];
//        } catch(Map<long,long,false>::NotFound&) {
//            long new_id=next_id();
//            ids().insert_new(sys.original_sys_id, new_id);
//            return new_id;
//        }
//    }
//
//    long BxpMinMaxCsp::get_id(const NormalizedSystem& sys) {
//        try {
//            return ids()[sys.original_sys_id];
//        } catch(Map<long,long,false>::NotFound&) {
//            long new_id=next_id();
//            ids().insert_new(sys.original_sys_id, new_id);
//            return new_id;
//        }
//    }

///* Cost functions for BxpMinMax class */
//
//    CellCostFmaxlb::CellCostFmaxlb(EvalMax &evalmax) : evalmax(evalmax) {
//    }
//
//    double CellCostFmaxlb::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::get_id(evalmax)]);
//        if (pBxp)
//            return pBxp->fmax.lb();
//        else
//            ibex_error("[ibex_BxpMinMax] null ptr in CellCostFmaxlb: cell properties have no BxpMinMax.");
//
////        auto pBxp_opt = dynamic_cast<const BxpMinMaxOpti*>(elem.prop[BxpMinMaxOpti::id]);
////        if (pBxp_opt)
////            return pBxp_opt->fmax.lb();
////        auto pBxp_csp = dynamic_cast<const BxpMinMaxCsp*>(elem.prop[BxpMinMaxCsp::id]);
////        if (pBxp_csp)
////            return pBxp_csp->fmax.lb();
////        else
////            ibex_error("[ibex_BxpMinMax] null ptr in CellCostFmaxlb: cell properties have no BxpMinMax.");
//}
//
//    void CellCostFmaxlb::set_minmax_data(Cell& c) {
//// ************ On le laisse dans le Light ************
////        double new_fmax_ub = y_heap->top1()->get<OptimData>().pf.ub(); // get the upper bound of max f(x,y_heap)
////        double new_fmax_lb = y_heap->top2()->get<OptimData>().pf.lb(); // get the lower bound of max f(x,y_heap)
////
////        //std::cout<<"new_fmax_ub: "<<new_fmax_ub<<std::endl<<"new_fmax_lb: "<<new_fmax_lb<<std::endl<<"fmax_lb (from found point): "<<data_x->fmax.lb()<<std::endl;
////
////        if (new_fmax_ub< new_fmax_lb) {
////            ibex_error("ibex_LightOptimMinMax: error, please report this bug.");
////        }
////
//////                std::cout<<"fmax ini: "<<data_x->fmax<<std::endl;
////        data_x->fmax &= Interval(new_fmax_lb, new_fmax_ub);
//    }
//
//    void CellCostFmaxlb::add_property(BoxProperties &map) {
////        c.prop.add(new BxpMinMax());.
//      if (!map[BxpMinMax::get_id(evalmax)])
//		map.add(new BxpMinMax(evalmax));
//    }
//
//
//        CellCostmaxFmaxub::CellCostmaxFmaxub(EvalMax &evalmax) : evalmax(evalmax) {
//    }
//
//    double CellCostmaxFmaxub::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::get_id(evalmax)]);
//        if (pBxp)
//            return -pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostmaxFmaxub: cell element has no BxpMinMax.");
//    }
//
//    void CellCostmaxFmaxub::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMax());
//    }
//
//    void CellCostmaxFmaxub::add_property(BoxProperties &map) {
//        if (!map[BxpMinMax::get_id(evalmax)])
//            map.add(new BxpMinMax(evalmax));
//    }
//
//        CellCostFmaxub::CellCostFmaxub(EvalMax &evalmax) : evalmax(evalmax) {
//    }
//
//    double CellCostFmaxub::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::get_id(evalmax)]);
//        if (pBxp)
//            return pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostFmaxub: cell element has no BxpMinMax.");
//    }
//
//    void CellCostFmaxub::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMax());
//    }
//
//    void CellCostFmaxub::add_property(BoxProperties &map) {
//      if (!map[BxpMinMax::get_id(evalmax)])
//		map.add(new BxpMinMax(evalmax));
//    }


//    /* Cost functions for BxpMinMaxOpti class */
//        CellCostFmaxlb_opt::CellCostFmaxlb_opt(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostFmaxlb_opt::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxOpti*>(elem.prop[BxpMinMaxOpti::get_id(sys)]);
//        if (pBxp)
//            return pBxp->fmax.lb();
//        else
//            ibex_error("[ibex_BxpMinMax] null ptr in CellCostFmaxlb_opt: cell properties have no BxpMinMaxOpti.");
//}
//
//    void CellCostFmaxlb_opt::set_minmax_data(Cell& c) {
////        c.prop.add(new BxpMinMaxOpti());
//    }
//
//    void CellCostFmaxlb_opt::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxOpti::get_id(sys)])
//		map.add(new BxpMinMaxOpti(sys));
//    }
//
//        CellCostmaxFmaxub_opt::CellCostmaxFmaxub_opt(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostmaxFmaxub_opt::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxOpti*>(elem.prop[BxpMinMaxOpti::get_id(sys)]);
//        if (pBxp)
//            return -pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostmaxFmaxub_opt: cell element has no BxpMinMaxOpti.");
//    }
//
//    void CellCostmaxFmaxub_opt::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMaxOpti());
//    }
//
//    void CellCostmaxFmaxub_opt::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxOpti::get_id(sys)])
//		map.add(new BxpMinMaxOpti(sys));
//    }
//
//        CellCostFmaxub_opt::CellCostFmaxub_opt(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostFmaxub_opt::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxOpti*>(elem.prop[BxpMinMaxOpti::get_id(sys)]);
//        if (pBxp)
//            return pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostFmaxub_opt: cell element has no BxpMinMaxOpti.");
//    }
//
//    void CellCostFmaxub_opt::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMaxOpti());
//    }
//
//    void CellCostFmaxub_opt::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxOpti::get_id(sys)])
//		map.add(new BxpMinMaxOpti(sys));
//    }
//
//
//    /* Cost functions for BxpMinMaxCsp class */
//        CellCostFmaxlb_csp::CellCostFmaxlb_csp(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostFmaxlb_csp::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxCsp*>(elem.prop[BxpMinMaxCsp::get_id(sys)]);
//        if (pBxp)
//            return pBxp->fmax.lb();
//        else
//            ibex_error("[ibex_BxpMinMax] null ptr in CellCostFmaxlb_csp: cell properties have no_csp BxpMinMaxCsp.");
//    }
//
//    void CellCostFmaxlb_csp::set_minmax_data(Cell& c) {
////        c.prop.add(new BxpMinMaxCsp());
//    }
//
//    void CellCostFmaxlb_csp::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxCsp::get_id(sys)])
//		map.add(new BxpMinMaxCsp(sys));
//    }
//
//        CellCostmaxFmaxub_csp::CellCostmaxFmaxub_csp(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostmaxFmaxub_csp::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxCsp*>(elem.prop[BxpMinMaxCsp::get_id(sys)]);
//        if (pBxp)
//            return -pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostmaxFmaxub_csp: cell element has no BxpMinMaxCsp.");
//    }
//
//    void CellCostmaxFmaxub_csp::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMaxCsp());
//    }
//
//    void CellCostmaxFmaxub_csp::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxCsp::get_id(sys)])
//		map.add(new BxpMinMaxCsp(sys));
//    }
//
//        CellCostFmaxub_csp::CellCostFmaxub_csp(ExtendedSystem &sys) : sys(sys) {
//    }
//
//    double CellCostFmaxub_csp::cost(const Cell& elem) const {
//        auto pBxp = dynamic_cast<const BxpMinMaxCsp*>(elem.prop[BxpMinMaxCsp::get_id(sys)]);
//        if (pBxp)
//            return pBxp->fmax.ub();
//        else
//            ibex_error("[ibex_BxpMinMax] -- null ptr in CellCostFmaxub_csp: cell element has no BxpMinMaxCsp.");
//    }
//
//    void CellCostFmaxub_csp::set_minmax_data(Cell &c) {
////        c.prop.add(new BxpMinMaxCsp());
//    }
//
//    void CellCostFmaxub_csp::add_property(BoxProperties &map) {
//      if (!map[BxpMinMaxCsp::get_id(sys)])
//		map.add(new BxpMinMaxCsp(sys));
//    }

} // end namespace ibex
