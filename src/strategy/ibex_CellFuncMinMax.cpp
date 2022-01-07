//============================================================================
//                                  I B E X                                   
// File        : ibex_CellFuncMinMax.h
// Author      : Jordan Ninin
// Copyright   : ENSTA Bretagne
// License     : See the LICENSE file
// Created     : Dec 12, 2021
//============================================================================

#include "ibex_CellFuncMinMax.h"
#include "ibex_EvalMax.h"


namespace ibex {



BxpMinMaxSub::BxpMinMaxSub(const EvalMax& evalmax) : Bxp(BxpMinMaxSub::get_id(evalmax)),maxfxy(Interval::all_reals()), feasible(0), evalmax(evalmax) {

}

BxpMinMaxSub::~BxpMinMaxSub() {};


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



CellCostMaxPFub_MinMax::CellCostMaxPFub_MinMax(const EvalMax& evalmax) : evalmax(evalmax) { }

void CellCostMaxPFub_MinMax::add_property(BoxProperties& map) {
	if (!map[BxpMinMaxSub::get_id(evalmax)])
		map.add(new BxpMinMaxSub(evalmax));
}

//
//    void CellCostMaxPFub_MinMax::set_optim_data(Cell& c) {
//        ((BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)])->compute_pf(*evalmax.xy_sys.goal,c.box);
//    }

double CellCostMaxPFub_MinMax::cost(const Cell& c) const {
	const BxpMinMaxSub *data = (BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)];
	if (data) {
		return -data->maxfxy.ub();
	} else {
		ibex_error("CellCostMaxPFub_MinMax::cost : invalid cost");
		return POS_INFINITY;
	}
}


CellCostPFlb_MinMax::CellCostPFlb_MinMax(const EvalMax& evalmax) : evalmax(evalmax) { }


void CellCostPFlb_MinMax::add_property(BoxProperties& map) {
	if (!map[BxpMinMaxSub::get_id(evalmax)])
		map.add(new BxpMinMaxSub(evalmax));
}

//    void CellCostPFlb_MinMax::set_optim_data(Cell& c) {
//        ((BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)])->compute_pf(*evalmax.xy_sys.goal,c.box);
//    }

double CellCostPFlb_MinMax::cost(const Cell& c) const {
	const BxpMinMaxSub *data = (BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)];
	if (data) {
		return  data->maxfxy.lb();
	} else {
		ibex_error("CellCostPFlb_MinMax::cost : invalid cost"); return POS_INFINITY;
	}
}



} // namespace ibex

