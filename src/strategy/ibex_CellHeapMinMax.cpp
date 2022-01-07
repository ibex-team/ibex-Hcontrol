//============================================================================
//                                  I B E X                                   
// File        : ibex_CellMinMaxHeap.h
// Author      : Gilles Chabert, Jordan Ninin
// Copyright   : ENSTA Bretagne
// License     : See the LICENSE file
// Created     : Dec 12, 2021
//============================================================================


#include "ibex_CellHeapMinMax.h"

#include "ibex_EvalMax.h"

namespace ibex {


CellHeapMinMax::CellHeapMinMax(const EvalMax& evalmax) :
				DoubleHeap<Cell>(*(evalmax.cost1), false, *(evalmax.cost2), false, evalmax.crit_heap),
						evalmax(evalmax) {
}


CellHeapMinMax::CellHeapMinMax(const CellHeapMinMax& dhcp, bool deep_copy) :
				DoubleHeap<Cell>(dhcp, deep_copy), evalmax(dhcp.evalmax) {
}


void CellHeapMinMax::add_property(const IntervalVector& init_box, BoxProperties& map) {
	// add data "pu" and "pf" (if required)
	cost2().add_property(map);
}


inline CellCostMaxPFub_MinMax& CellHeapMinMax::cost1()      { return (CellCostMaxPFub_MinMax&) heap1->costf; }

inline CellCostPFlb_MinMax& CellHeapMinMax::cost2()      { return (CellCostPFlb_MinMax&) heap2->costf; }

} // namespace ibex

