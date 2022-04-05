	//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderXTaylor.cpp
// Author      : Gilles Chabert, Ignacio Araya, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 12, 2012
// Last Update : Jul 09, 2017
//============================================================================

#include "ibex_LoupFinderMinMax.h"

using namespace std;

namespace ibex {

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderMinMax::LoupFinderMinMax(EvalMax& eval, const System& sys) :
		sys(sys),
		evalgoal(eval) {

}

void LoupFinderMinMax::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	evalgoal.add_property(init_box,prop);
}

std::pair<IntervalVector, double> LoupFinderMinMax::find(const IntervalVector& box, const IntervalVector& old_loup_point, double current_loup, BoxProperties& prop) {


	throw NotFound();
}

} /* namespace ibex */
