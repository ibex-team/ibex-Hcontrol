//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderMinMax.h
// Author      : Jordan Ninin
// Copyright   : ENSTA Bretagne (France)
// License     : See the LICENSE file
// Created     : Jan 11, 2022
// Last Update : Jan 11, 2022
//============================================================================

#ifndef __IBEX_LOUP_FINDER_MINMAX_H__
#define __IBEX_LOUP_FINDER_MINMAX_H__

#include "ibex_LoupFinder.h"
#include "ibex_EvalMax.h"

namespace ibex {
/**
 * \ingroup optim
 *
 * \brief Upper-bounding algorithm for MinMax problem.
 *
 * \note ancien LightLocalSolver
 */

class LoupFinderMinMax : public LoupFinder {
public:

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * \param evalgoal    - The goal evaluation function.
	 * \param sys         - The constraints of the problem.
	 */
	LoupFinderMinMax(EvalMax& evalgoal, const System& sys);
	LoupFinderMinMax(EvalMax& evalgoal);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup, BoxProperties& prop);

	/**
	 * \brief Add properties required by this loup finder.
	 */
	virtual void add_property(const IntervalVector& init_box, BoxProperties& prop);


protected:
	/**
	 * \brief The constraints of the problem..
	 */
	const System& sys;

	/**
	 * \brief The eval function of the objective function
	 */
	EvalMax& evalgoal;
};

/*============================================ inline implementation ============================================ */

inline std::pair<IntervalVector, double> LoupFinderMinMax::find(const IntervalVector& box, const IntervalVector& loup_point, double loup) {
	BoxProperties prop(box);
	return find(box, loup_point, loup, prop);
}

} /* namespace ibex */



#endif /* __IBEX_LOUP_FINDER_X_TAYLOR_H__ */
