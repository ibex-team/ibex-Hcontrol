//============================================================================
//                                  I B E X                                   
// File        : ibex_CellFuncMinMax.h
// Author      : Jordan Ninin
// Copyright   : ENSTA Bretagne
// License     : See the LICENSE file
// Created     : Dec 12, 2021
//============================================================================

#ifndef __IBEX_CELLFUNC_MINMAX_H__
#define __IBEX_CELLFUNC_MINMAX_H__


#include "ibex_Cell.h"
#include "ibex_Heap.h"
#include "ibex_Bxp.h"


namespace ibex {
class EvalMax;



/**
 * \brief Data required for ordering the elements in y_heap of BxpMinMax
 */
class BxpMinMaxSub : public Bxp {
public:
	/**
	 * \brief Constructor for the root node (followed by a call to init_root)
	 */
	explicit BxpMinMaxSub(const EvalMax& evalmax);

	/**
	 * \brief Delete *this
	 */
	~BxpMinMaxSub() override;

	/**
	 * \brief Create a copy.
	 */
	virtual Bxp* copy(const IntervalVector& box, const BoxProperties& prop) const;

	/**
	 * \brief Update the property upon box modification.
	 *
	 */
	void update(const BoxEvent& event, const BoxProperties& prop) override;


	static long get_id(const EvalMax& evalmax);


	/**
	 * \brief Value of "max f(x,y)"
	 *
	 * Image of the objective on the current box
	 */
	Interval maxfxy;

	/**
	 * \brief to know if  all constraints are satisfed
	 *
	 * 0 unknown and 1 feasible i.e. all constraints satisfied.
	 */
	bool feasible;

protected:

	/**
	 * \brief The EvalMax (previously light optim MinMax).
	 */
	const EvalMax& evalmax;
	static Map<long, long, false> &ids();

	/**
	 * \brief Constructor by copy.
	 */
	BxpMinMaxSub(const BxpMinMaxSub& e);
};







// Cost function to sort the heap with the Maximum of the value of "pf.ub()"
class CellCostMaxPFub_MinMax: public CostFunc<Cell> {
public:
	explicit CellCostMaxPFub_MinMax(const EvalMax& evalmax);

	/**
	 * \brief Add BxpData
	 */
	virtual void add_property(BoxProperties& map);

	/** The "cost" of a element. */
	virtual	double cost(const Cell& c) const;

private:

	/**
	 * \brief The EvalMax (previously light optim MinMax).
	 */
	const EvalMax& evalmax;
};




// Cost function to sort the heap with the Minimum of the value of "pf.lb()"
class CellCostPFlb_MinMax: public CostFunc<Cell> {
public:
	CellCostPFlb_MinMax(const EvalMax& evalmax);

	/**
	 * \brief Add BxpData
	 */
	virtual void add_property(BoxProperties& map);


	/** The "cost" of a element. */
	virtual	double cost(const Cell& c) const;

private:

	/**
	 * \brief The EvalMax (previously light optim MinMax).
	 */
	const EvalMax& evalmax;
};



/*================================== inline implementations ========================================*/


/*================================== inline implementations ========================================*/


inline Bxp* BxpMinMaxSub::copy(const IntervalVector& box, const BoxProperties& prop) const {
	return new BxpMinMaxSub(*this);
}

inline void BxpMinMaxSub::update(const BoxEvent& event, const BoxProperties& prop) {
	// TODO: we should call compute_pf and compute_pu here or create some
	// "is_up_to_date" flag.
	// This is actually done from the CellCostFunc classes
	// and makes no problem so far as this property is not used elsewhere.
}


} // namespace ibex

#endif // __IBEX_CELL_DOUBLE_HEAP_H__
