//============================================================================
//                                  I B E X                                   
// File        : ibex_CellMinMaxHeap.h
// Author      : Gilles Chabert, Jordan Ninin
// Copyright   : ENSTA Bretagne
// License     : See the LICENSE file
// Created     : Dec 12, 2021
//============================================================================

#ifndef __IBEX_CELL_MINMAX_HEAP_H__
#define __IBEX_CELL_MINMAX_HEAP_H__

#include "ibex_DoubleHeap.h"
#include "ibex_Cell.h"
#include "ibex_CellBufferOptim.h"
#include "ibex_Bxp.h"
#include "ibex_EvalMax.h"


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
	 * \brief Casado criterion
	 *
	 * Image of the objective on the current box
	 */
	Interval pf;

	/**
	 * \brief Casado criterion
	 *
	 * Constraint factor of the current box : between 0 infeasible and 1 for all constraints satisfied.
	 */
	double pu;

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







class CellCostMaxPFub_MinMax: public CostFunc<Cell> {
public:
	explicit CellCostMaxPFub_MinMax(const EvalMax& evalmax);

	/**
	 * \brief Add BxpData
	 */
	virtual void add_property(BoxProperties& map);

	//        /**
	//         * \brief Set "pf" in BxpData in the cell
	//         */
	//        virtual void set_optim_data(Cell& c);

	/** The "cost" of a element. */
	virtual	double cost(const Cell& c) const;

protected:

	/**
	 * \brief The EvalMax (previously light optim MinMax).
	 */
	const EvalMax& evalmax;
};




class CellCostPFlb_MinMax: public CostFunc<Cell> {
public:
	CellCostPFlb_MinMax(const EvalMax& evalmax);

	/**
	 * \brief Add BxpData
	 */
	virtual void add_property(BoxProperties& map);
	//
	//        /**
	//         * \brief Set "pf" in BxpData in the cell
	//         */
	//        virtual void set_optim_data(Cell& c);

	/** The "cost" of a element. */
	virtual	double cost(const Cell& c) const;

protected:

	/**
	 * \brief The EvalMax (previously light optim MinMax).
	 */
	const EvalMax& evalmax;
};



/**
 * \ingroup optim
 *
 * \brief Double-heap buffer (for global optimization of MinMax problem)
 *
 * This is a double-heap buffer where the first heap criterion
 * is LB (lower bound of the objective domain) and the second
 * heap criterion is set by the user (default is UB).
 *
 */
class CellMinMaxHeap : public DoubleHeap<Cell>, public CellBufferOptim {

public:

	/**
	 * \brief Create the buffer.
	 *
	 * \param evalmax    - the EvalMax
	 * \param crit   - second criterion in node selection (the first criterion is the
	 *                 minimum of the objective estimate). default value CellHeapOPtim::UB.
	 */
	CellMinMaxHeap(const EvalMax& evalmax, int crit2_pr=50);

	/**
	 * \brief Copy constructor.
	 */
	explicit CellMinMaxHeap(const CellMinMaxHeap& dhcp, bool deep_copy=false);

	/**
	 * \brief Delete *this.
	 */
	~CellMinMaxHeap();

	/**
	 * \brief Add properties required by this buffer.
	 */
	virtual void add_property(const IntervalVector& init_box, BoxProperties& map);

	/**
	 * \brief Flush the buffer.
	 *
	 * All the remaining cells will be *deleted*
	 */
	void flush();

	/** \brief Return the size of the buffer. */
	unsigned int size() const;

	/** \brief Return true if the buffer is empty. */
	bool empty() const;

	/** \brief Push a new cell on the stack. */
	void push(Cell* cell);

	/** \brief Pop a cell from the stack and return it.*/
	Cell* pop();

	/** \brief Return the next box (but does not pop it).*/
	Cell* top() const;


	std::ostream& print(std::ostream& os) const;

	/**
	 * \brief Return the minimum value of the heap
	 *
	 */
	virtual double minimum() const;

	/**
	 * \brief Contract the heap
	 *
	 * Removes (and deletes) from the heap all the cells
	 * with a cost (according to the cost function of the
	 * first heap) greater than \a loup.
	 */
	virtual void contract(double loup);

	/**
	 * \brief Cost function of the first heap
	 */
	CellCostMaxPFub_MinMax& cost1();

	/**
	 * \brief Cost function of the second heap
	 */
	CellCostPFlb_MinMax& cost2();

protected:

	/**
	 * The system
	 */
	const EvalMax& evalmax;
};

/*================================== inline implementations ========================================*/


inline CellCostMaxPFub_MinMax::CellCostMaxPFub_MinMax(const EvalMax& evalmax) : evalmax(evalmax) { }

inline void CellCostMaxPFub_MinMax::add_property(BoxProperties& map) {
	if (!map[BxpMinMaxSub::get_id(evalmax)])
		map.add(new BxpMinMaxSub(evalmax));
}

//
//    void CellCostMaxPFub_MinMax::set_optim_data(Cell& c) {
//        ((BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)])->compute_pf(*evalmax.xy_sys.goal,c.box);
//    }

inline double CellCostMaxPFub_MinMax::cost(const Cell& c) const {
	const BxpMinMaxSub *data = (BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)];
	if (data) {
		return -data->pf.ub();
	} else {
		ibex_error("CellCostMaxPFub_MinMax::cost : invalid cost");
		return POS_INFINITY;
	}
}


inline CellCostPFlb_MinMax::CellCostPFlb_MinMax(const EvalMax& evalmax) : evalmax(evalmax) { }


inline void CellCostPFlb_MinMax::add_property(BoxProperties& map) {
	if (!map[BxpMinMaxSub::get_id(evalmax)])
		map.add(new BxpMinMaxSub(evalmax));
}

//    void CellCostPFlb_MinMax::set_optim_data(Cell& c) {
//        ((BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)])->compute_pf(*evalmax.xy_sys.goal,c.box);
//    }

inline double CellCostPFlb_MinMax::cost(const Cell& c) const {
	const BxpMinMaxSub *data = (BxpMinMaxSub*) c.prop[BxpMinMaxSub::get_id(evalmax)];
	if (data) {
		return  data->pf.lb();
	} else {
		ibex_error("CellCostPFlb_MinMax::cost : invalid cost"); return POS_INFINITY;
	}
}


inline CellMinMaxHeap::CellMinMaxHeap(const EvalMax& evalmax, int crit2_pr) :
				DoubleHeap<Cell>(*new CellCostMaxPFub_MinMax(evalmax), false,
						*new CellCostPFlb_MinMax(evalmax), false, crit2_pr),
						evalmax(evalmax) {
}


inline CellMinMaxHeap::CellMinMaxHeap(const CellMinMaxHeap& dhcp, bool deep_copy) :
				DoubleHeap<Cell>(dhcp, deep_copy), evalmax(dhcp.evalmax) {
}


inline CellMinMaxHeap::~CellMinMaxHeap() {
	flush();
	delete &cost1();
	delete &cost2();
}

inline void CellMinMaxHeap::contract(double new_loup) {

	DoubleHeap<Cell>::contract(new_loup);
}

inline CellCostMaxPFub_MinMax& CellMinMaxHeap::cost1()      { return (CellCostMaxPFub_MinMax&) heap1->costf; }

inline CellCostPFlb_MinMax& CellMinMaxHeap::cost2()      { return (CellCostPFlb_MinMax&) heap2->costf; }

inline void CellMinMaxHeap::add_property(const IntervalVector& init_box, BoxProperties& map) {
	// add data "pu" and "pf" (if required)
	cost2().add_property(map);
}

inline void CellMinMaxHeap::flush()               { DoubleHeap<Cell>::flush(); }

inline unsigned int CellMinMaxHeap::size() const  { return DoubleHeap<Cell>::size(); }

inline bool CellMinMaxHeap::empty() const         { return DoubleHeap<Cell>::empty(); }

inline void CellMinMaxHeap::push(Cell* cell) {
	// we know cost1() does not require OptimData
	//cost2().set_optim_data(*cell); // TODO a verifier qu'on peut l'enlever

	// the cell is put into the 2 heaps
	DoubleHeap<Cell>::push(cell);


}

inline Cell* CellMinMaxHeap::pop()                { return DoubleHeap<Cell>::pop(); }
inline Cell* CellMinMaxHeap::top() const          { return DoubleHeap<Cell>::top(); }

inline double CellMinMaxHeap::minimum() const     { return DoubleHeap<Cell>::minimum(); }

inline std::ostream& CellMinMaxHeap::print(std::ostream& os) const {
	os << "==============================================================================\n";
	if (empty()) {
		return os << " EMPTY heap" << std::endl;
	} else {
		os << " first heap " << " size " << heap1->size() << " top " << heap1->top()->box << std::endl;
		os << " second heap " << " size " << heap2->size() << " top " << heap2->top()->box ;
		return  os << std::endl;
	}
}


} // namespace ibex

#endif // __IBEX_CELL_DOUBLE_HEAP_H__
