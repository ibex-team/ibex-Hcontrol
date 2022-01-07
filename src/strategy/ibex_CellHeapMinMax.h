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
#include "ibex_CellFuncMinMax.h"



namespace ibex {
class EvalMax;
//class CellCostMaxPFub_MinMax;
//class CellCostPFlb_MinMax;

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
class CellHeapMinMax : public DoubleHeap<Cell>, public CellBufferOptim {

public:

	/**
	 * \brief Create the buffer.
	 *
	 * \param evalmax    - the EvalMax
	 * \param crit   - second criterion in node selection (the first criterion is the
	 *                 minimum of the objective estimate). default value CellHeapOPtim::UB.
	 */
	CellHeapMinMax(const EvalMax& evalmax);

	/**
	 * \brief Copy constructor.
	 */
	explicit CellHeapMinMax(const CellHeapMinMax& dhcp, bool deep_copy=false);

	/**
	 * \brief Delete *this.
	 */
	~CellHeapMinMax();

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

private:

	/**
	 * The system
	 */
	const EvalMax& evalmax;

};

/*================================== inline implementations ========================================*/


inline CellHeapMinMax::~CellHeapMinMax() {
	flush();
}

inline void CellHeapMinMax::contract(double new_loup) {

	DoubleHeap<Cell>::contract(new_loup);
}

inline void CellHeapMinMax::flush()               { DoubleHeap<Cell>::flush(); }

inline unsigned int CellHeapMinMax::size() const  { return DoubleHeap<Cell>::size(); }

inline bool CellHeapMinMax::empty() const         { return DoubleHeap<Cell>::empty(); }

inline void CellHeapMinMax::push(Cell* cell) {
	// we know cost1() does not require OptimData
	//cost2().set_optim_data(*cell); // TODO a verifier qu'on peut l'enlever

	// the cell is put into the 2 heaps
	DoubleHeap<Cell>::push(cell);


}

inline Cell* CellHeapMinMax::pop()                { return DoubleHeap<Cell>::pop(); }
inline Cell* CellHeapMinMax::top() const          { return DoubleHeap<Cell>::top(); }

inline double CellHeapMinMax::minimum() const     { return DoubleHeap<Cell>::minimum(); }

inline std::ostream& CellHeapMinMax::print(std::ostream& os) const {
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
