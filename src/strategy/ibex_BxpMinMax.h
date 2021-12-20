//
// Created by joris on 17/11/2021.
//

#ifndef IBEX_HCONTROL_BXPMINMAX_H
#define IBEX_HCONTROL_BXPMINMAX_H

#include "ibex_Bxp.h"
#include "ibex_Interval.h"
#include "ibex_Cell.h"
#include "ibex_DoubleHeap.h"
#include "ibex_CellCostFunc.h"
#include "ibex_EvalMax.h"

namespace ibex {
    class EvalMax;

    class feasible_point {
    public:
        feasible_point(const Vector& point,const Interval& eval);
        feasible_point(const feasible_point& pt);
        ~feasible_point();

        Vector point;
        Interval eval;
    };


    class CellCostMaxPFub_MinMax: public CostFunc<Cell> {
    public:
        explicit CellCostMaxPFub_MinMax(const EvalMax& evalmax);

        /**
         * \brief Add BxpData
         */
        virtual void add_property(BoxProperties& map);

        /**
         * \brief Set "pf" in BxpData in the cell
         */
        virtual void set_optim_data(Cell& c);

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

        /**
         * \brief Set "pf" in BxpData in the cell
         */
        virtual void set_optim_data(Cell& c);

        /** The "cost" of a element. */
        virtual	double cost(const Cell& c) const;

    protected:

        /**
         * \brief The EvalMax (previously light optim MinMax).
         */
        const EvalMax& evalmax;
    };


    /**
     * \brief Data required for EvalMax
     */
    class BxpMinMax : public Bxp {
    public:
        /**
         * \brief Constructor for the root node (followed by a call to init_root)
         */
        explicit BxpMinMax(EvalMax& evalmax);

        /**
         * \brief Delete *this
         */
        ~BxpMinMax() override;

        /**
         * \brief clear list of feasible by deleting contained pointer on feasible_point object, call only when x_cell discard because not solution. /!\ not when x_cell bisected
         */
        void clear_fsbl_list();

        /**
         * remove pointers from fsbl_point_list that point is not contained in x_box, object pointed at is also deleted if strong_del is true
         */
        void clear_notin_point(const IntervalVector& x_box, bool strong_del);

        /**
         * \brief Duplicate the structure into the left/right nodes
         */
        std::pair<Bxp*, Bxp*> down();

        /**
         * \brief Create a copy.
         */
        virtual Bxp* copy(const IntervalVector& box, const BoxProperties& prop) const;

        /**
         * \brief Update the property upon box modification.
         *
         */
        void update(const BoxEvent& event, const BoxProperties& prop) override;

        /**
         * Enclosure of maximum of the objective function
         */
        Interval fmax;
        IntervalVector *best_sol;

        /**
         * y_heap inherited from father of box
         */
        DoubleHeap<Cell>* y_heap;

        std::vector<feasible_point> fsbl_pt_list;
        long int nb_bisect;

        /**
         * Cost function of the heap to store the element of the light solver
         */
        CellCostMaxPFub_MinMax y_heap_costf1;
        CellCostPFlb_MinMax y_heap_costf2;

        static long get_id(const EvalMax& evalmax);

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
        BxpMinMax(const BxpMinMax& e);
    };


//    /* Inherited classes from DataMinMax */
//
//    class BxpMinMaxOpti : public BxpMinMax {
//    public:
//        explicit BxpMinMaxOpti(ExtendedSystem& sys) : BxpMinMax() {}
//        static long get_id(const NormalizedSystem& sys);
//
//
//    protected:
//        static Map<long, long, false> &ids();
//        Bxp* copy() const {return new BxpMinMaxOpti(*this);};
//    };
//
//    class BxpMinMaxCsp : public BxpMinMax {
//    public:
//        explicit BxpMinMaxCsp(ExtendedSystem& sys) : BxpMinMax(sys) {}
//        static long get_id(const NormalizedSystem& sys);
//
//    protected:
//        Bxp* copy() const {return new BxpMinMaxCsp(*this);};
//        static Map<long, long, false> &ids();
//    };

/*================================== inline implementations ========================================*/

    inline Bxp* BxpMinMax::copy(const IntervalVector& box, const BoxProperties& prop) const {
        return new BxpMinMax(*this);
    }

    inline void BxpMinMax::update(const BoxEvent& event, const BoxProperties& prop) {
        // TODO: we should call compute_pf and compute_pu here or create some
        // "is_up_to_date" flag.
        // This is actually done from the CellCostFunc classes
        // and makes no problem so far as this property is not used elsewhere.
    }


///* Cost function definition for BxpMinMax class */
//
//    class CellCostFmaxlb : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxlb(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };
//
//    class CellCostmaxFmaxub : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostmaxFmaxub(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };
//
//    class CellCostFmaxub : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxub(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };


///* Cost function definition for BxpMinMaxOpti class */
//
//    class CellCostFmaxlb_opt : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxlb_opt(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };
//
//    class CellCostmaxFmaxub_opt : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostmaxFmaxub_opt(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };
//
//    class CellCostFmaxub_opt : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxub_opt(EvalMax& evalmax);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        EvalMax& evalmax;
//    };


///* Cost function definition for BxpMinMaxCsp class */
//
//    class CellCostFmaxlb_csp : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxlb_csp(ExtendedSystem& sys);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        ExtendedSystem& sys;
//    };
//
//    class CellCostmaxFmaxub_csp : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostmaxFmaxub_csp(ExtendedSystem& sys);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        ExtendedSystem &sys;
//    };
//
//    class CellCostFmaxub_csp : public CostFunc<Cell> { // element are sorted from the lowest lb of the evaluation of the objective function to the greatest
//    public:
//
//        explicit CellCostFmaxub_csp(ExtendedSystem& sys);
//
//        /**
//         * \brief Return the cost associated to the cell
//         */
//        double cost(const Cell& elem) const override;
//
//        /**
//         * \brief Add backtrackable data required by this cost function
//         *
//         * This function is called for the root cell (before a strategy is executed).
//         *
//         * Does nothing by default.
//         */
//        void add_property(BoxProperties& map);
//
//        /**
//         * \brief Set data in OptimData in the cell
//         *
//         * This function is called right after "contract_and_bound" in the Optimizer.
//         *
//         * The data required depends on the cost function.
//         *
//         * Does nothing (by default).
//         */
//        static void set_minmax_data(Cell& c);
//
//        ExtendedSystem &sys;
//    };


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
         * \brief Duplicate the structure into the left/right nodes
         */
        std::pair<Bxp*, Bxp*> down();

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

        void compute_pf(const Function& goal, const IntervalVector& box);

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


} // end namespace ibex

#endif //IBEX_HCONTROL_BXPMINMAX_H
