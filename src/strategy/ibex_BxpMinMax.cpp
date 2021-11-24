//
// Created by joris on 17/11/2021.
//

#include "ibex_BxpMinMax.h"

namespace ibex {


feasible_point::feasible_point(const Vector& box,const Interval& eval) : point(box), eval(eval) {}
feasible_point::feasible_point(const feasible_point& pt) = default;
feasible_point::~feasible_point() = default;

const long BxpMinMax::id = next_id();

BxpMinMax::BxpMinMax() : Bxp(id),
                         y_heap(new DoubleHeap<Cell>(y_heap_costf1, false, y_heap_costf2, false)),
                         nb_bisect(0), pu(0) {

}

BxpMinMax::~BxpMinMax() {
    if (y_heap != nullptr) {
//        std::cout<<" flushing y_heap, which contains"<<y_heap->size()<< "elements"<<std::endl;
        y_heap->flush();
//        std::cout<<"y_heap flushed"<<std::endl;
        delete y_heap;
    }
//    std::cout<<"yheap flushed, delete best sol..."<<std::endl;
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

std::pair<Bxp*, Bxp*> BxpMinMax::down() {
    return {new BxpMinMax(*this), new BxpMinMax(*this)};
}

BxpMinMax::BxpMinMax(const BxpMinMax &e) : Bxp(e.id), y_heap(e.y_heap), nb_bisect(e.nb_bisect), pu(e.pu) {

}

//long BxpMinMax::get_id() {
//    auto *pBxpOpti = dynamic_cast<BxpMinMaxOpti*>(this);
//    if (pBxpOpti)
//        return BxpMinMaxOpti::id;
//    else
//        return BxpMinMaxCsp::id;
//}

    double CellCostFmaxlb::cost(const Cell& elem) const {
        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::id]);
        if (pBxp)
            return pBxp->fmax.lb();
        else
            ibex_error("ibex_BxpMinMax -- null ptr in CellCostFmaxlb: cell element has no BxpMinMax.");
}

    void CellCostFmaxlb::set_optim_data(Cell& c) {
        c.prop.add(new BxpMinMax());
    }

    double CellCostmaxFmaxub::cost(const Cell& elem) const {
        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::id]);
        if (pBxp)
            return -pBxp->fmax.ub();
        else
            ibex_error("ibex_BxpMinMax -- null ptr in CellCostmaxFmaxub: cell element has no BxpMinMax.");
    }

    void CellCostmaxFmaxub::set_optim_data(Cell &c) {
        c.prop.add(new BxpMinMax());
    }


    double CellCostFmaxub::cost(const Cell& elem) const {
        auto pBxp = dynamic_cast<const BxpMinMax*>(elem.prop[BxpMinMax::id]);
        if (pBxp)
            return pBxp->fmax.ub();
        else
            ibex_error("ibex_BxpMinMax -- null ptr in CellCostFmaxub: cell element has no BxpMinMax.");
    }

    void CellCostFmaxub::set_optim_data(Cell &c) {
        c.prop.add(new BxpMinMax());
    }

} // end namespace ibex