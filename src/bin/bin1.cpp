#include "ibex.h"
#include "ibex_EvalMax.h"
#include "ibex_SystemFactory.h"
#include "ibex_CtcIdentity.h"

using namespace std;
using namespace ibex;

int main (int argc, char *argv[]) {

    double y_prec = 1e-4;
    double stop_prec = 1e-2;

    Variable x(1), y(1);
    IntervalVector x_ini(1,Interval(-20,20));
    IntervalVector y_ini(1,Interval(-20,20));
    Function func(x,y,(pow(x,2)-pow(y,2)));

    SystemFactory x_fac;
    x_fac.add_var(x, x_ini);
    System x_sys(x_fac);
    CtcIdentity x_ctc(x_ini.size());

    SystemFactory xy_fac;
    xy_fac.add_var(x, x_ini);
    xy_fac.add_var(y, y_ini);
    xy_fac.add_goal(func);
    System xy_sys(xy_fac);
    CtcIdentity xy_ctc(x_ini.size()+y_ini.size());
    cout << xy_sys << endl;
    //EvalMax ex1(y_ini,xy_sys, xy_ctc);
    EvalMax ex1(y_ini,func);

    ex1.set_timeout(100);
    ex1.set_monitor(true);
    ex1.set_max_nb_iter(10);
   // ex1.set_visit_all(true);
    ex1.set_prec_y(y_prec);
    ex1.set_prec_y(0);
    ex1.set_goal_rel_prec(stop_prec);

    Cell *x_cell=new Cell(x_ini);
    Interval res1 = ex1.eval(x_cell->box,x_cell->prop);
    cout << "result: " << res1 << endl;

    LargestFirst bsc;

    std::pair<Cell*,Cell*> subcells_pair = bsc.bisect(*x_cell);
    delete x_cell;
    Interval res2 = ex1.eval(subcells_pair.first->box,subcells_pair.first->prop);
    cout << subcells_pair.first->box <<"  result: " << res2 << endl;

    delete subcells_pair.first;
    Interval res3 = ex1.eval(subcells_pair.second->box,subcells_pair.second->prop);
    cout << subcells_pair.second->box <<"  result: " << res3 << endl;
    delete subcells_pair.second;
    return 0;
}

