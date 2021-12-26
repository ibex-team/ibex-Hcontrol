#include "ibex.h"
#include "ibex_EvalMax.h"
#include "ibex_SystemFactory.h"
#include "ibex_CtcIdentity.h"

using namespace std;
using namespace ibex;

int main (int argc, char *argv[]) {

    double x_prec = 1e-6;
    double y_prec = 1e-4;
    double stop_prec = 1e-2;

    Variable x(1), y(1);
    IntervalVector x_ini(1,Interval(-20,20));
    IntervalVector y_ini(1,Interval(10,20));
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

    EvalMax ex1(y_ini,xy_sys, xy_ctc);
    ex1.timeout = 100;
  //  ex1.prec_y = y_prec;
    Interval res = ex1.eval(x_sys.box);
    cout << "result: " << res << endl;

    return 0;
}

