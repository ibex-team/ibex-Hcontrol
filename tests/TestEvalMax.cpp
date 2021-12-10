/* ============================================================================
 * I B E X - EvalMax Tests
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Joris Tillet
 * Created     : Nov 29, 2021
 * ---------------------------------------------------------------------------- */

#include "TestEvalMax.h"
#include "ibex_EvalMax.h"
#include "ibex_SystemFactory.h"
#include "ibex_CtcIdentity.h"

using namespace std;

namespace ibex {


void TestEvalMax::ex1() {

	double x_prec = 1e-6;
	double y_prec = 1e-4;
	double stop_prec = 1e-2;

	Variable x(1), y(1);
	IntervalVector x_ini(1,Interval(-20,20));
	IntervalVector y_ini(1,Interval(10,20));
	Function func(x,y,(pow(x,2)-pow(y,2)));

	SystemFactory x_fac;
	x_fac.add_var(x, x_ini);
	NormalizedSystem x_sys(x_fac);
	CtcIdentity x_ctc(x_ini.size());

	SystemFactory xy_fac;
	xy_fac.add_var(x, x_ini);
	xy_fac.add_var(y, y_ini);
	xy_fac.add_goal(func);
	NormalizedSystem xy_sys(xy_fac);
	CtcIdentity xy_ctc(x_ini.size()+y_ini.size());

    EvalMax ex1(xy_sys, 1, 1, xy_ctc);
    ex1.timeout = 100;
    ExtendedSystem sys(xy_sys);
    auto bxpties = BoxProperties(x_ini);
    BxpMinMaxOpti bxpmm(sys);
    auto bxp = dynamic_cast<Bxp*>(&bxpmm);
    if (!bxp) ibex_error("casting error");
//    bxpties.add(bxp);
    auto res = ex1.eval(x_sys.box, bxpties, 1000);
    cout << "result: " << res << endl;

    CPPUNIT_ASSERT(res.contains(-100)); // TODO

//	CPPUNIT_ASSERT(res==Optim::SUCCESS);
//	CPPUNIT_ASSERT(res.loup<=0.5286);
//	CPPUNIT_ASSERT(res.nb_cells<=9000);

}


} // end namespace
