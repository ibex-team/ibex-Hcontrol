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

#ifndef __TEST_EVALMAX_H__
#define __TEST_EVALMAX_H__

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "utils.h"

namespace ibex {

class TestEvalMax : public CppUnit::TestFixture {

public:

	CPPUNIT_TEST_SUITE(TestEvalMax);


	CPPUNIT_TEST(ex1);
	CPPUNIT_TEST_SUITE_END();

	static void ex1();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestEvalMax);


} // namespace ibex
#endif // __TEST_EVALMAX_H__
