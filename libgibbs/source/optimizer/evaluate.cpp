/*=============================================================================
	Copyright (c) 2012-2013 Richard Otis

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

// evaluate.cpp -- evaluate energies from a Database

#include "libgibbs/include/libgibbs_pch.hpp"
#include "libtdb/include/database.hpp"
#include "libtdb/include/structure.hpp"
#include "libgibbs/include/optimizer/optimizer.hpp"
#include "libgibbs/include/optimizer/opt_Gibbs.hpp"
#include "external/coin/IpIpoptApplication.hpp"
#include "external/coin/IpSolveStatistics.hpp"
#include <iostream>
#include <fstream>

using namespace Ipopt;

void evaluate(const Database &DB, const evalconditions &conditions) {
	Phase_Collection phase_col;
	// TODO: temporary code for suspending all phases but FCC and liquid
	for (auto i = DB.get_phase_iterator(); i != DB.get_phase_iterator_end(); ++i) {
		if (i->first == "FCC_A1" || i->first == "LIQUID") {
			phase_col[i->first] = i->second;
		}
	}

	const Phase_Collection::const_iterator phase_iter = phase_col.cbegin();
	const Phase_Collection::const_iterator phase_end = phase_col.cend();

	// TODO: check validity of conditions

	// Create an instance of your nlp...
	SmartPtr<TNLP> mynlp = new GibbsOpt(phase_iter, phase_end, conditions);

	// Create an instance of the IpoptApplication
	//
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	app->Options()->SetStringValue("derivative_test","first-order");
	app->Options()->SetStringValue("hessian_approximation","limited-memory");
	//app->Options()->SetIntegerValue("print_level",12);
	app->Options()->SetStringValue("derivative_test_print_all","yes");

	// Initialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
		std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
		return;
	}

	status = app->OptimizeTNLP(mynlp);

	if (status == Solve_Succeeded) {
		// Retrieve some statistics about the solve
		Index iter_count = app->Statistics()->IterationCount();
		std::cout << std::endl << std::endl << "*** The problem solved in " << iter_count << " iterations!" << std::endl;

		Number final_obj = app->Statistics()->FinalObjective();
		std::cout << std::endl << std::endl << "*** The final value of the objective function is " << final_obj << '.' << std::endl;
	}
}