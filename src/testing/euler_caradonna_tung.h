#ifndef __EULER_CARADONNA_TUNG_H__
#define __EULER_CARADONNA_TUNG_H__

#include <deal.II/grid/manifold_lib.h>

#include "tests.h"
#include "dg/dg_base.hpp"
#include "physics/physics.h"
#include "parameters/all_parameters.h"

namespace PHiLiP {
namespace Tests {

/// Performs grid convergence for various polynomial degrees.
template <int dim, int nstate>
class CaradonnaTung: public TestsBase
{
public:
    /// Constructor. Deleted the default constructor since it should not be used
    CaradonnaTung () = delete;
    /// Constructor.
    /** Simply calls the TestsBase constructor to set its parameters = parameters_input
     */
    CaradonnaTung(const Parameters::AllParameters *const parameters_input,
                  const dealii::ParameterHandler &parameter_handler_input);


    /// Parameter handler for storing the .prm file being ran
    const dealii::ParameterHandler &parameter_handler;

    /// Grid convergence on Euler Caradonna Tung
    /** Will run the a grid convergence test for various p
     *  on multiple grids to determine the order of convergence.
     *
     *  Expecting the solution to converge at p+1. and output to converge at 2p+1.
     *  Note that the output solution currently convergens slightly suboptimally
     *  depending on the case (around 2p). The implementation of the boundary conditions
     *  play a large role on this adjoint consistency.
     *
     *  Want to see entropy go to 0.
     */
    int run_test () const;

protected:

    //  // Integrate entropy over the entire domain to use as a functional.
    //  double integrate_entropy_over_domain(DGBase<dim,double> &dg) const;
};


//   /// Manufactured grid convergence
//   /** Currently the main function as all my test cases simply
//    *  check for optimal convergence of the solution
//    */
//   template<int dim>
//   int manufactured_grid_convergence (Parameters::AllParameters &parameters);

} // Tests namespace
} // PHiLiP namespace
#endif
