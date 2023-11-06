#include "entropy_adjoint.hpp"

namespace PHiLiP {

template <int dim,int nstate,typename real,typename MeshType>
EntropyFunctional<dim,nstate,real,MeshType>
::EntropyFunctional(
    std::shared_ptr<DGBase<dim,real,MeshType>> dg_input)
    : Functional<dim,nstate,real,MeshType>(dg_input)
    , euler_fad_fad(dynamic_cast< Physics::Euler<dim,dim+2,FadFadType> &>(*(this->physics_fad_fad)))
{}

template <int dim,int nstate,typename real,typename MeshType>
real EntropyFunctional<dim,nstate,real,MeshType>
::evaluate_functional(const bool compute_dIdW, 
                      const bool compute_dIdX, 
                      const bool compute_d2I)
{
    double value = Functional<dim,nstate,real,MeshType>::evaluate_functional(compute_dIdW, compute_dIdX, compute_d2I);
    this->pcout << "Entropy value: " << value << "\n";

    return value;
}

#if PHILIP_DIM!=1
template class EntropyFunctional <PHILIP_DIM, PHILIP_DIM+2, double, dealii::parallel::distributed::Triangulation<PHILIP_DIM>>;
#endif

} // PHiLiP namespace