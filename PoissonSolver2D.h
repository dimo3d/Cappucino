#pragma once
#include <openvdb/tools/PoissonSolver.h>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tools {
namespace poisson {
	template<
		typename PreconditionerType,
		typename TreeType,
		typename BoundaryOp,
		typename SurfaceNormalTreeType,
		typename Interrupter>
		inline typename TreeType::Ptr
		solveWithBoundaryConditionsAndPreconditioner2D(
			const TreeType&,
			const BoundaryOp&,
			const SurfaceNormalTreeType&,
			math::pcg::State&,
			Interrupter&,
			bool staggered = false);

	
	template<
		typename PreconditionerType,
		typename TreeType,
		typename BoundaryOp, typename SurfaceNormalTreeType,
		typename Interrupter>
		inline typename TreeType::Ptr
		solveWithBoundaryConditionsAndPreconditioner2D(
			const TreeType& inTree,
			const BoundaryOp& boundaryOp, const SurfaceNormalTreeType& surfaceNormalTree,
			math::pcg::State& state,
			Interrupter& interrupter,
			bool staggered)
	{
		return solveWithBoundaryConditionsAndPreconditioner2D<PreconditionerType>(
			/*source=*/inTree, /*domain mask=*/inTree, boundaryOp, surfaceNormalTree, state, interrupter, staggered);
	};
	
	template<
		typename PreconditionerType,
		typename TreeType,
		typename DomainTreeType,
		typename BoundaryOp,
		typename SurfaceNormalTreeType,
		typename Interrupter>
		inline typename TreeType::Ptr
		solveWithBoundaryConditionsAndPreconditioner2D(
			const TreeType& inTree,
			const DomainTreeType& domainMask,
			const BoundaryOp& boundaryOp,
			const SurfaceNormalTreeType& surfaceNormalTree,
			math::pcg::State& state,
			Interrupter& interrupter,
			bool staggered)
	{
		using TreeValueT = typename TreeType::ValueType;
		using VecValueT = LaplacianMatrix::ValueType;
		using VectorT = typename math::pcg::Vector<VecValueT>;
		using VIdxTreeT = typename TreeType::template ValueConverter<VIndex>::Type;
		using MaskTreeT = typename TreeType::template ValueConverter<bool>::Type;

		// 1. Create a mapping from active voxels of the input tree to elements of a vector.
		typename VIdxTreeT::ConstPtr idxTree = createIndexTree(domainMask);

		// 2. Populate a vector with values from the input tree.
		typename VectorT::Ptr b = createVectorFromTree<VecValueT>(inTree, *idxTree);

		// 3. Create a mask of the interior voxels of the input tree (from the densified index tree).
		/// @todo Is this really needed?
		typename MaskTreeT::Ptr interiorMask(
			new MaskTreeT(*idxTree, /*background=*/false, TopologyCopy()));
		tools::erodeVoxels(*interiorMask, /*iterations=*/1, tools::NN_FACE);

		// 4. Create the Laplacian matrix.
		LaplacianMatrix::Ptr laplacian = createISLaplacianWithBoundaryConditions(
			*idxTree, *interiorMask, boundaryOp, *b, staggered);

		// 5. Solve the Poisson equation.
		laplacian->scale(-1.0); // matrix is negative-definite; solve -M x = -b
		b->scale(-1.0);
		typename VectorT::Ptr x(new VectorT(b->size(), zeroVal<VecValueT>()));
		typename math::pcg::Preconditioner<VecValueT>::Ptr precond(
			new PreconditionerType(*laplacian));
		if (!precond->isValid()) {
			precond.reset(new math::pcg::JacobiPreconditioner<LaplacianMatrix>(*laplacian));
		}

		state = math::pcg::solve(*laplacian, *b, *x, *precond, interrupter, state);

		// 6. Populate the output tree with values from the solution vector.
		/// @todo if (state.success) ... ?
		return createTreeFromVector<TreeValueT>(*x, *idxTree, /*background=*/zeroVal<TreeValueT>());
	}

	
} // namespace poisson
} // namespace tools
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

