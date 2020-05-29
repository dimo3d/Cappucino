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
	

	template<typename BoolTreeType, typename BoundaryOp, typename SurfaceNormalTreeType>
	inline LaplacianMatrix::Ptr
		createISLaplacianWithBoundaryConditions2D(
			const typename BoolTreeType::template ValueConverter<VIndex>::Type& idxTree,
			const BoolTreeType& interiorMask,
			const BoundaryOp& boundaryOp,
			const SurfaceNormalTreeType& surfaceNormalTree,
			typename math::pcg::Vector<LaplacianMatrix::ValueType>& source,
			bool staggered)
	{
		using VIdxTreeT = typename BoolTreeType::template ValueConverter<VIndex>::Type;
		using VIdxLeafMgrT = typename tree::LeafManager<const VIdxTreeT>;

		// The number of active voxels is the number of degrees of freedom.
		const Index64 numDoF = idxTree.activeVoxelCount();

		// Construct the matrix.
		LaplacianMatrix::Ptr laplacianPtr(
			new LaplacianMatrix(static_cast<math::pcg::SizeType>(numDoF)));
		LaplacianMatrix& laplacian = *laplacianPtr;

		// Populate the matrix using a second-order, 7-point CD stencil.
		VIdxLeafMgrT idxLeafManager(idxTree);
		if (staggered) {
			idxLeafManager.foreach(internal::ISStaggeredLaplacianOp<BoolTreeType, BoundaryOp>(
				laplacian, idxTree, interiorMask, boundaryOp, source));
		}
		else {
			idxLeafManager.foreach(internal::ISLaplacianOp2D<VIdxTreeT, BoundaryOp, SurfaceNormalTreeType>(
				laplacian, idxTree, boundaryOp, surfaceNormalTree, source));
		}

		return laplacianPtr;
	}
	namespace internal {
		/// Functor for use with LeafManager::foreach() to populate a sparse %Laplacian matrix
		template<typename VIdxTreeT, typename BoundaryOp, typename SurfaceNormalTreeType>
		struct ISLaplacianOp2D
		{
			using VIdxLeafT = typename VIdxTreeT::LeafNodeType;
			using ValueT = LaplacianMatrix::ValueType;
			using VectorT = typename math::pcg::Vector<ValueT>;

			LaplacianMatrix* laplacian;
			const VIdxTreeT* idxTree;
			const BoundaryOp boundaryOp;
			const SurfaceNormalTreeType* surfaceNormalTree;
			VectorT* source;

			ISLaplacianOp2D(LaplacianMatrix& m, const VIdxTreeT& idx, const BoundaryOp& op, const SurfaceNormalTreeType& snt, VectorT& src) :
				laplacian(&m), idxTree(&idx), boundaryOp(op), surfaceNormalTree(&snt), source(&src) {}

			void operator()(const VIdxLeafT& idxLeaf, size_t /*leafIdx*/) const
			{
				typename tree::ValueAccessor<const VIdxTreeT> vectorIdx(*idxTree);
				typename tree::ValueAccessor<const SurfaceNormalTreeType> surfNormalIdx(*surfaceNormalTree);
				const int kNumOffsets = 6;
				const Coord ijkOffset[kNumOffsets] = {
	#if OPENVDB_TOOLS_POISSON_LAPLACIAN_STENCIL == 1
					Coord(-1,0,0), Coord(1,0,0), Coord(0,-1,0), Coord(0,1,0), Coord(0,0,-1), Coord(0,0,1)
	#else
					Coord(-2,0,0), Coord(2,0,0), Coord(0,-2,0), Coord(0,2,0), Coord(0,0,-2), Coord(0,0,2)
	#endif
				};

				// For each active voxel in this leaf...
				for (typename VIdxLeafT::ValueOnCIter it = idxLeaf.cbeginValueOn(); it; ++it) {
					assert(it.getValue() > -1);

					const Coord ijk = it.getCoord();
					const math::pcg::SizeType rowNum = static_cast<math::pcg::SizeType>(it.getValue());
					SurfaceNormalTreeType::ValueType normal = surfNormalIdx.getValue(ijk);
					normal.normalize();
					LaplacianMatrix::RowEditor row = laplacian->getRowEditor(rowNum);

					ValueT modifiedDiagonal = 0.f;

					// For each of the neighbors of the voxel at (i,j,k)...
					for (int dir = 0; dir < kNumOffsets; ++dir) {
						const Coord neighbor = ijk + ijkOffset[dir];
						Vec3s offsetDir = Vec3s(ijkOffset[dir].x(), ijkOffset[dir].y(), ijkOffset[dir].z());
						offsetDir.normalize();
						const double length_of_projection = normal.length() > 0 ? (offsetDir-offsetDir.projection(normal)).length():1.0f;
						VIndex column;
						// For collocated vector grids, the central differencing stencil requires
						// access to neighbors at a distance of two voxels in each direction
						// (-x, +x, -y, +y, -z, +z).
#if OPENVDB_TOOLS_POISSON_LAPLACIAN_STENCIL == 1
						const bool ijkIsInterior = (vectorIdx.probeValue(neighbor + ijkOffset[dir], column)
							&& vectorIdx.isValueOn(neighbor));
#else
						const bool ijkIsInterior = vectorIdx.probeValue(neighbor, column);
#endif
						if (ijkIsInterior) {
							// If (i,j,k) is sufficiently far away from the exterior,
							// set its weight to one and adjust the center weight accordingly.
							row.setValue(column, length_of_projection);
							modifiedDiagonal -= length_of_projection;
						}
						else {
							// If (i,j,k) is adjacent to or one voxel away from the exterior,
							// invoke the boundary condition functor.
							boundaryOp(ijk, neighbor, source->at(rowNum), modifiedDiagonal);
						}
					}
					// Set the (possibly modified) weight for the voxel at (i,j,k).
					row.setValue(rowNum, modifiedDiagonal);
				}
			}
		};
	}//namespaces
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
		LaplacianMatrix::Ptr laplacian = createISLaplacianWithBoundaryConditions2D(
			*idxTree, *interiorMask, boundaryOp, surfaceNormalTree, *b, staggered);

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

