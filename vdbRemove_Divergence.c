#include "vdbRemove_Divergence.h"
#include <limits.h>
#include <SYS/SYS_Math.h>


#include <UT/UT_Interrupt.h>

#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>

#include <PRM/PRM_Include.h>
#include <CH/CH_LocalVariable.h>

#include <OP/OP_AutoLockInputs.h>

#include <GU/GU_PrimVDB.h>
#include <openvdb/math/ConjGradient.h> // for JacobiPreconditioner
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/LevelSetUtil.h> // for tools::sdfInteriorMask()
#include <openvdb/tools/PoissonSolver.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Prune.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/openvdb.h>
#include <ParmFactory.h>
#include <Utils.h>
using namespace VdbCappucino;
namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;
using ColorAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Color_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using VelocityAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Velocity_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using GradientAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Gradient_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;

// label node inputs, 0 corresponds to first input, 1 to the second one
const char *
SOP_VdbRemove_Divergence::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "velocity";
	case 1: return "gradient of distancefield";
	default: return "external divergence";
	}
};




// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbRemove_Divergence::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbRemove_Divergence(net, name, op);
}

SOP_VdbRemove_Divergence::SOP_VdbRemove_Divergence(OP_Network *net, const char *name, OP_Operator *op) : openvdb_houdini::SOP_NodeVDB(net, name, op) {}

SOP_VdbRemove_Divergence::~SOP_VdbRemove_Divergence() {}



namespace {


	template<typename TreeType>
	struct CorrectVelocityOp
	{
		typedef typename TreeType::LeafNodeType LeafNodeType;
		typedef typename TreeType::ValueType    ValueType;

		CorrectVelocityOp(LeafNodeType** velocityNodes, const LeafNodeType** gradientOfPressureNodes, double dx)
			: mVelocityNodes(velocityNodes), mGradientOfPressureNodes(gradientOfPressureNodes), mVoxelSize(dx)
		{
		}

		void operator()(const tbb::blocked_range<size_t>& range) const {

			typedef typename ValueType::value_type ElementType;
			const ElementType scale = ElementType(mVoxelSize * mVoxelSize);

			for (size_t n = range.begin(), N = range.end(); n < N; ++n) {

				LeafNodeType& velocityNode = *mVelocityNodes[n];
				ValueType* velocityData = velocityNode.buffer().data();

				const ValueType* gradientOfPressureData = mGradientOfPressureNodes[n]->buffer().data();

				for (typename LeafNodeType::ValueOnIter it = velocityNode.beginValueOn(); it; ++it) {
					const openvdb::Index pos = it.pos();
					velocityData[pos] -= scale * gradientOfPressureData[pos];
				}
			}
		}

		LeafNodeType       * const * const mVelocityNodes;
		LeafNodeType const * const * const mGradientOfPressureNodes;
		double                       const mVoxelSize;
	}; // class CorrectVelocityOp


	   /// Constant boundary condition functor
	struct DirichletOp {
		inline void operator()(const openvdb::Coord&,
			const openvdb::Coord&, double&, double& diag) const {
			diag -= 1;
		}
	};


	template<typename VectorGridType>
	inline bool
		removeDivergence(VectorGridType& velocityGrid, openvdb::Vec3SGrid::ConstPtr gradient_grid, openvdb::FloatGrid::ConstPtr external_divergencegrid, const int iterations, hvdb::Interrupter& interrupter)
	{
		typedef typename VectorGridType::TreeType       VectorTreeType;
		typedef typename VectorTreeType::LeafNodeType   VectorLeafNodeType;
		typedef typename VectorGridType::ValueType      VectorType;
		typedef typename VectorType::ValueType          VectorElementType;

		typedef typename VectorGridType::template ValueConverter<VectorElementType>::Type   ScalarGrid;
		typedef typename ScalarGrid::TreeType                                               ScalarTree;

		//openvdb::tools::Divergence<VectorGridType> divergenceOp(velocityGrid);
		//typename ScalarGrid::Ptr divGrid = divergenceOp.process();
		openvdb::VectorGrid::Ptr divGrid = velocityGrid.deepCopy();
		openvdb::tools::foreach(divGrid->beginValueOn(), Diverge(divGrid->transform(), velocityGrid, gradient_grid), false, 0);
		openvdb::FloatGrid::Ptr
			sourceGrid = external_divergencegrid->deepCopy();
		openvdb::FloatGrid::Ptr	targetGrid = divGrid->deepCopy();
		// Get the source and target grids' index space to world space transforms.
		const openvdb::math::Transform
			&sourceXform = sourceGrid->transform(),
			&targetXform = targetGrid->transform();
		// Compute a source grid to target grid transform.
		// (For this example, we assume that both grids' transforms are linear,
		// so that they can be represented as 4 x 4 matrices.)
		openvdb::Mat4R xform =
			sourceXform.baseMap()->getAffineMap()->getMat4() *
			targetXform.baseMap()->getAffineMap()->getMat4().inverse();
		// Create the transformer.
		openvdb::tools::GridTransformer transformer(xform);
		// Resample using trilinear interpolation.
		transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(
			*sourceGrid, *targetGrid);
		
		// Prune the target tree for optimal sparsity.
		targetGrid->tree().prune();

		openvdb::FloatGrid::Ptr diffDivergence = divGrid->deepCopy();
		// Define a local function that subtracts two floating-point values.
		struct Local {
			static inline void diff(const float& a, const float& b, float& result) {
				result = a - b;
			}
		};
		
		// Compute the difference between corresponding voxels of aGrid and bGrid
		// and store the result in aGrid, leaving bGrid empty.
		diffDivergence->tree().combine(targetGrid->tree(), Local::diff);
		openvdb::math::pcg::State state = openvdb::math::pcg::terminationDefaults<VectorElementType>();
		state.iterations = iterations;
		state.relativeError = state.absoluteError = openvdb::math::Delta<VectorElementType>::value();

		typedef openvdb::math::pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix> PCT;

		typename ScalarTree::Ptr pressure =
			openvdb::tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
				diffDivergence->tree(), DirichletOp(), state, interrupter);

		typename ScalarGrid::Ptr pressureGrid = ScalarGrid::create(pressure);
		pressureGrid->setTransform(velocityGrid.transform().copy());

		openvdb::tools::Gradient<ScalarGrid> gradientOp(*pressureGrid);
		typename VectorGridType::Ptr gradientOfPressure = gradientOp.process();

		struct ProjectVectorToSurface {
			openvdb::Vec3SGrid::ConstPtr grad_grid;
			openvdb::math::Transform trans;
			ProjectVectorToSurface(
				openvdb::Vec3SGrid::ConstPtr grad_g, openvdb::math::Transform tr
				) : grad_grid(grad_g), trans(tr)
			{}

			inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
				std::unique_ptr<GradientAccessor> gradientAccessor;
				std::unique_ptr<Gradient_fastSampler> gradient_fastSampler;
								
				gradientAccessor.reset(new GradientAccessor(grad_grid->getConstAccessor()));
				gradient_fastSampler.reset(new Gradient_fastSampler(*gradientAccessor, grad_grid->transform()));
			
				
				openvdb::Vec3f normal = gradient_fastSampler->wsSample(trans.indexToWorld(iter.getCoord()));
				normal.normalize();	
				openvdb::Vec3f projected_Velocity = iter.getValue().projection(normal);
				iter.setValue(iter.getValue()-projected_Velocity);
			}
		};


		// Iterate over all active values.
		openvdb::tools::foreach(gradientOfPressure->beginValueOn(), ProjectVectorToSurface( gradient_grid, gradientOfPressure->transform()), true, 0);

		{
			std::vector<VectorLeafNodeType*> velocityNodes;
			velocityGrid.tree().getNodes(velocityNodes);

			std::vector<const VectorLeafNodeType*> gradientNodes;
			gradientNodes.reserve(velocityNodes.size());
			gradientOfPressure->tree().getNodes(gradientNodes);

			const double dx = velocityGrid.transform().voxelSize()[0];

			tbb::parallel_for(tbb::blocked_range<size_t>(0, velocityNodes.size()),
				CorrectVelocityOp<VectorTreeType>(&velocityNodes[0], &gradientNodes[0], dx));
			
			openvdb::tools::foreach(velocityGrid.beginValueOn(), ProjectVectorToSurface(gradient_grid, velocityGrid.transform()), true, 0);

			
		}

		return state.success;
	}


} // unnamed namespace

OP_ERROR
SOP_VdbRemove_Divergence::cookMySop(OP_Context& context)
{
	try {
		hutil::ScopedInputLock lock(*this, context);
		duplicateSourceStealable(0, context);

		const fpreal time = context.getTime();

		hvdb::Interrupter boss("Removing Divergence");

		const int iterations = evalInt("iterations", 0, time);

		UT_String velocityGroupStr;
		evalString(velocityGroupStr, "velocitygroup", 0, time);
		const GA_PrimitiveGroup* velocityGroup = matchGroup(*gdp, velocityGroupStr.toStdString());

		
		const GU_Detail* gradientGdp = inputGeo(1, context);
		const GU_Detail* divergenceGdp = inputGeo(2, context);

		UT_String gradientGroupStr;
		evalString(gradientGroupStr, "gradientgroup", 0, time);
		const GA_PrimitiveGroup* gradientGroup = matchGroup(const_cast<GU_Detail&>(*gradientGdp), gradientGroupStr.toStdString());

		UT_String divergenceGroupStr;
		evalString(divergenceGroupStr, "divergencegroup", 0, time);
		const GA_PrimitiveGroup* divergenceGroup = matchGroup(const_cast<GU_Detail&>(*divergenceGdp), divergenceGroupStr.toStdString());

		bool processedVDB = false;
		//get gradient
		hvdb::VdbPrimCIterator gIt(gradientGdp, gradientGroup);
		const GU_PrimVDB *gradientPrim = *gIt; 

		if (!gradientPrim) {
			addError(SOP_MESSAGE, "Missing gradient grid");
			return error();
		}
		if (gradientPrim->getStorageType() != UT_VDB_VEC3F) {
			addError(SOP_MESSAGE, "Expected gradient grid to be of type Vec3f");
			return error();
		}
		openvdb::GridBase::ConstPtr gradient_baseGrid;
		openvdb::Vec3SGrid::ConstPtr gradient_grid;
		gradient_baseGrid = gradientPrim->getConstGridPtr();
		gradient_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(gradient_baseGrid);
		if (!gradient_grid) {
			addError(SOP_MESSAGE, "Missing gradient grid");
			return error();
		}
		//get external divergence
		hvdb::VdbPrimCIterator dIt(divergenceGdp, divergenceGroup);
		const GU_PrimVDB *divergencePrim = *dIt;

		if (!divergencePrim) {
			addError(SOP_MESSAGE, "Missing divergence grid");
			return error();
		}
		if (divergencePrim->getStorageType() != UT_VDB_FLOAT) {
			addError(SOP_MESSAGE, "Expected divergence grid to be of type Skalar");
			return error();
		}
		openvdb::GridBase::ConstPtr divergence_baseGrid;
		openvdb::FloatGrid::ConstPtr divergence_grid;
		divergence_baseGrid = divergencePrim->getConstGridPtr();
		divergence_grid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(divergence_baseGrid);
		if (!divergence_grid) {
			addError(SOP_MESSAGE, "Missing divergence grid");
			return error();
		}

		//process
		for (hvdb::VdbPrimIterator vdbIt(gdp, velocityGroup); vdbIt; ++vdbIt) {

			if (boss.wasInterrupted()) break;

			if (vdbIt->getGrid().type() == openvdb::Vec3fGrid::gridType()) {

				processedVDB = true;

				vdbIt->makeGridUnique();

				openvdb::Vec3fGrid& grid = static_cast<openvdb::Vec3fGrid&>(vdbIt->getGrid());

				if (!removeDivergence(grid, gradient_grid, divergence_grid, iterations, boss) && !boss.wasInterrupted()) {
					const std::string msg = grid.getName() + " did not fully converge.";
					addWarning(SOP_MESSAGE, msg.c_str());
				}
			}
		}

		if (!processedVDB && !boss.wasInterrupted()) {
			addWarning(SOP_MESSAGE, "No Vec3f VDBs found.");
		}

	}
	catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
	}

	return error();
}