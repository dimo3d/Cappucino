#pragma once

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
#include <SOP/SOP_Node.h>
#include <UT/UT_Interrupt.h>
#include <SOP_NodeVDB.h>
#include <GU/GU_PrimVDB.h>
#include <openvdb/openvdb.h>
#include <openvdb/Types.h>

#include <openvdb/tools/GridOperators.h>
#include <openvdb/math/ConjGradient.h> // for JacobiPreconditioner

#include <openvdb/tools/LevelSetUtil.h> // for tools::sdfInteriorMask()
#include <openvdb/tools/PoissonSolver.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Prune.h>
#include <openvdb/tools/Composite.h>

#include <ParmFactory.h>
#include <Utils.h>
#include <Diverge.h>

typedef openvdb::BoolGrid   ColliderMaskGrid; ///< @todo really should derive from velocity grid
typedef openvdb::BBoxd      ColliderBBox;
typedef openvdb::Coord      Coord;
typedef openvdb::GridBase           Grid;
typedef openvdb::GridBase::Ptr      GridPtr;
typedef openvdb::GridBase::ConstPtr GridCPtr;
typedef openvdb::GridBase&          GridRef;
typedef const openvdb::GridBase&    GridCRef;




namespace VdbCappucino {
	class SOP_VdbCrossProduct: public openvdb_houdini::SOP_NodeVDB
	{
	public:
		// node contructor for HDK
		static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

		// parameter array for Houdini UI
		static PRM_Template myTemplateList[];

	protected:
		// constructor, destructor
		SOP_VdbCrossProduct(OP_Network *net, const char *name, OP_Operator *op);

		virtual ~SOP_VdbCrossProduct();

		// labeling node inputs in Houdini UI
		virtual const char *inputLabel(unsigned idx) const;

		// main function that does geometry processing
		virtual OP_ERROR cookMySop(OP_Context &context);


	

	};


}
