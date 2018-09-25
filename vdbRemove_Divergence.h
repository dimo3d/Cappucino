#pragma once
#include <SOP/SOP_Node.h>
#include <openvdb/openvdb.h>
#include <openvdb/Types.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/math/ConjGradient.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <UT/UT_Interrupt.h>
#include <SOP_NodeVDB.h>
#include <GU/GU_PrimVDB.h>

typedef openvdb::BoolGrid   ColliderMaskGrid; ///< @todo really should derive from velocity grid
typedef openvdb::BBoxd      ColliderBBox;
typedef openvdb::Coord      Coord;
typedef openvdb::GridBase           Grid;
typedef openvdb::GridBase::Ptr      GridPtr;
typedef openvdb::GridBase::ConstPtr GridCPtr;
typedef openvdb::GridBase&          GridRef;
typedef const openvdb::GridBase&    GridCRef;


enum ColliderType { CT_NONE, CT_BBOX, CT_STATIC, CT_DYNAMIC };
/// @brief Wrapper class that adapts a Houdini @c UT_Interrupt object
/// for use with OpenVDB library routines
/// @sa openvdb/util/NullInterrupter.h
class Interrupter
{
public:
	Interrupter(const char* title = NULL) :
		mUTI(UTgetInterrupt()), mRunning(false)
	{
		if (title) mUTI->setAppTitle(title);
	}
	~Interrupter() { if (mRunning) this->end(); }

	/// @brief Signal the start of an interruptible operation.
	/// @param name  an optional descriptive name for the operation
	void start(const char* name = NULL) { if (!mRunning) { mRunning = true; mUTI->opStart(name); } }
	/// Signal the end of an interruptible operation.
	void end() { if (mRunning) { mUTI->opEnd(); mRunning = false; } }

	/// @brief Check if an interruptible operation should be aborted.
	/// @param percent  an optional (when >= 0) percentage indicating
	///     the fraction of the operation that has been completed
	bool wasInterrupted(int percent = -1) { return mUTI->opInterrupt(percent); }

private:
	UT_Interrupt* mUTI;
	bool mRunning;
};


namespace VdbCappucino {
	class SOP_VdbRemove_Divergence : public openvdb_houdini::SOP_NodeVDB
	{
	public:
		// node contructor for HDK
		static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

		// parameter array for Houdini UI
		static PRM_Template myTemplateList[];

	protected:
		// constructor, destructor
		SOP_VdbRemove_Divergence(OP_Network *net, const char *name, OP_Operator *op);

		virtual ~SOP_VdbRemove_Divergence();

		// labeling node inputs in Houdini UI
		virtual const char *inputLabel(unsigned idx) const;

		// main function that does geometry processing
		virtual OP_ERROR cookMySop(OP_Context &context);

	private:
		// helper function for returning value of parameter
		int DEBUG() { return evalInt("debug", 0, 0); }
		float DT() { return evalFloat("dt", 0, 0); }

	};


}
