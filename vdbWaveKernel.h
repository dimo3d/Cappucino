#pragma once
#include <SOP/SOP_Node.h>

namespace VdbCappucino {
	class SOP_VdbWaveKernel : public SOP_Node
	{
	public:
		// node contructor for HDK
		static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

		// parameter array for Houdini UI
		static PRM_Template myTemplateList[];

	protected:
		// constructor, destructor
		SOP_VdbWaveKernel(OP_Network *net, const char *name, OP_Operator *op);

		virtual ~SOP_VdbWaveKernel();

		// labeling node inputs in Houdini UI
		virtual const char *inputLabel(unsigned idx) const;

		// main function that does geometry processing
		virtual OP_ERROR cookMySop(OP_Context &context);

	private:
		// helper function for returning value of parameter
		int DEBUG() { return evalInt("debug", 0, 0); }
		fpreal DIM(int t) { return evalInt("dim", 0, t); }
		fpreal DK(fpreal t) { return evalFloat("dk", 0, t); }
		fpreal ENDK(fpreal t) { return evalFloat("endk", 0, t); }
		fpreal SIGMA(fpreal t) { return evalFloat("sigma", 0, t); }

	};


}
