#include "vdbWave.h"
#include "vdbWaveKernel.h"
#include "vdbCpt.h"
#include "vdbConvolve.h"
#include "vdbReact.h"
#include "vdbDivergence.h"
#include "vdbRemove_Divergence.h"
#include "vdbProjectVector.h"
#include <limits.h>
#include <SYS/SYS_Math.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>

#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>

#include <PRM/PRM_Include.h>
#include <CH/CH_LocalVariable.h>

#include <OP/OP_AutoLockInputs.h>

#include <GU/GU_PrimVDB.h>

#include <openvdb/openvdb.h>
#include <ParmFactory.h>
namespace hutil = houdini_utils;

using namespace VdbCappucino;

// register the operator in Houdini, it is a hook ref for Houdini
void newSopOperator(OP_OperatorTable *table)
{
    OP_Operator *op_wave;
	OP_Operator *op_waveKernel;
	OP_Operator *op_cpt;
	OP_Operator *op_convolve;
	OP_Operator *op_react;
	OP_Operator *op_divergence;
	OP_Operator *op_remove_divergence;
	OP_Operator *op_projectVectorToSurface;
	op_wave = new OP_Operator(
    		"vdbWave",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
    		"VDB Wave",                   // UI name
    		SOP_VdbWave::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
    		SOP_VdbWave::myTemplateList,    // my parameters - An array of PRM_Template objects defining the parameters to this operator
    		2,                                            // min # of sources
    		2);                                           // max # of sources

    // place this operator under the VDB submenu in the TAB menu.
    op_wave->setOpTabSubMenuPath("VDB");

    // after addOperator(), 'table' will take ownership of 'op'
    table->addOperator(op_wave);
	
	op_waveKernel = new OP_Operator(
		"vdbWaveKernel",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB Wave Kernel",                   // UI name
		SOP_VdbWaveKernel::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		SOP_VdbWaveKernel::myTemplateList,    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		1,                                            // min # of sources
		1);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_waveKernel->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_waveKernel);

	hutil::ParmList parms_cpt;
	parms_cpt.add(hutil::ParmFactory(PRM_STRING, "group", "Group")
		.setHelpText("Specify grids to process")
		.setChoiceList(&hutil::PrimGroupMenuInput1));
	
	parms_cpt.add(hutil::ParmFactory(PRM_STRING, "cptGroup", "CPTGroup")
		.setHelpText("Specify CPT Grid")
		.setChoiceList(&hutil::PrimGroupMenuInput2));

	parms_cpt.add(hutil::ParmFactory(PRM_STRING, "distanceGroup", "DistanceGroup")
		.setHelpText("Specify Distance Grid")
		.setChoiceList(&hutil::PrimGroupMenuInput3));

	parms_cpt.add(hutil::ParmFactory(PRM_TOGGLE, "debug", "Debug")
		.setDefault(PRMzeroDefaults)
		.setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 1));

	parms_cpt.add(hutil::ParmFactory(PRM_INT_J, "maxcells", "maxcells")
		.setDefault(2)
		.setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 50));

	parms_cpt.add(hutil::ParmFactory(PRM_TOGGLE, "doworldpos", "doworldpos")
		.setDefault(PRMzeroDefaults)
		.setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 1));

	op_cpt = new OP_Operator(
		"vdbcpt",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB CPT",                   // UI name
		SOP_VdbCpt::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		parms_cpt.get(),    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		3,                                            // min # of sources
		3);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_cpt->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_cpt);

	op_convolve = new OP_Operator(
		"vdbconvolve",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB Convolve",                   // UI name
		SOP_VdbConvolve::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		SOP_VdbConvolve::myTemplateList,    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		2,                                            // min # of sources
		2);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_convolve->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_convolve);
	
	op_react= new OP_Operator(
		"vdbreact",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB React",                   // UI name
		SOP_VdbReact::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		SOP_VdbReact::myTemplateList,    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		2,                                            // min # of sources
		2);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_react->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_react);

	hutil::ParmList parms_divergence;
	parms_divergence.add(hutil::ParmFactory(PRM_STRING, "velocitygroup", "Velocity Group")
		.setHelpText("Specify velocity vector grids to process")
		.setChoiceList(&hutil::PrimGroupMenuInput1));

	parms_divergence.add(hutil::ParmFactory(PRM_STRING, "gradientgroup", "Gradient Group")
		.setHelpText("Specify gradient of distance field")
		.setChoiceList(&hutil::PrimGroupMenuInput2));



	op_divergence = new OP_Operator(
		"vdbdivergence",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB Divergence",                   // UI name
		SOP_VdbDivergence::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		parms_divergence.get(),    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		2,                                            // min # of sources
		2);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_divergence->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_divergence);
	
	
	hutil::ParmList parms;
	parms.add(hutil::ParmFactory(PRM_STRING, "velocitygroup", "Velocity Group")
		.setHelpText("Specify velocity vector grids to process")
		.setChoiceList(&hutil::PrimGroupMenuInput1));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gradientgroup", "Gradient Group")
		.setHelpText("Specify gradient of distance field")
		.setChoiceList(&hutil::PrimGroupMenuInput2));

	parms.add(hutil::ParmFactory(PRM_STRING, "divergencegroup", "Divergence Group")
		.setHelpText("Specify external divergence grids")
		.setChoiceList(&hutil::PrimGroupMenuInput3));

	parms.add(hutil::ParmFactory(PRM_INT_J, "iterations", "Iterations")
		.setDefault(50)
		.setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 100));




	op_remove_divergence = new OP_Operator(
		"vdbremove_divergence",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB Remove Divergence",                   // UI name
		SOP_VdbRemove_Divergence::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		parms.get(),    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		3,                                            // min # of sources
		3);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_remove_divergence->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_remove_divergence);
	

	hutil::ParmList parms_projectVectorToSurface;
	parms_projectVectorToSurface.add(hutil::ParmFactory(PRM_STRING, "velocitygroup", "Velocity Group")
		.setHelpText("Specify velocity vector grids to process")
		.setChoiceList(&hutil::PrimGroupMenuInput1));

	parms_projectVectorToSurface.add(hutil::ParmFactory(PRM_STRING, "gradientgroup", "Gradient Group")
		.setHelpText("Specify gradient of distance field")
		.setChoiceList(&hutil::PrimGroupMenuInput2));



	op_projectVectorToSurface = new OP_Operator(
		"vdbprojectVector",                      // internal name, needs to be unique in OP_OperatorTable (table containing all nodes for a network type - SOPs in our case, each entry in the table is an object of class OP_Operator which basically defines everything Houdini requires in order to create nodes of the new type)
		"VDB Project Vector",                   // UI name
		SOP_VdbProjectVector::myConstructor,     // how to build the node - A class factory function which constructs nodes of this type
		parms_projectVectorToSurface.get(),    // my parameters - An array of PRM_Template objects defining the parameters to this operator
		2,                                            // min # of sources
		2);                                           // max # of sources

													  // place this operator under the VDB submenu in the TAB menu.
	op_projectVectorToSurface->setOpTabSubMenuPath("VDB");

	// after addOperator(), 'table' will take ownership of 'op'
	table->addOperator(op_projectVectorToSurface);

}
