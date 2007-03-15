// Plugin definition for the algorithm

#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06AlignmentAlgorithm.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmPluginFactory.h"

#include "PluginManager/ModuleDef.h"

DEFINE_SEAL_MODULE();
DEFINE_SEAL_PLUGIN( AlignmentAlgorithmPluginFactory,
					CSA06AlignmentAlgorithm, "CSA06AlignmentAlgorithm" );

