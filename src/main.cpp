//
//  main.cpp
//  AdsorptionToolBox
//
//  Created by aladshaw3 on 2/23/15.
//  Copyright (c) 2015 Georgia Institute of Technology. All rights reserved.
//

#include <fstream>
#include "flock.h"
#include "school.h"
#include "sandbox.h"
#include "yaml_tests.h"
#include "Trajectory.h"

int main(int argc, const char * argv[])
{	
	int success = 0;
	
	//------------------------------Command Line Interface------------------------
	
	std::cout << argv[0] << std::endl;	//Name of executable with path
	if (argc > 1)						//Next array of arguments
		std::cout << argv[1] << std::endl;
	if (argc > 2)						//Next array of arguments
		std::cout << argv[2] << std::endl;
	
	//Assume argv[2] is a file to open
	if (argc > 2)
	{
		// We assume argv[2] is a filename to open
		std::ifstream the_file ( argv[2] );
		// Always check to see if file opening succeeded
		if ( !the_file.is_open() )
			std::cout<<"Could not open file\n";
		else
		{
			char x;
			// the_file.get ( x ) returns false if the end of the file
			//  is reached or an error occurs
			while ( the_file.get ( x ) )
				std::cout << x << std::endl;
			
			//Note this will get the file even if said file is in another location 
		}
	}
	
	//------------------------------Scenario Suite--------------------------------
	
	//success =  MAGPIE_SCENARIOS("CO2_H2S_C3H8_Input.txt","Talu_All.txt");
  	//success =  MAGPIE_SCENARIOS("FMOFZn_Kr_Xe_MAGPIE_Input.txt","FMOFZn_Kr_Xe_MAGPIE_Scenario.txt");
  	//success =  MAGPIE_SCENARIOS("FMOFCu_Kr_Xe_Neg_Input.txt","FMOFCu_Kr_Xe_Neg_Scenario.txt");
	//success =  MAGPIE_SCENARIOS("FMOFCu_Kr_Xe_Pos_Input.txt","FMOFCu_Kr_Xe_Pos_Scenario.txt");
	//success = gsta_optimize("SU_H2O_AgZ_Input.txt");
	//success = gsta_optimize("SU_H2O_Ag0Z_Input.txt");
	//success = gsta_optimize("SU_H2O_AgZ_Combined.txt");
	//success = gsta_optimize("SU_H2O_AgZ_Combined_Int.txt");
	//success = gsta_optimize("SU_H2O_AgZ_ReducedOnly.txt");
	//success = MAGPIE_SCENARIOS("SU_H2O_AgZ_Ext_Param.txt", "SU_H2O_AgZ_Ext_Scene.txt");
	//success = MAGPIE_SCENARIOS("SU_H2O_AgZ_Int_Param.txt", "SU_H2O_AgZ_Int_Scene.txt");
	//success = MAGPIE_SCENARIOS("SU_H2O_AgZ_All_Param.txt", "SU_H2O_AgZ_All_Scene.txt");
	//success = MAGPIE_SCENARIOS("SU_H2O_AgZ_SUOnly_Param.txt", "SU_H2O_AgZ_SUOnly_Scene.txt");
	//success = SKUA_SCENARIOS("TestScenario.txt", "TestAdsorbent.txt", "TestComponents.txt", "TestAdsorbate.txt");
	//success = SKUA_OPTIMIZE("H2O_MS3A_OPTSCENE.txt", "H2O_MS3A_SORBENT_SKUA.txt", "H2O_MS3A_COMPONENTS.txt", "H2O_MS3A_SORBATE.txt", "H2O_MS3A_DATA.txt");
	//success = SCOPSOWL_OPTIMIZE("H2O_MS3A_OPTSCENE.txt", "H2O_MS3A_SORBENT.txt", "H2O_MS3A_COMPONENTS.txt", "H2O_MS3A_SORBATE.txt", "H2O_MS3A_DATA.txt");
    
	//------------------------------Testing Suite--------------------------------
	
    //success = MACAW_TESTS();
	//success = LARK_TESTS();
	//success = FINCH_TESTS();
	//success = MESH_TESTS();
	//success = SKUA_TESTS();
	//success = EGRET_TESTS();
	//success = SCOPSOWL_TESTS();
	//success = EEL_TESTS();
	//success = MOLA_TESTS();
	//success = DOGFISH_TESTS();
	//success = MONKFISH_TESTS();
	//success = RUN_SANDBOX();
	//success = SHARK_TESTS();
	//success = YAML_TEST01();
	//success = YAML_TEST02();
	//success = YAML_TEST03();
	//success = YAML_WRAPPER_TESTS();
	success = Run_Trajectory();
	
	
	std::cout << "Exit Code:\t" << success << std::endl;
	return success;
}

