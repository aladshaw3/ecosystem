##   @package ControlledTWC_data_processing
#
#    @brief Script to read in all folders of Co-Optima controlled experiments
#
#    @details Python script to read in Co-Optima transient data for
#                for all sets of folders of the controlled experiments
#                performed without the presence of HCs.
#
#    @author Austin Ladshaw
#
#    @date 07/14/2021
#
#    @copyright This software was designed and built at the Oak Ridge National
#                    Laboratory (ORNL) National Transportation Research Center
#                    (NTRC) by Austin Ladshaw for research in the catalytic
#                    reduction of NOx. Copyright (c) 2021, all rights reserved.

from transient_data import TransientData
import os, sys, getopt
import time

## Function to read in a single file, do basic data pre-processing, and return transient object for post-processing
def readCoOptimaFile(run_name):

    # Use the fuel name as the folder to store processed data
    #       Sub-directory in that folder will be named after date and run #
    print("Reading file: " + run_name + "...")

    base_folder = run_name.split("/")[0]
    sub_folder =  run_name.split("/")[0] + "/" + run_name.split("/")[1].split("-")[0] + "_run" + run_name.split("/")[1].split("-")[-1]
    print(sub_folder)
    print(base_folder)
    sub_folder = sub_folder.split("/")[1]

    run = TransientData(run_name)

    # Next, the Ethanol and CO columns have different tolerances, but with different units.
    # To deal with this, we need to first do some column manipulation to create new columns based
    # on the ppm tolerances.
    run.mathOperation('CO% (1)',"*",10000,True,'CO (10000)')

    run.deleteColumns(['CO% (1)'])

    # Now we can compress the columns for each AND delete columns we don't need...
    run.compressColumns()

    #       NOTE: Many of these items to retain may change from file to file, which is why we have to do this manually
    #run.displayColumnNames()
    #exit()

    # CO2_H2O
    if base_folder == "CO2_H2O":
        list = ['Elapsed Time (min)', 'N2O (100,200,300)', 'NO (350,3000)',
                'NO2 (150)', 'H2O% (20)', 'CO2% (20)', 'AI 2', 'TC top sample in (C)',
                'TC top sample mid 2 (C)', 'TC top sample out (C)',
                'P tup in (bar)', 'DP top (bar)', 'P top out (bar)',
                'NH3 (300,3000)', 'CO (500,10000)']

    run.retainOnlyColumns(list)

    # If not errors are reported, then both the bypass and the run will have the same columns/names at this point

    # Next, remove any negative values from observations since they are not physically relavent
    run.removeNegatives(list)

    # Next, append columns to the run using the average inlet values
    if base_folder == "CO2_H2O":
        frame_list = [0.]*len(run.getTimeFrames())
        run.appendColumnByFrame('O2% inlet', frame_list)

        frame_list = [13.]*len(run.getTimeFrames())
        run.appendColumnByFrame('CO2% inlet', frame_list)

        frame_list = [13.]*len(run.getTimeFrames())
        run.appendColumnByFrame('H2O% inlet', frame_list)

        frame_list = [0.]*len(run.getTimeFrames())
        run.appendColumnByFrame('CO ppm inlet', frame_list)

        frame_list = [0.]*len(run.getTimeFrames())
        run.appendColumnByFrame('H2 ppm inlet', frame_list)

        frame_list = [0.]*len(run.getTimeFrames())
        run.appendColumnByFrame('NO ppm inlet', frame_list)


    #run.displayColumnNames()
    #exit()

    H2_name = ""

    if base_folder == "CO2_H2O":
        H2_name = "AI 2"

    max_val = run.getMaximum(H2_name)
    min_val = run.getMinimum(H2_name)
    run.mathOperation(H2_name,"-",min_val,append_new=True,append_name="H2 normal signal")
    run.mathOperation("H2 normal signal","/",max_val)

    #run.displayColumnNames()
    #exit()

    return run

## Function to read in a specific folder
def readCoOptimaPureFuelFolder(folder):
    # Read in the bypass and run files separately
    run = []

    # CO2_H2O (2)
    if folder == "CO2_H2O":
        run_name =    "20210504-CPTMA-MalibuTWC-SGDI-30k-13pctCO2-13pctH2O-5Cramp-lambda0_999-1"
        run.append(readCoOptimaFile(folder+"/"+run_name))
        avg_run = readCoOptimaFile(folder+"/"+run_name)

        run_name =    "20210510-CPTMA-MalibuTWC-SGDI-30k-13pctCO2-13pctH2O-5Cramp-lambda0_999-1"
        run.append(readCoOptimaFile(folder+"/"+run_name))


    #After the nested if statements, the run list should be holding all pre-processed data
    # Loop through this data to produce plots, charts, output, and post-processing
    i=0
    for obj in run:
        #Start aggregating runs to create the average
        if i > 0:
            for item in avg_run.data_map:
                j=0
                for value in avg_run.data_map[item]:
                    try:
                        avg_run.data_map[item][j] += obj.data_map[item][j]
                    except:
                        avg_run.data_map[item][j] += obj.data_map[item][-1]
                    j+=1

        # At this point, we can automatically create and save some plots for visualization
        # NOTE: The time frame indexed by 1 represents the temperature ramp
        base_name = obj.input_file_name.split(".")[0]

        base_folder = base_name.split("-")[5]
        if base_folder == "1":
            base_folder += base_name.split("-")[6]
        sub_folder = base_name.split("-")[0] + "_run" + base_name.split("-")[-1]

        i+=1


    #Output information for the average run
    # First, loop to average the sum
    for item in avg_run.data_map:
        j=0
        for value in avg_run.data_map[item]:
            avg_run.data_map[item][j] = avg_run.data_map[item][j]/len(run)
            j+=1


    #Next create the plots and calculate some specific rate info
    base_name = avg_run.input_file_name.split(".")[0]

    base_folder = base_name.split("-")[5]
    if base_folder == "1":
        base_folder += base_name.split("-")[6]
    sub_folder = base_name.split("-")[0] + "_avg"


    # Lastly, we will compress the rows and print the data to a file
    i=0
    for obj in run:
        base_name = obj.input_file_name.split(".")[0]

        base_folder = base_name.split("-")[5]
        if base_folder == "1":
            base_folder += base_name.split("-")[6]
        sub_folder = base_name.split("-")[0] + "_run" + base_name.split("-")[-1]

        obj.compressRows(10)
        obj.printAlltoFile(base_folder+"-output/"+sub_folder+"/"+base_name+"_output.dat")
        i+=1

    #Output for the averaged data
    base_name = avg_run.input_file_name.split(".")[0]

    base_folder = base_name.split("-")[5]
    if base_folder == "1":
        base_folder += base_name.split("-")[6]
    sub_folder = base_name.split("-")[0] + "_avg"

    avg_run.compressRows(10)
    avg_run.printAlltoFile(base_folder+"-output/"+sub_folder+"/"+base_name+"_Avg_output.dat")
    return


## Main function
def main(argv):
    # These data sets have a bypass run, but we cannot use info from that bypass
    # very effectively, thus, we will just read in the signal responses and
    # condense that information into an easier to digest format.

    readCoOptimaPureFuelFolder("CO2_H2O")

    return

##Directs python to call the main function
if __name__ == "__main__":
   main(sys.argv[1:])
