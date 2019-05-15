## Python script to read cardinal control files and make changes ##
## Run python scripts using Python 3.5 or newer ##

''' Control file script:
    ----------------
    Object-Oriented approach to reading the Cardinal control files
    This object will use the previously established python yaml reader
    to read the control files, then make changes to the python yaml map
    and record those changes in the new control file.
    
    Author:     Austin Ladshaw
    Date:       05/14/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.
    '''
import sys, os
sys.path.insert(0, '../yaml-reader')
import yaml_reader as yr
from enum import Enum

class ControlParam(Enum):
    Yield  = 0
    Height = 1
    Level = 2
    Wind = 3
    Mass = 4
    FissExtent = 5
    FissYield = 6
    PartDist = 7
    Casing = 8
    FracMod = 9

class ControlFile(yr.YAML):
    def __init__(self):
        yr.YAML.__init__(self)

    def printMap(self):
        print(self.map)

    def editValue_PercentChange(self, e_num, per):
        if e_num == ControlParam.Yield:
            self.map["Simulation_Conditions"]["bomb_yield"] = self.map["Simulation_Conditions"]["bomb_yield"]*(1.0 + (per/100.0))
            if self.map["Simulation_Conditions"]["bomb_yield"] < 1.0:
                self.map["Simulation_Conditions"]["bomb_yield"] = 1.0
        if e_num == ControlParam.Height:
            if self.map["Simulation_Conditions"]["burst_height"] == 0.0:
                self.map["Simulation_Conditions"]["burst_height"] = per
            else:
                self.map["Simulation_Conditions"]["burst_height"] = self.map["Simulation_Conditions"]["burst_height"]*(1.0 + (per/100.0))
        if e_num == ControlParam.Level:
            if self.map["Simulation_Conditions"]["ground_level"] == 0.0:
                self.map["Simulation_Conditions"]["ground_level"] = per
            else:
                self.map["Simulation_Conditions"]["ground_level"] = self.map["Simulation_Conditions"]["ground_level"]*(1.0 + (per/100.0))
        if e_num == ControlParam.Wind:
            for alt in self.map["Wind_Profile"]:
                self.map["Wind_Profile"][alt]["vx"] = self.map["Wind_Profile"][alt]["vx"]*(1.0 + (per/100.0))
                self.map["Wind_Profile"][alt]["vy"] = self.map["Wind_Profile"][alt]["vy"]*(1.0 + (per/100.0))
        if e_num == ControlParam.Mass:
            self.map["Weapon_Info"]["total_mass"] = self.map["Weapon_Info"]["total_mass"]*(1.0 + (per/100.0))
            if self.map["Weapon_Info"]["total_mass"] < 1.0:
                self.map["Weapon_Info"]["total_mass"] = 1.0
        if e_num == ControlParam.FissExtent:
            self.map["Weapon_Info"]["fission_extent"] = self.map["Weapon_Info"]["fission_extent"]*(1.0 + (per/100.0))
            if self.map["Weapon_Info"]["fission_extent"] > 100:
                self.map["Weapon_Info"]["fission_extent"] = 100.0
            if self.map["Weapon_Info"]["fission_extent"] < 1:
                self.map["Weapon_Info"]["fission_extent"] = 1.0
        if e_num == ControlParam.FissYield:
            self.map["Weapon_Info"]["fission_yield"] = self.map["Weapon_Info"]["fission_yield"]*(1.0 + (per/100.0))
            if self.map["Weapon_Info"]["fission_yield"] > self.map["Simulation_Conditions"]["bomb_yield"]:
                self.map["Weapon_Info"]["fission_yield"] = self.map["Simulation_Conditions"]["bomb_yield"]
            if self.map["Weapon_Info"]["fission_yield"] < 0.0:
                self.map["Weapon_Info"]["fission_yield"] = 1.0
        if e_num == ControlParam.PartDist:
            if self.map["Simulation_Conditions"]["part_dist"]["useCustom"] == True:
                self.map["Simulation_Conditions"]["part_dist"]["min_dia"] = self.map["Simulation_Conditions"]["part_dist"]["min_dia"]*(1.0 + (per/100.0))
                self.map["Simulation_Conditions"]["part_dist"]["max_dia"] = self.map["Simulation_Conditions"]["part_dist"]["max_dia"]*(1.0 + (per/100.0))
                self.map["Simulation_Conditions"]["part_dist"]["mean_dia"] = self.map["Simulation_Conditions"]["part_dist"]["mean_dia"]*(1.0 + (per/100.0))
                self.map["Simulation_Conditions"]["part_dist"]["std_dia"] = self.map["Simulation_Conditions"]["part_dist"]["std_dia"]*(1.0 + (per/100.0))
        if e_num == ControlParam.Casing:
            self.map["Weapon_Info"]["casing_thickness"] = self.map["Weapon_Info"]["casing_thickness"]*(1.0 + (per/100.0))
            if self.map["Weapon_Info"]["casing_thickness"] < 1.0:
                self.map["Weapon_Info"]["casing_thickness"] = 1.0
        if e_num == ControlParam.FracMod:
            if self.map["Weapon_Info"]["fractionation_model"] == "modified-freiling-tompkins":
                self.map["Weapon_Info"]["fractionation_model"] = "freiling-tompkins"
            elif self.map["Weapon_Info"]["fractionation_model"] == "freiling-tompkins":
                self.map["Weapon_Info"]["fractionation_model"] = "modified-freiling-tompkins"
            if self.map["Weapon_Info"]["fractionation_model"] == "freiling":
                self.map["Weapon_Info"]["fractionation_model"] = "modified-freiling"
            elif self.map["Weapon_Info"]["fractionation_model"] == "modified-freiling":
                self.map["Weapon_Info"]["fractionation_model"] = "freiling"

    def editValue_LinearChange(self, e_num, value):
        if e_num == ControlParam.Yield:
            self.map["Simulation_Conditions"]["bomb_yield"] = self.map["Simulation_Conditions"]["bomb_yield"] + value
            if self.map["Simulation_Conditions"]["bomb_yield"] < 1.0:
                self.map["Simulation_Conditions"]["bomb_yield"] = 1.0
        if e_num == ControlParam.Height:
            if self.map["Simulation_Conditions"]["burst_height"] == 0.0:
                self.map["Simulation_Conditions"]["burst_height"] = value
            else:
                self.map["Simulation_Conditions"]["burst_height"] = self.map["Simulation_Conditions"]["burst_height"] + value
        if e_num == ControlParam.Level:
            if self.map["Simulation_Conditions"]["ground_level"] == 0.0:
                self.map["Simulation_Conditions"]["ground_level"] = value
            else:
                self.map["Simulation_Conditions"]["ground_level"] = self.map["Simulation_Conditions"]["ground_level"] + value
        if e_num == ControlParam.Wind:
            for alt in self.map["Wind_Profile"]:
                self.map["Wind_Profile"][alt]["vx"] = self.map["Wind_Profile"][alt]["vx"] + value
                self.map["Wind_Profile"][alt]["vy"] = self.map["Wind_Profile"][alt]["vy"] + value
        if e_num == ControlParam.Mass:
            self.map["Weapon_Info"]["total_mass"] = self.map["Weapon_Info"]["total_mass"] + value
            if self.map["Weapon_Info"]["total_mass"] < 1.0:
                self.map["Weapon_Info"]["total_mass"] = 1.0
        if e_num == ControlParam.FissExtent:
            self.map["Weapon_Info"]["fission_extent"] = self.map["Weapon_Info"]["fission_extent"] + value
            if self.map["Weapon_Info"]["fission_extent"] > 100:
                self.map["Weapon_Info"]["fission_extent"] = 100.0
            if self.map["Weapon_Info"]["fission_extent"] < 1:
                self.map["Weapon_Info"]["fission_extent"] = 1.0
        if e_num == ControlParam.FissYield:
            self.map["Weapon_Info"]["fission_yield"] = self.map["Weapon_Info"]["fission_yield"] + value
            if self.map["Weapon_Info"]["fission_yield"] > self.map["Simulation_Conditions"]["bomb_yield"]:
                self.map["Weapon_Info"]["fission_yield"] = self.map["Simulation_Conditions"]["bomb_yield"]
            if self.map["Weapon_Info"]["fission_yield"] < 0.0:
                self.map["Weapon_Info"]["fission_yield"] = 1.0
        if e_num == ControlParam.PartDist:
            if self.map["Simulation_Conditions"]["part_dist"]["useCustom"] == True:
                self.map["Simulation_Conditions"]["part_dist"]["min_dia"] = self.map["Simulation_Conditions"]["part_dist"]["min_dia"] + value
                self.map["Simulation_Conditions"]["part_dist"]["max_dia"] = self.map["Simulation_Conditions"]["part_dist"]["max_dia"] + value
                self.map["Simulation_Conditions"]["part_dist"]["mean_dia"] = self.map["Simulation_Conditions"]["part_dist"]["mean_dia"] + value
                self.map["Simulation_Conditions"]["part_dist"]["std_dia"] = self.map["Simulation_Conditions"]["part_dist"]["std_dia"] + value
        if e_num == ControlParam.Casing:
            self.map["Weapon_Info"]["casing_thickness"] = self.map["Weapon_Info"]["casing_thickness"] + value
            if self.map["Weapon_Info"]["casing_thickness"] < 1.0:
                self.map["Weapon_Info"]["casing_thickness"] = 1.0
        if e_num == ControlParam.FracMod:
            if self.map["Weapon_Info"]["fractionation_model"] == "modified-freiling-tompkins":
                self.map["Weapon_Info"]["fractionation_model"] = "freiling-tompkins"
            elif self.map["Weapon_Info"]["fractionation_model"] == "freiling-tompkins":
                self.map["Weapon_Info"]["fractionation_model"] = "modified-freiling-tompkins"
            if self.map["Weapon_Info"]["fractionation_model"] == "freiling":
                self.map["Weapon_Info"]["fractionation_model"] = "modified-freiling"
            elif self.map["Weapon_Info"]["fractionation_model"] == "modified-freiling":
                self.map["Weapon_Info"]["fractionation_model"] = "freiling"

    def clearMap(self):
        self.map.clear()

## Testing ##
'''
control = ControlFile()
control.readFile("1979-Test-Case.txt")
control.editValue_LinearChange(ControlParam.Mass, 10.0)
control.print2file("1979-Test-Case-Mod.txt")
'''

