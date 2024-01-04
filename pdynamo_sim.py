import sys                                                                                
sys.path.append("/home/igorchem/pDynamo3_scripts")

import os,glob

from pBabel             import *                                     
from pCore              import *                                     
from pMolecule          import *                    
from SimulationProject import SimulationProject
from TrajectoryAnalysis import *
from ReactionCoordinate import * 

#------------------------------------------
local = os.getcwd()
#------------------------------------------

def Def_MM_Sys():
	'''
	'''
	proj = SimulationProject.From_Force_Field("sys01amber.top","sys01amber.crd",os.path.join( local,"MM_SetUp") )
	parameters_a = {"simulation_type":"Geometry_Optimization","maxIterations":10000,"rmsGradient":0.1 }	
	proj.Run_Simulation(parameters_a)
	proj.SaveSystem()

#================================
# 414+438+79+40+61+628+529+447+715+572+397+330+256
def QCMM_OPT(method="am1"):
	'''
	'''
	
	proj = SimulationProject.From_PKL(os.path.join(local,"MM_SetUp","sys01amber.pkl"),os.path.join(local,"OPT_QMMM",method) )
	
	co2_a  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:*")
	co2_b  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.438:*")
	cat1 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.79:*")
	cat2 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.40:*")
	cat3 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.61:*")
	wat1 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:*")
	wat2 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.529:*")
	wat3 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.447:*")
	wat4 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:*")
	wat5 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:*")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:*")
	aco2 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.330:*")
	aco3 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.256:*")
	
	selections = [ co2_a,co2_b,cat1,cat2,cat3,wat1,wat2,wat3,wat4,wat5,aco,aco2,aco3 ]
	_parameters_b = {"method_class":"SMO","Hamiltonian":"am1","QCcharge":0,"multiplicity":1,"region":selections}
	proj.DEBUG = True
	proj.Set_QC_Method(_parameters_b)

	'''
	parameters = {"maxIterations":2500 						,
				"rmsGradient":0.5   						,
				"log_frequency":10 							,
				"simulation_type":"Geometry_Optimization"	,
				"save_frequency" : 20 						,
				"save_format":".dcd"						,
				"trajectory_name":"opt.ptGeo"				,
				"Debug":True								,
				"save_pdb": True 							}

	proj.Run_Simulation(parameters)
	
	_parameters_b = {"method_class":"SMO","Hamiltonian":method,"QCcharge":0,"multiplicity":1,"region":selections}
	proj.DEBUG = True
	proj.Set_QC_Method(_parameters_b)
	
	parameters = {"maxIterations":2500 						,
				"rmsGradient":0.1   						,
				"log_frequency":10 							,
				"simulation_type":"Geometry_Optimization"	,
				"save_frequency" : 10 						,
				"save_format":".dcd"						,
				"trajectory_name":"opt_Steep.ptGeo"			,
				"Debug":True								,
				"save_pdb": True							}
	proj.Run_Simulation(parameters)
	
	'''	
	xsi = len( glob.glob( os.path.join( local,"OPT_QMMM",method,"opt_Steep.ptGeo","*.pkl") ) )
	_path = os.path.join(local,"OPT_QMMM",method,"opt_Steep.ptGeo")	
	trajAn = TrajectoryAnalysis(_path,proj.system,xsi)
	trajAn.CalculateRG_RMSD(qc_mm=True)
	trajAn.PlotRG_RMS()
	
	

#=======================================================================
def MD_runs(method):
	'''
	'''
	base_pkl = os.path.join(local,"OPT_QMMM","sys01amber.pkl")
	fst_path = os.path.join(local,"MD_QMMM",method,"system_heated.pkl")
	scn_path = os.path.join(local,"MD_QMMM",method,"sys_equilibrated.pkl")
	trd_path = os.path.join(local,"MD_QMMM",method,"production.ptGeo")
	
	if not os.path.exists( os.path.join(local,"MD_QMMM",method) ): 
		os.makedirs(os.path.join(local,"MD_QMMM",method))
	
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"MD_QMMM",method) )
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":method,"QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	proj.SaveSystem(_cname="sys_QC_saved")
	if not os.path.exists(fst_path):
		#----------------------------------
		parameters = {"protocol":"heating"                    ,
						"nsteps":10000                        ,
						"MD_method":"Verlet",
						"temperature_scale_option":"linear"   ,
						"temperature_scale":  		     15   ,
						"start_temperature":             20   ,
						"sampling_factor":               5000 ,
						"log_frequency":				 50   ,
						"simulation_type":"Molecular_Dynamics",
						"temperature":               298.15   }
		proj.Run_Simulation(parameters)
		#----------------------------------
		proj.SaveSystem(_cname="system_heated")
	if not os.path.exists(scn_path):
		_path = os.path.join( local, "MD_QMMM", method, "system_heated.pkl" ) 
		proj = SimulationProject.From_PKL(_path,os.path.join(local,"MD_QMMM",method) )
		#----------------------------------
		parameters = {"protocol":"sampling"                   ,
						"nsteps":50000                        ,
						"MD_method":"LeapFrog"                ,
						"pressure":20.0                       ,
						"pressure_coupling":True              ,
						"sampling_factor":               5000 ,
						"trajectory_name":"equilibration"     ,
						"log_frequency":				   10 ,
						"simulation_type":"Molecular_Dynamics",
						"temperature":               298.15   }
		proj.Run_Simulation(parameters)
		proj.SaveSystem(_cname="sys_equilibrated")

	#------------------------------------
	_path = os.path.join( local, "MD_QMMM", method, "sys_equilibrated.pkl" ) 
	proj = SimulationProject.From_PKL(_path,os.path.join(local,"MD_QMMM",method) )
	parameters = {"protocol":"sampling"                   ,
					"nsteps":200000                       ,
					"MD_method":"LeapFrog"                ,
					"pressure":20.0                       ,
					"pressure_coupling":True              ,
					"sampling_factor":               1000 ,
					"log_frequency":				 20   ,
					"save_format":".dcd.",
					"trajectory_name":"production"        ,
					"simulation_type":"Molecular_Dynamics",
					"temperature":               298.15   }
	proj.Run_Simulation(parameters)
	#------------------------------------
	proj.SaveSystem(_cname="sys_production")


def Analysis(method="",n=0):
	'''
	'''
	proj = SimulationProject.From_PKL(os.path.join(local,"MD_QMMM",method,"system_heated.pkl"),os.path.join(local,"MD_QMMM") )

	C_414  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	C_438  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.438:C")
	OW_628 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	OW_715 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:O")
	OW_572 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:O")
	OW_447 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.447:O")
	HW_628 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	HW_715 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:H1")
	HW_572 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:H1")
	HW_447 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.447:H1")
	ACO_397 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	ACO_256 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.256:O1")
	ACO_256b = AtomSelection.FromAtomPattern(proj.system,"*:ACO.256:O2")
	ACO_330 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.330:O1")

	atoms_a = [HW_628, ACO_397]
	atoms_b = [C_414, OW_628]
	
	atoms_c = [HW_715, ACO_397]
	atoms_d = [C_438, OW_715]
	
	atoms_e = [HW_628, ACO_256]
	atoms_f = [C_438, OW_628]
	
	atoms_g = [HW_715, ACO_256]
	atoms_h = [C_414, OW_715]
	
	atoms_i = [HW_572, ACO_256]
	atoms_j = [C_414, OW_572]
	
	atoms_k = [HW_572, ACO_256]
	atoms_l = [C_438, OW_572]
	
	atoms_m = [HW_447, ACO_256b]
	atoms_n = [OW_447, HW_715]
	
	
	rc1 = ReactionCoordinate(atoms_a,False)	
	rc2 = ReactionCoordinate(atoms_b,False)
	rc3 = ReactionCoordinate(atoms_c,False)
	rc4 = ReactionCoordinate(atoms_d,False)
	rc5 = ReactionCoordinate(atoms_e,False)
	rc6 = ReactionCoordinate(atoms_f,False)
	rc7 = ReactionCoordinate(atoms_g,False)
	rc8 = ReactionCoordinate(atoms_h,False)
	rc9 = ReactionCoordinate(atoms_i,False)
	rc10 = ReactionCoordinate(atoms_j,False)
	rc11 = ReactionCoordinate(atoms_k,False)
	rc12 = ReactionCoordinate(atoms_l,False)
	
	rc13 = ReactionCoordinate(atoms_m,False)
	rc14 = ReactionCoordinate(atoms_n,False)
	
	
	rc1.GetRCLabel(proj.system)	
	rc2.GetRCLabel(proj.system)
	rc3.GetRCLabel(proj.system)
	rc4.GetRCLabel(proj.system)
	rc5.GetRCLabel(proj.system)
	rc6.GetRCLabel(proj.system)
	rc7.GetRCLabel(proj.system)
	rc8.GetRCLabel(proj.system)
	rc9.GetRCLabel(proj.system)
	rc10.GetRCLabel(proj.system)
	rc11.GetRCLabel(proj.system)
	rc12.GetRCLabel(proj.system)
	rc13.GetRCLabel(proj.system)
	rc14.GetRCLabel(proj.system)
	
	rcs_ = [ rc13, rc6 ]
	rcs_a = [ rc1, rc2 ]
	rcs_b = [ rc3, rc4 ]
	rcs_c = [ rc5, rc6 ]
	rcs_d = [ rc7, rc8 ]
	rcs_e = [ rc9, rc10 ]
	rcs_f = [ rc11, rc12 ]
	rcs_g = [ rc13,rc14 ]
	rcs_h = [ rc8, rc14 ] 
	
	
	_path = os.path.join(local,"MD_QMMM",method,"production.ptGeo")
	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.Save_DCD()
	trajAn.DistancePlots(rcs_,False)
	'''
	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.DistancePlots(rcs_b,False)
	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.DistancePlots(rcs_c,False)

	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.DistancePlots(rcs_d,False)
	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.DistancePlots(rcs_g,False)	
	trajAn = TrajectoryAnalysis(_path,proj.system,n)
	trajAn.DistancePlots(rcs_h,False)
	#trajAn.ExtractFrames()
	'''

	
#=======================================================================
def Run_Scan_1D(name="SCAN1D",coord="1"):
	'''
	'''
	
	base_pkl = os.path.join(local,"OPT_QMMM","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,name))
	proj.system.coordinates3 = Unpickle( os.path.join(local,"OPT_QMMM","opt_Steep.ptGeo","frame11.pkl") )[0]
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":"pm6","QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	

	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.438:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	OW_715 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	H1_715  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	atoms1 = [ aco, H1  ] 
	atoms2 = [ OW, H1_715 ]
	atoms3 = [ OW, H1, aco ] 
	atoms4 = [ OW_715, H1_715,OW ]
	atoms5 = [ H1_715,OW_715,C ] 
	
	atoms = None 
	
	if coord == "1": atoms = atoms1
	elif coord == "2": atoms = atoms2
	elif coord == "3": atoms = atoms3
	elif coord == "4": atoms = atoms4
	elif coord == "5": atoms = atoms5
			
	parameters = { "ATOMS_RC1":atoms,
					"dincre_RC1":0.1,
					"nsteps_RC1":20  ,
					"ndim": 1,
					"force_constant_1":2000.0,
					"save_format":".dcd",
					"MC_RC1":		True  ,
					"log_frequency":50    ,
					"simulation_type":"Relaxed_Surface_Scan",
					"NmaxThreads":        1}
					
	proj.Run_Simulation(parameters)
	proj.SaveSystem(_cname="scan_initial")
	rc1 = ReactionCoordinate(atoms,True)
	rc1.GetRCLabel(proj.system)
	
	log_path = os.path.join(local,name,"ScanTraj.log")
	parameters = {"xsize":20,
				  "type":"1D",
				  "log_name":log_path,
				  "crd1_label":rc1.label,
				  "analysis_type":"Energy_Plots"}
				  
	proj.Run_Analysis(parameters)
	
#-----------------------------------------------------------------------
def Run_Scan_2D(name,analysis="False"):
	'''
	'''

	base_pkl = os.path.join(local,"OPT_QMMM","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,name))
	proj.system.coordinates3 = Unpickle( os.path.join(local,"OPT_QMMM","opt_Steep.ptGeo","frame11.pkl") )[0]
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":"pm6","QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	
	
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	OW_715 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	H1_715  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	
	atoms1 = [ OW, H1, aco ] 
	atoms2 = [ H1, OW, C ]
	
	'''
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	H2_529 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.529:H2")
	OW_529 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.529:O")
	OW_572 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:O")
	H2_572 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:H2")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O2")
	
	
		
	atoms1 = [OW_572, H2_572, aco ]
	atoms2 = [OW_529, H2_529, OW_572] 
	'''

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	parameters = { "ATOMS_RC1":atoms1	,
					"ATOMS_RC2":atoms2	,
					"dincre_RC1":0.05,
					"dincre_RC2":0.1,  
					"nsteps_RC1":10  ,
					"nsteps_RC2":22  , 
					"rmsGradient":0.12,
					"ndim": 2 		 ,
					"force_constant_1":4000.0,
					"force_constant_2":3200.0,
					"MC_RC1":		True  ,
					"MC_RC2":		True ,
					"log_frequency":50    ,
					"simulation_type":"Relaxed_Surface_Scan",
					"NmaxThreads":        4}
	
	if analysis == "False":	proj.Run_Simulation(parameters)
	
	log_path = os.path.join(local,name,"ScanTraj.log")
	parameters = {"ysize":22,
				  "xsize":10,
				  "type":"2D",
				  "log_name":log_path,
				  "crd1_label":rc1.label,
				  "crd2_label":rc2.label,
				  "contour_lines":10,
				  "analysis_type":"Energy_Plots"}
				  
	proj.Run_Analysis(parameters)
	
def convert():
	'''
	'''
	base_pkl = os.path.join(local,"OPT_QMMM","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"SCAN2D_a"))
	_path = os.path.join(local,"SCAN2D_a","trajA.ptGeo")
	trajAn = TrajectoryAnalysis(_path,proj.system,16)
	trajAn.Save_DCD()

#-----------------------------------------------------------------------
if __name__ == "__main__":
	if 		sys.argv[1] == "opt" : QCMM_OPT(sys.argv[2])
	elif 	sys.argv[1] == "md_runs"  : MD_runs(sys.argv[2])
	elif 	sys.argv[1] == "analysis" : Analysis(sys.argv[2],int(sys.argv[3]))
	elif 	sys.argv[1] == "scan1d"     : Run_Scan_1D(sys.argv[2],sys.argv[3])
	elif 	sys.argv[1] == "scan2d"     : Run_Scan_2D(sys.argv[2],sys.argv[3])
	elif 	sys.argv[1] == "scan_analysis"     : Scan_Analysis()
	else: convert()
	

