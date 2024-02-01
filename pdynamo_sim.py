import sys                                                                                
sys.path.append("/home/igorchem/pDynamo3_scripts")
sys.path.append("/home/igorchem/programs/pDynamo3_scripts")

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
	parameters_a = {"simulation_type":"Geometry_Optimization"	,
					"maxIterations":10000						,
					"rmsGradient":1.0 							,
					"save_frequency": 20						,
					"save_format":".dcd"						, 
					"trajectory_name":"opt_first.ptGeo"				}
					
	proj.Run_Simulation(parameters_a)
	
	xsi = len( glob.glob( os.path.join( local,"MM_SetUp","opt_first.ptGeo","*.pkl") ) )
	_path = os.path.join(local,"MM_SetUp","opt_first.ptGeo")
	trajAn = TrajectoryAnalysis(_path,proj.system,xsi)
	trajAn.CalculateRG_RMSD()
	trajAn.Save_DCD()
	trajAn.PlotRG_RMS()
	
	_pattern = "*:CO2.414:C"
	proj.Spherical_Pruning(_pattern,30.0)
	proj.Setting_Free_Atoms(_pattern,20.0)	
	proj.SaveSystem()
	
#-----------------------------------------------------------------------
def Def_MM_Sys2():
	

	proj = SimulationProject.From_PKL(os.path.join( local,"MM_SetUp","sys01amber.pkl"), os.path.join( local,"MM_SetUp2") )
	
	parameters_a = {"simulation_type":"Geometry_Optimization"	,
					"maxIterations":10000						,
					"rmsGradient":0.1 							,
					"save_frequency": 20						,
					"Debug":True								,
					"save_format":".dcd"						, 
					"trajectory_name":"opt_second.ptGeo"		}
					
	proj.Run_Simulation(parameters_a)
		
	xsi = len( glob.glob( os.path.join( local,"MM_SetUp2","opt_second.ptGeo","*.pkl") ) )
	_path = os.path.join(local,"MM_SetUp2","opt_second.ptGeo")
	trajAn = TrajectoryAnalysis(_path,proj.system,xsi)
	trajAn.CalculateRG_RMSD()
	trajAn.Save_DCD()
	trajAn.PlotRG_RMS()	
	proj.SaveSystem()


#================================
# 414+438+79+40+61+628+529+447+715+572+397+330+256
def QCMM_OPT(method="am1"):
	'''
	'''
	
	proj = SimulationProject.From_PKL(os.path.join(local,"MM_SetUp2","sys01amber.pkl"),os.path.join(local,"OPT_QMMM",method) )
	
	co2_a  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:*")
	co2_b  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.438:*")
	#cat1 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.79:*")
	#cat2 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.40:*")
	cat3 = AtomSelection.FromAtomPattern(proj.system,"*:c4c.61:*")
	wat1 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:*")
	wat2 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.529:*")
	#wat3 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.447:*")
	#wat4 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.715:*")
	wat5 = AtomSelection.FromAtomPattern(proj.system,"*:WAT.572:*")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:*")
	#aco2 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.330:*")
	#aco3 = AtomSelection.FromAtomPattern(proj.system,"*:ACO.256:*")
	
	
	selections = [ co2_a,co2_b,cat3,wat1,wat2,wat5,aco ]
	_parameters_b = {"method_class":"SMO","Hamiltonian":method,"QCcharge":0,"multiplicity":1,"region":selections}
	proj.DEBUG = True
	proj.Set_QC_Method(_parameters_b)
	proj.Energy
		
	parameters = {"maxIterations":2500 						,
				"rmsGradient":0.1   						,
				"log_frequency":10 							,
				"simulation_type":"Geometry_Optimization"	,
				"optmizer":"Stepeest_Descent"				,
				"save_frequency" :1 						,
				"save_format":".dcd"						,
				"trajectory_name":"opt.ptGeo"				,
				"Debug":True								,
				"save_pdb": True 							}

	proj.Run_Simulation(parameters)	
	proj.SaveSystem()
	
	xsi = len( glob.glob( os.path.join( local,"OPT_QMMM",method,"opt.ptGeo","*.pkl") ) )
	if xsi > 0:
		_path = os.path.join(local,"OPT_QMMM",method,"opt.ptGeo")	
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
def Run_Scan_1D(name="SCAN1D",coord="1",method="am1"):
	'''
	'''
	
	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,name))
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":method,"QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	proj.Energy
	

	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	

	atoms1 = [ OW[0], H1[0], aco[0]  ]
	atoms2 = [ OW[0], C[0] ]
	

	atoms 	= None 
	dincre 	= 0.05
	fc 		= 3000.0
	steps   = 12 
	if coord == "1": 
		atoms = atoms1
	elif coord == "2": 
		dincre = -0.1
		atoms = atoms2
		steps = 24
		fc    = 2000.0
				
	parameters = { "ATOMS_RC1":atoms							,
					"dincre_RC1":dincre							, 
					"nsteps_RC1":steps							,
					"ndim": 1									,
					"force_constant_1":fc						,
					"save_format":".dcd"						,
					"MC_RC1":		True  						,
					"optmizer":"SteepestDescent"				,
					"log_frequency":50    						,
					"simulation_type":"Relaxed_Surface_Scan"	,
					"NmaxThreads":        1						}
					
	proj.Run_Simulation(parameters)
	proj.SaveSystem(_cname="scan_initial")
	
	rc1 = ReactionCoordinate(atoms,True)
	rc1.GetRCLabel(proj.system)
	
	log_path = os.path.join(local,name,"ScanTraj.log")
	parameters = {"xsize":steps,
				  "type":"1D",
				  "log_name":log_path,
				  "crd1_label":rc1.label,
				  "analysis_type":"Energy_Plots"}
				  
	proj.Run_Analysis(parameters)
	
#-----------------------------------------------------------------------
def Run_Scan_2D(name,analysis="False"):
	'''
	'''

	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,name))
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":name,"QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	proj.Energy
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	
	atoms1 = [ OW[0], H1[0], aco[0] ] 
	atoms2 = [ OW[0], C[0] ]

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	stepsx= 12
	stepsy= 24
	
	parameters = { "ATOMS_RC1":atoms1	,
					"ATOMS_RC2":atoms2	,
					"dincre_RC1":0.05,
					"dincre_RC2":-0.1,  
					"nsteps_RC1":stepsx,
					"nsteps_RC2":stepsy, 
					"rmsGradient":0.1,
					"optmizer":"SteepestDescent",
					"ndim": 2 		 ,
					"force_constant_1":3000.0,
					"force_constant_2":2000.0,
					"MC_RC1":		True  ,
					"MC_RC2":		True ,
					"log_frequency":50    ,
					"simulation_type":"Relaxed_Surface_Scan",
					"NmaxThreads":        8}
	
	if analysis == "False":	proj.Run_Simulation(parameters)
	
	log_path = os.path.join(local,name,"ScanTraj.log")
	parameters = {"xsize":stepsx,
				  "ysize":stepsy,
				  "type":"2D",
				  "xlim_list":[-0.6,0.6],
				  "ylim_list":[3.5,1.30],
				  "log_name":log_path,
				  "crd1_label":rc1.label,
				  "crd2_label":rc2.label,
				  "contour_lines":16,
				  "analysis_type":"Energy_Plots"}
				  
	proj.Run_Analysis(parameters)

#-----------------------------------------------------------------------
def convert():
	'''
	'''
	base_pkl = os.path.join(local,"OPT_QMMM","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"SCAN2D_a"))
	_path = os.path.join(local,"SCAN2D_a","trajA.ptGeo")
	trajAn = TrajectoryAnalysis(_path,proj.system,16)
	trajAn.Save_DCD()
	
	
#-----------------------------------------------------------------------
def Refine_MOPAC(run="True",folder="Refine",cut=0.0):
	'''
	'''
	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,folder))
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	atoms1 = [ OW[0], H1[0], aco[0] ] 
	atoms2 = [ OW[0], C[0] ]

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	cqr = False
	if cut > 0.001: cqr = True
	print(run)
		
	methods = ["am1","pm3","pm6","rm1","pddgpm3","am1dphot"]
	_path = os.path.join( os.path.join(local,"SCANS2D","pm3","ScanTraj.ptGeo") )
	parameters = { "xnbins":12			,
				   "ynbins":24			,
				   "mopac_keywords":["ITRY=5000"] ,
				   "source_folder":_path,
				   "folder":os.path.join(local, folder),
				   "charge":0		    ,
				   "multiplicity":1 	,
				   "change_qc_region":cqr                   ,
				   "center": [26.732,7.702,29.268]           ,
				   "radius": cut ,
				   "methods_lists":methods,	
				   "NmaxThreads":4		,
				   "simulation_type":"Energy_Refinement",
				   "Software":"pDynamo"	}
	#---------------------------------------------
	if run == "True": proj.Run_Simulation(parameters)	
	parameters= {"xsize":12,
				 "ysize":24,
				 "ndim":2,
				 "contour_lines":14,
				 "xlim_list":[-0.6,0.6],
				 "ylim_list":[3.5,1.30],
				 "log_name":os.path.join(local,folder,"energy.log"),
				 "crd2_label":rc1.label,
				 "crd1_label":rc1.label,"multiple_plot":"log_names",
				 "analysis_type":"Energy_Plots","type":"2DRef" }
	#--------------------------------------------
	proj.Run_Analysis(parameters)
#-----------------------------------------------------------------------
def Refine_ORCA():
	'''
	'''	
	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"ORCA_REF"))
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	atoms1 = [ OW[0], H1[0], aco[0] ] 
	atoms2 = [ OW[0], C[0] ]

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	cqr = False
	if cut > 0.001: cqr = True
	print(run)
	_path = os.path.join( os.path.join(local,"SCANS2D","pm3","ScanTraj.ptGeo") )
	#---------------------------------------------
	parameters = { "xnbins":12			                                           ,
				   "ynbins":24			                                           ,
				   "source_folder":_path                                           ,
				   "orca_method":"b3lyp"                                           ,
				   "basis":"6-31G*"                                                 , 
				   "folder":os.path.join(local,"ORCA_REF2")	        		,
				   "charge":0		                                               ,
				   "multiplicity":1 	                                           ,
				   "restart":False                                                 ,                                             
				   "NmaxThreads":4	                                   ,
				   "simulation_type":"Energy_Refinement"                           ,
				   "Software":"ORCA"	                                           }
				   
	#---------------------------------------------
	proj.Run_Simulation(parameters)

	parameters= {"xsize":12,"ysize":24,
				 "log_name":os.path.join(local,"ORCA_REF2","energy.log"),
				 "crd1_label":rc1.label,
				 "crd2_label":rc2.label,
				 "xlim_list":[-0.6,0.6],
				 "ylim_list":[3.5,1.30],
				 "contour_lines":16,
				 "analysis_type":"Energy_Plots","type":"2DRef" }
	proj.Run_Analysis(parameters)
#-----------------------------------------------------------------------
def Traj_1D(folder,cut=0.0,program="pDynamo",np=4,run="True"):
	'''
	'''
	
	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"Traj1D"))
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	atoms1 = [ OW[0], H1[0], aco[0] ] 
	atoms2 = [ OW[0], C[0] ]

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	_path = os.path.join( os.path.join(local,"Traj1D.ptGeo") )
	methods = ["am1","pm3","pm6","rm1","pddgpm3"]
	if program == "mopac": methods = ["am1","pm3","pm6","rm1","pm7"]
	parameters = { "xnbins":28			,
				   "mopac_keywords":["ITRY=5000"] ,
				   "source_folder":_path,
				   "folder":os.path.join(local, folder),
				   "charge":0		    ,
				   "multiplicity":1 	,   
				   "NmaxThreads":np		,
				   "simulation_type":"Energy_Refinement",
				   "Software":program	}
	if not program == "ORCA":
		parameters["methods_lists"] = methods
		if cut > 0.0:
			parameters["change_qc_region"] = True
			parameters["radius"] = cut
			parameters["center"] = [26.732,7.702,29.268] 
	elif program == "ORCA":
		parameters["orca_method"]="b3lyp"
		parameters["basis"]="6-311+G*"   
		
		
	#---------------------------------------------
	if run == "True": proj.Run_Simulation(parameters)	
	parameters= {"xsize":28,
				 "ndim":1,
				 "log_name":os.path.join(local,folder,"energy.log"),
				 "crd1_label":"frames","multiple_plot":"log_names",
				 "analysis_type":"Energy_Plots","type":"1DRef" }				
	#--------------------------------------------
	proj.Run_Analysis(parameters)
	trajAn = TrajectoryAnalysis(_path,proj.system,28)
	trajAn.Save_DCD()
	
#-----------------------------------------------------------------------
def Free_energy(NmaxThreads=8,run="True"):
	'''
	'''
	base_pkl = os.path.join(local,"OPT_QMMM","am1","sys01amber.pkl")
	proj = SimulationProject.From_PKL(base_pkl,os.path.join(local,"freeenergy"))
	
	selections = proj.system.qcState.pureQCAtoms
	_parameters_b = {"method_class":"SMO","Hamiltonian":"pm3","QCcharge":0,"multiplicity":1,"region":selections}
	proj.Set_QC_Method(_parameters_b)
	proj.Energy
	
	C  = AtomSelection.FromAtomPattern(proj.system,"*:CO2.414:C")
	OW = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:O")
	H1  = AtomSelection.FromAtomPattern(proj.system,"*:WAT.628:H1")
	aco = AtomSelection.FromAtomPattern(proj.system,"*:ACO.397:O1")
	
	atoms1 = [ OW[0], H1[0], aco[0] ] 
	atoms2 = [ OW[0], C[0] ]

	rc1 = ReactionCoordinate(atoms1,True)
	rc1.GetRCLabel(proj.system)
	rc2 = ReactionCoordinate(atoms2,True)
	rc2.GetRCLabel(proj.system)
	
	_path = os.path.join( os.path.join(local,"Traj1D.ptGeo") )
	
	USparameters = { "ATOMS_RC1":atoms1				,
				   "ATOMS_RC2":atoms2				,
				   "ndim":2 						,
				   "sampling_factor":500			,
				   "equilibration_nsteps":10000 	,
				   "production_nsteps":50000		,
				   "source_folder":_path 			,
				   "pressure":20.0					,
				   "pressure_coupling":True			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "MC_RC2":True					,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":NmaxThreads		}

	if run == "True": proj.Run_Simulation(USparameters)

	_path = os.path.join( scratch_path, "freeenergy")
	PMFparameters = { "source_folder":_path,
				   "xnbins":10           ,
				   "ynbins":10           ,
				   "ywindows":0          ,
				   "xwindows":28         ,
				   "crd1_label":rc1.label,
				   "crd2_label":rc2.label,
				   "oneDimPlot":True     ,
				   "analysis_type":"PMF_Analysis",
				   "temperature":300.15	 }
	proj.Run_Analysis(PMFparameters)
#-----------------------------------------------------------------------
if __name__ == "__main__":
	if 		sys.argv[1] == "MM": Def_MM_Sys()
	if 		sys.argv[1] == "MM2": Def_MM_Sys2()
	elif 	sys.argv[1] == "opt" : QCMM_OPT(sys.argv[2])
	elif 	sys.argv[1] == "md_runs"  : MD_runs(sys.argv[2])
	elif 	sys.argv[1] == "analysis" : Analysis(sys.argv[2],int(sys.argv[3]))
	elif 	sys.argv[1] == "scan1d"     : Run_Scan_1D(name=sys.argv[2],coord=sys.argv[3],method=sys.argv[4])
	elif 	sys.argv[1] == "scan2d"     : Run_Scan_2D(sys.argv[2],sys.argv[3])
	elif 	sys.argv[1] == "scan_analysis"     : Scan_Analysis()
	elif 	sys.argv[1] == "orca"     : Refine_ORCA()
	elif 	sys.argv[1] == "mopac"     : Refine_MOPAC(sys.argv[2],sys.argv[3],float(sys.argv[4]))
	elif 	sys.argv[1] == "traj1d": Traj_1D(sys.argv[2],float(sys.argv[3]),sys.argv[4],int(sys.argv[5]),sys.argv[6])
	elif 	sys.argv[1] == "FE":     Free_energy( int(sys.argv[2]),sys.argv[3] )
		
		
	else: convert()
	

