pro sami2,VERBOSE=_Verbose,DATAPATH=_DataPath,GUI=_Gui
!PATH="./common/:"+!PATH
@common/common_app_inc.pro
@common/common_in_data_inc.pro


	if _DataPath EQ !NULL THEN _DataPath=""
	
	common_app_Verbose = KEYWORD_SET(_Verbose)
	common_app_DataPath= _DataPath
	common_app_gui=KEYWORD_SET(_Gui)
	Trace,"Verbose mode is on."
	Trace,"Path to data: "+common_app_DataPath
	Trace,"GUI is "+(common_app_gui?"ON":"OFF")
    
	common_in_data_LoadData,common_app_DataPath

	
	


end