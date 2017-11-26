PRO doc_widget1_event, ev

  IF ev.SELECT THEN WIDGET_CONTROL, ev.TOP, /DESTROY

END

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
set_plot,'X'
DeVICE,decomposed=1
;DEVICE,GET_PIXEL_DEPTH=b
;DEVICE,SET_RESOLUTION=[1920,1280]
;DEVICE,SET_PIXEL_DEPTH=24
;LOADCT,39
;LOADCT,13
	!P.background='FFFFFF'xL
	!P.Color='000000'xL
;	TVLCT, 0,0,0, !P.COLOR
;	TVLCT, 255,255,255, !P.Background
;	TVLCT, 0,0,0, !P.Background
;	TVLCT, 255,255,255, !P.COLOR
window,xsize=7680,ysize=4320

window,0,xsize=640,ysize=480
print,!D.window
;print,b
GraphProperties={graph_properties_struct}

GraphProperties.CHARSIZE=2

GraphProperties.CONTOUR_AREA=[0.1,0.1,0.95,0.8]
GraphProperties.LEGEND_AREA=[.15,.85,.15,.95]
	GraphProperties.VIEW_AREA=[min(glat),min(zalt),max(glat),max(zalt)]
	GraphProperties.VIEW_AREA=[0,0,nz,nf]
	;GraphProperties.VIEW_AREA=[15,0,25,3500]
grphProp=GraphProperties
print,grphProp
;contour_plot,dene(*,*,10),glat,zalt,grphProp
contour_plot,indgen(nz,nf),bindgen(nz) # replicate(1b,nf),transpose(bindgen(nf) # replicate(1b,nz)),grphProp

	imagePngTrueParam=1
	imageToSave=tvrd(TRUE=imagePngTrueParam)
	filename=filepath('tmp.png' ,$
	ROOT_DIR='/home/shulindinae',$
	SUBDIRECTORY=['Pictures'])
	TVLCT, r, g, b, /Get
	write_png,filename,imageToSave,r,g,b;,/ORDER




;a=tvrd(/true)
;set_plot,'X'
;	!P.background=255
;	!P.Color=0
;	TVLCT, 0,0,0, !P.COLOR
;	TVLCT, 255,255,255, !P.Background


;window,0,xsize=640,ysize=480
;TV,a
;xloadct
if(common_app_gui) then begin	
	base = WIDGET_BASE(/COLUMN)

  button = WIDGET_BUTTON(base, value='Done')

  WIDGET_CONTROL, base, /REALIZE

  XMANAGER, 'doc_widget1', base
endif

end