pro common_in_data_LoadData,DataPath
COMMON common_in_data
	nf=52
	nz=51
	nion=7
    
	time=read_time_dat(DataPath)
	nt=n_ELEMENTS(time)
	glat=READ_GLATU_DAT(DataPath,nz,nf)
	zalt=READ_ZALTU_DAT(DataPath,nz,nf)
	deni=READ_DENIU_DAT(DataPath,nz,nf,nion,nt)
	dene=total(deni,3)




end


pro common_in_data
end