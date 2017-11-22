pro common_in_data_LoadData,DataPath
COMMON common_in_data
	time=read_time_dat(DataPath)
	nt=n_ELEMENTS(time)
	glat=READ_GLATU_DAT(DataPath,nz,nf)
	zalt=READ_ZALTU_DAT(DataPath,nz,nf)
	deni=READ_DENIU_DAT(DataPath,nz,nf,nion,nt)
	dene=total(deni,3)




end


pro common_in_data
COMMON common_in_data
COMMON block_constants

    ions_names=['H+','O+','NO+','O2+','He+','N2','N']
    ions_out=[1,1,0,0,1,0,0]
    
    	nf=52
	nz=51
	nion=7
    


end