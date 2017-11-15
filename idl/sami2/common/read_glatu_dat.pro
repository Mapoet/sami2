FUNCTION READ_GLATU_DAT,DataDir,NZ,NF
filename=DataDir+'glatu.dat'
glat = fltarr(nz,nf)
OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/f77_unformatted
readu,LUN,glat
FREE_LUN,LUN
RETURN,GLat
END

