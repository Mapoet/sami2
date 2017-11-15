FUNCTION READ_ZALTU_DAT,DataDir,NZ,NF
filename=DataDir+'zaltu.dat'
zalt = fltarr(nz,nf)

OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/f77_unformatted
readu,LUN,zalt
FREE_LUN,LUN

RETURN,zalt
END

