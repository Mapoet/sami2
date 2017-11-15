FUNCTION READ_GLONU_DAT,DataDir,NZ,NF
filename=DataDir+'glonu.dat'
glon = fltarr(nz,nf)

OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/f77_unformatted
readu,LUN,glon
FREE_LUN,LUN

RETURN,GLON
END

