FUNCTION READ_DENIU_DAT,DataDir,NZ,NF,NION,NT
filename=DataDir+'deniu.dat'

deni = fltarr(nz,nf,nion,nt)
OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/F77_UNFORMATTED

  denitmp = fltarr(nz,nf,nion)
  for i = 0,nt-1 do begin
    readu,LUN,denitmp
    deni(*,*,*,i) = denitmp
  endfor
FREE_LUN,LUN


RETURN,deni
END

