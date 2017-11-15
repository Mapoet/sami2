FUNCTION READ_TIME_DAT_1,DataDir
filename=DataDir+'time.dat'

timerow=Lindgen(4)
countofrows=0
OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/F77_UNFORMATTED
while not eof(LUN) do begin
READU,LUN,timerow
countofrows=countofrows+1
endwhile
FREE_LUN,LUN

time=Lindgen(4,countofrows)
countofrows=0

OPENR,LUN,filename,/GET_LUN,/swap_IF_LITTLE_ENDIAN,/F77_UNFORMATTED
while not eof(LUN) do begin
READu,LUN,timerow
time(*,countofrows)=timerow
countofrows=countofrows+1
endwhile
FREE_LUN,LUN

RETURN, time
END

function TimeTitleUT,timep
	hr  = timep.HH
  	phr  = string(format='(i2)',hr)
  	while(((i=strpos(phr,' '))) ne -1) do strput,phr,'0',i
  	min  = timep.mm
	sec  = timep.ss
  	pmin = string(format='(i2)',min)
 return, phr+':'+pmin+' UT'
end


FUNCTION READ_TIME_DAT,DataDir
time=READ_TIME_DAT_1(DataDir)
void={timepoint_struct,$
	NTM:0,$
    HH:0.0,$
    mm:0.0,$
    ss:0.0,$
	TITLEUT:"",$
	TITLELT:""$
}

len=(SIZE(time,/DIMENSIoNS))[1]
timep={timepoint_struct}
timep.NTM=time[0,0]
timep.HH=time[1,0]
timep.mm=time[2,0]
timep.ss=time[3,0]
timep.TITLEUT=TimeTitleUT(timep)

ntms=[timep]


for i = 1,Len-1 do begin
timep.NTM=time[0,i]
timep.HH=time[1,i]
timep.mm=time[2,i]
timep.ss=time[3,i]
timep.TITLEUT=TimeTitleUT(timep)

ntms=[ntms,timep]
end


RETURN, ntms

END