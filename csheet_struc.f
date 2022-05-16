C******************************************************************************
	SUBROUTINE csheet_struc(ZNS3,R,theta,phi,XJSO,YJSO,ctime)
C     !Version 2004 08 15
C	INPUT:  R,theta,phi  ->in System III
C             XJSO,YJSO    ->in JSO
C             ctime        ->UNIX cline time
C     OUTPUT: ZNS3         ->hight of the current sheet in System III
C
      IMPLICIT REAL*8(A-H,O-Z)


	DATA XH /47/ 
	DATA dLamVIP4   /5.92068/
	DATA dThetaVIP4 /0.16615534478986/

	RHO=R*dsin(theta)

	CALL JSun(ctime,stheta,sphi,phase)

	RLT=datan2(YJSO,XJSO) !zero at noon


	dlamprime=waveDelay(RLT,RHO)+bendBack(RLT,RHO)+dLamVIP4
	rho1=dsqrt((XH*dtanh(dabs(XJSO/XH)))**2+YJSO**2)
	term2=rho1*dtan(dThetaVIP4)
	ZNS3=term2*dcos(phi-dlamprime)+XJSO*(1.0-dtanh(dabs(XH/XJSO)))
     ~	*dtan(stheta)


	RETURN
	END
C******************************************************************************
C******************************************************************************
	FUNCTION waveDelay(RLT,RHO)
      IMPLICIT REAL*8(A-H,O-Z)


	Data B0 / 3.79406064748764D-002/
      Data C1 / 1.34045565128326 / !C1=Domegaj*rho2/V0
								 !with Domegaj=36.26125*pi/180
	Data rho2/ 83.3625411987305/	 
	Data B1/ 0.748919927647621/,	 Phi1/ -1.04848085782547/
	Data B2/ 0.614373903493367/,	 Phi2/ -1.95757783830202/
	Data B3/ 0.541435668088830/,	 Phi3/ -1.25514203423947/
	Data B4/ 7.02600475148873D-002/, Phi4/ -1.25542899066018/
	Data B5/ 26.1889910697937/

	rhoor2=RHO/rho2
	rhooB5=RHO/B5

	waveDelay=B0-C1*dlog(dcosh(rhoor2))*
     ~(1+dtanh(rhooB5)*(B1*dcos(RLT-Phi1)+B2*dcos(2*RLT-Phi2)+
     ~B3*dcos(3*RLT-Phi3)+B4*dcos(4*RLT-Phi4)))


	RETURN
	END
C******************************************************************************
C******************************************************************************
	FUNCTION bendBack(RLT,RHO)
      IMPLICIT REAL*8(A-H,O-Z)


	Data rho0   /33.0/
	Data A1  / 1.062528141129624e-004 /!divided by rho1 (1.6*rho0=52.8)
	Data A2  / 5.950493765539412e-005 /!divided by rho1 (1.6*rho0=52.8)
	Data A3  / 6.788463593277820e-005 /!divided by rho1 (1.6*rho0=52.8)
	Data A4  / 3.009364679056859e-005 /!divided by rho1 (1.6*rho0=52.8)
	Data Phi1 / 1.79696364159234/
	Data Phi2 / 0.41408648226213/
	Data Phi3 /-0.22365844912454/
	Data Phi4 /-0.42846193975833/


      Data Cb1n   /-0.18591516533271      /!Cb1=Cb1*rho0
	Data Cb1logc/-8.534743746501581e-005/!Cb1n*dlog(dcosh(dabs(1/rho0)))
	Data ddK2   /-0.01893939393939      /!dK2=-1/(1.6*rho0)

	
	ddK1=A1*dcos(  RLT-Phi1)
	ddK3=A2*dcos(2*RLT-Phi2)
	ddK5=A3*dcos(3*RLT-Phi3)
	ddK7=A4*dcos(4*RLT-Phi4)

	T3=dexp(ddK2*rho)/(ddK2**2)*(ddK2*rho-1) -
     ~   dexp(ddK2    )/(ddK2**2)*(ddK2    -1)
	T4=ddK1/2+ddK3/2+ddK5/2+ddK7/2
	rhoSq=rho**2

	bendBack=Cb1n*dlog(dcosh(dabs(rho/rho0)))
     ~	-Cb1logc+T4*rhoSq+T3*(-ddK1-ddK3-ddK5-ddK7)-T4


	RETURN
	END
C******************************************************************************
c******************************************************************************
	Subroutine JSun(ctime,stheta,sphi,phase)
c	INPUT:  ctime of the data point 
C	OUTPUTS: stheta, sphi, latitude and longitude  (in radians) of the Sun in system III (RH).  
C     OUTPUTS: phase, Orbital phase angle of Jupiter (in radians) 
c	The equations are written in etime the J200 time convention followed by PDS. 
c	We first convert ctime to etime

c	Last updated August 12, 2002.
c
C	The program first calculates the direction of the sun in non-rotating  
C	coordinates.from equations of the type:
c	theta=a1*cos(omegay*t)+a2*sin(omegay*t)+a3*cos(2.*omegay*t)+
c	a4*sin(2.*omegay*t)+a5
c	fphi=b1*cos(omegay*t)+b2*sin(omegay*t)+b3*cos(2.*omegay*t)+
c	b4*sin(2.*omegay*t)+b5
c	Then we rotate into the System III coordinates

c	ctime is in double precision. theta and phi are single precision 
c	+variables.
c

c	
c	omega is jupiter's rotation rate
c	omegay is jupiter's yearly orbital rate

c	Initialize variables
      IMPLICIT REAL*8(A-H,O-Z)
	Dimension aa(7),bb(7), x(7)
	Real*8 ctime,etime,etime1,fphi,yrjup,omega,D360,omegay,t,year
	Real*8 stheta,sphi,phase
	Parameter (PI=3.1415927,twopi=2.*PI,radian=PI/180.,degree=180./PI)
	Parameter (yrjup=.1185652502D2*86400d0*.36525d3)
	Parameter (omega=870.536/86400.0,omegay=2.*PI/yrjup)
	Parameter (year=86400D0*.36525D3,D360=360.0)
	Parameter (etime1=-8.25767955817479d8)
	Parameter (three=3.123*radian)
	Parameter (tan3=0.054560676,sin3=0.054479647,cos3=0.99851488)
	Data aa /0.14347029,3.1145815,-0.12025561,0.093909436,
     ~-0.39321884e-5,0.10194945e-3,-0.12799464/
	Data bb /-4.5467523,3.1848875,-0.16329986,-0.09776818,
     ~0.17556527e-3,-0.01978317,44.55915/
	

c	First calculate the latitude and longitude in non-rotating Jupiter coordinates.

c	Calculate the best fit theta and fphi
	t=etime(ctime)-etime1
	x(1)=dcos(omegay*t)
	x(2)=dsin(omegay*t)
	x(3)=dcos(2.*omegay*t)
	x(4)=dsin(2.*omegay*t)
	x(5)=(t/year)**2
	x(6)=t/year
	x(7)=1.0
	stheta=0.
c	fphi is phi in Jupiter fixed (non-rotating) coordinate
	fphi=0.0
	Do 1 j=1,7
	fphi=fphi+bb(j)*x(j)
    1 stheta=stheta+aa(j)*x(j)
C	Now rotate the longitude to Jupiter System III
c	First Add the rotation of Jupiter around the Sun.
c	fphi is the phi of the Sun as unspinning Jupiter goes around the sun
	fphi=DMod(fphi+t/yrjup*360.d0, D360)
c	Next add the rotation of Jupiter around its axis.
	sphi=DMod(fphi-t*omega, D360)
	if (sphi .lt. 0.0) sphi = sphi+360.0
	sphi=sphi*radian
	stheta=stheta*radian
c	Now compute the orbital phase (called phi2 or phase here)
c	There are two solutions to the problem. Only one that is close to phi2b is correct.
	if(stheta .ge. three)stheta=(three)
	if(-stheta .ge. three)stheta=-(three)
      phi21=(sphi+PI)+dacos(dtan(stheta)/tan3)
	phi22=(sphi+PI)-dacos(dtan(stheta)/tan3)
	phi21=dMOD(phi21,twopi)
	phi22=dMOD(phi22,twopi)
	if(phi22 .lt. 0) phi22=phi22+twopi
	dphi=fphi+48.23012
	phi2b=sphi-dphi*radian
	phi2b=dMOD(phi2b,twopi)
	if(phi2b .lt. 0) phi2b=phi2b+twopi
	phase=phi21
	dif2=dabs(phi22-phi2b)
	if(dif2 .gt. 350*radian) dif2=twopi-dif2
	dif1=dabs(phi21-phi2b)
	if(dif1 .gt. 350*radian) dif1=twopi-dif1
	if(dif2 .lt. dif1) phase=phi22
c	phase=phi2b
	phase=dMOD(phase,twopi)
c	Write(*,*)stheta*degree,phase*degree,phi2b*degree
	Return
      end 

c*******************************************************************
c*******************************************************************
	Function etime(ctime)
	Real*8 ctime, etime
	Tcor=0.
	if(ctime .ge. 189302400.000) Tcor = Tcor+10	
	if(ctime .ge. 205027200.000) Tcor = Tcor+1	
	if(ctime .ge. 220924800.000) Tcor = Tcor+1
	if(ctime .ge. 252460800.000) Tcor = Tcor+1
	if(ctime .ge. 283996800.000) Tcor = Tcor+1
	if(ctime .ge. 315532800.000) Tcor = Tcor+1
	if(ctime .ge. 347155200.000) Tcor = Tcor+1
	if(ctime .ge. 378691200.000) Tcor = Tcor+1
	if(ctime .ge. 410227200.000) Tcor = Tcor+1
	if(ctime .ge. 441763200.000) Tcor = Tcor+1
	if(ctime .ge. 489024000.000) Tcor = Tcor+1
	if(ctime .ge. 520560000.000) Tcor = Tcor+1
	if(ctime .ge. 552096000.000) Tcor = Tcor+1
	if(ctime .ge. 615254400.000) Tcor = Tcor+1
	if(ctime .ge. 694224000.000) Tcor = Tcor+1
	if(ctime .ge. 757382400.000) Tcor = Tcor+1
	if(ctime .ge. 788918400.000) Tcor = Tcor+1
	if(ctime .ge. 836179200.000) Tcor = Tcor+1
	if(ctime .ge. 867715200.000) Tcor = Tcor+1
	if(ctime .ge. 899251200.000) Tcor = Tcor+1
	if(ctime .ge. 946684800.000) Tcor = Tcor+1
	if(ctime .ge. 993945600.000) Tcor = Tcor+1
	if(ctime .ge. 1041379200.000) Tcor = Tcor+1
	etime=ctime+dble(Tcor)-.1072958367816D10
	Return
	End


C********************************************************************

	Function ctimer(iyr,imon,iday,ihr,imin,sec)
	integer iyr,imon,iday,ihr,imin
	integer*4 ndays,doy
	Real*8 sec,ctimer

c	First calculate the number of days from Jan 1, 1966
		ndays=0
		If	(iyr .ge. 1966) then
		Do 1 i=1966,iyr-1
		ndays=ndays+365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or. 
	+	(mod(i,400)) .eq. 0) ndays=ndays+1
   1		continue
c	Now add the number of days of the current year
    	ndays=ndays+doy(iyr,imon,iday)-1	
    	ctimer = dble(ndays)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
	+dble(sec)
	go to 4
	end if

c	Calculate the seconds for the negative years
	If	(iyr .lt. 1966) then
	Do 2 i=iyr, 1965
	ndays=ndays-365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or. 
	+	(mod(i,400)) .eq. 0) ndays=ndays-1
    2 continue
c	Now subtract the number of days of the current year	
	ndays=ndays+doy(iyr,imon,iday)
	end if
    	ctimer = dble(ndays-1)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
	+dble(sec)
    4 continue
	idyr=doy(iyr,imon,iday)
	Write(*,*) " doy is ",idyr
	Write(*,*)" ndays is ",ndays
	Return 
	End
c*******************************************************************
	Function doy(iyr,imon,iday)
	Integer mon(12),doy
	data mon /31,28,31,30,31,30,31,31,30,31,30,31/
	doy=0
	Do 1 i = 2,imon
    	doy=doy+mon(i-1)
c	Add an extra day for February
	if(i .eq. 3) then
		if(((mod(iyr,4) .eq. 0) .and. (mod(iyr,100) .ne. 0)) .or. 
	+	(mod(iyr,400)) .eq. 0) doy=doy+1
	end if

    1 continue
	doy=doy+iday
	Return
	End
c*******************************************************************
