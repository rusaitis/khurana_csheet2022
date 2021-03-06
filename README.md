# Khurana Jupiter Current Sheet Structure Model 2022
[![DOI](https://zenodo.org/badge/493023766.svg)](https://zenodo.org/badge/latestdoi/493023766)


## Description

The program calculates the current sheet structure at Jupiter at a particular point in space and time.

## Instructions

Call the Fortran subroutine 'csheet_struc' with the input parameters of the desired location in space, given seperately as 'R', 'theta', 'phi' in System III coordinates, and 'XJSO' and 'YJSO' in Jupiter-Sun-Orbit coordinates, as well as time argument 'ctime' in UNIX time. The output is 'ZNS3', the height of the current sheet in System III.

## Main Fortran Subroutine

``` Fortran
SUBROUTINE csheet_struc(ZNS3,R,theta,phi,XJSO,YJSO,ctime)

Arguments:
	INPUT:  R,theta,phi : position in System III
		XJSO,YJSO   : position in JSO
		ctime       : UNIX time
	OUTPUT: ZNS3        : hight of the current sheet in System III
```

## Other Useful Information
Helpful information on the Jupiter coordinate systems, such as System III and JSO:
https://lasp.colorado.edu/home/mop/files/2015/02/CoOrd_systems7.pdf  
F. Bagenal & R. J. Wilson, LASP – University of Colorado

## Author
Krishan Khurana, UCLA.  
05/16/2022

## Contacts
For any questions, please email Krishan Khurana.

Krishan Khurana  
kkhurana@igpp.ucla.edu  
Institute of Geophysics and Planetary Physics  
University of California at Los Angeles  
Los Angeles, CA, 90095  
