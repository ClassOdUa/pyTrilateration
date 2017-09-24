#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy

def f_tralateration(p_LatA, p_LngA, p_DistA, p_LatB, p_LngB, p_DistB, p_LatC, p_LngC, p_DistC):

	#assuming elevation = 0
	l_earthR = 6372795

	#using authalic sphere
	#if using an ellipsoid this step is slightly different
	
	#TODO: find and change the computation for ellipsoid
	
	#Convert geodetic Lat/Long to ECEF xyz
	#   1. Convert Lat/Long to radians
	#   2. Convert Lat/Long(radians) to ECEF
	l_xA = l_earthR *(math.cos(math.radians(p_LatA)) * math.cos(math.radians(p_LngA)))
	l_yA = l_earthR *(math.cos(math.radians(p_LatA)) * math.sin(math.radians(p_LngA)))
	l_zA = l_earthR *(math.sin(math.radians(p_LatA)))

	l_xB = l_earthR *(math.cos(math.radians(p_LatB)) * math.cos(math.radians(p_LngB)))
	l_yB = l_earthR *(math.cos(math.radians(p_LatB)) * math.sin(math.radians(p_LngB)))
	l_zB = l_earthR *(math.sin(math.radians(p_LatB)))

	l_xC = l_earthR *(math.cos(math.radians(p_LatC)) * math.cos(math.radians(p_LngC)))
	l_yC = l_earthR *(math.cos(math.radians(p_LatC)) * math.sin(math.radians(p_LngC)))
	l_zC = l_earthR *(math.sin(math.radians(p_LatC)))

	l_P1 = numpy.array([l_xA, l_yA, l_zA])
	l_P2 = numpy.array([l_xB, l_yB, l_zB])
	l_P3 = numpy.array([l_xC, l_yC, l_zC])

	#from wikipedia
	#transform to get circle 1 at origin
	#transform to get circle 2 on x axis
	l_ex = (l_P2 - l_P1)/(numpy.linalg.norm(l_P2 - l_P1))
	l_i = numpy.dot(l_ex, l_P3 - l_P1)
	l_ey = (l_P3 - l_P1 - l_i*l_ex)/(numpy.linalg.norm(l_P3 - l_P1 - l_i*l_ex))
	l_ez = numpy.cross(l_ex,l_ey)
	l_d = numpy.linalg.norm(l_P2 - l_P1)
	l_j = numpy.dot(l_ey, l_P3 - l_P1)

	#from wikipedia
	#plug and chug using above values
	l_x = (pow(p_DistA,2) - pow(p_DistB,2) + pow(l_d,2))/(2*l_d)
	l_y = ((pow(p_DistA,2) - pow(p_DistC,2) + pow(l_i,2) + pow(l_j,2))/(2*l_j)) - ((l_i/l_j)*l_x)

	# only one case shown here
	l_z = numpy.sqrt(abs(pow(p_DistA,2) - pow(l_x,2) - pow(l_y,2)))

	#l_triPt is an array with ECEF x,y,z of trilateration point
	l_triPt = l_P1 + l_x*l_ex + l_y*l_ey + l_z*l_ez

	#convert back to lat/long from ECEF
	#convert to degrees
	l_lat = math.degrees(math.asin(l_triPt[2] / l_earthR))
	l_lng = math.degrees(math.atan2(l_triPt[1],l_triPt[0]))

	return {'lat':l_lat, 'lng':l_lng}

#example data
p_LatA = 50.37954
p_LngA = 30.85162
p_DistA = 918
p_LatB = 50.37896 
p_LngB = 30.85699
p_DistB = 564
p_LatC = 50.37803 
p_LngC = 30.86295
p_DistC = 322

print f_tralateration(	p_LatA, p_LngA, p_DistA, 
			p_LatB, p_LngB, p_DistB, 
			p_LatC, p_LngC, p_DistC
		     )
