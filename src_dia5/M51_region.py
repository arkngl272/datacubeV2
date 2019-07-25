#!/usr/bin/env python
# -*- coding: utf8 -*-

##############################################################################
#                                                                            #
#           This routine calculates Total mass (HI, H2 and Dust)             #
#                     and Gas-Dust Rate in a M51 region                      #
#                                                                            #
##############################################################################

Gasmolmass=3.05e+08
Gasatmass=7.39e+07
Dmass=3.31e+07


Tmass=1.36*(Gasmolmass+Gasatmass)+Dmass
Rate=1.36*(Gasmolmass+Gasatmass)/Dmass

print "Total Mass:","{:.2e}".format(Tmass),"[M"u'\u2609'"]"
print "Gas-Dust Rate:",Rate
