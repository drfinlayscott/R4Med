# install_FLR_packages.R
# Copyright 2013 Finlay Scott and Chato Osio
# Maintainer: Finlay Scott, JRC, finlay.scott@jrc.ec.europa.eu
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#--------------------------------------------------------------------------
# This small script installs the latest FLR packages.
# They are for R >= 3.0

# Some of the FLR packages have dependencies on other non-FLR packages.
# So you should install:
# triangle
# copula
# plyr
# ggplot2
# using the normal install.packages() command

# To install FLR packages use the command:

# One package at a time
install.packages(pkgs = "FLCore", repos="http://flr-project.org/R")

# Or all of them together
install.packages(pkgs = c("FLCore","FLash","FLAssess","FLXSA","ggplotFL","FLBRP","FLa4a"), repos="http://flr-project.org/R")

# If you are not using RStudio you can just do:
install.packages(repos="http://flr-project.org/R")
# and select the packages from window that opens

# The packages to install:
# FLCore
# FLash
# FLAssess
# FLXSA
# ggplotFL
# FLBRP
# FLa4a

