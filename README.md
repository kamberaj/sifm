# sifm
Symbolic Information Flow Measurement
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Id: Readme.txt,v 1.0 19-03-2018, IBU
#
#                This source code is part of
#
#   Symbolic Information Flow Measure Code for styding the
#   information flow in dynamical systems
#
#                        VERSION 1.0
#
# Written by Hiqmet Kamberaj.
# Copyright (C) 2018 Hiqmet Kamberaj.
# Check out h.kamberaj@gmail.com for more information.
#
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software Foundation; 
# GPL-3.0
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program; 
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Set up system configuration file
# sysrun.exe 10001 1 20000 0 2 1 1 1 0.0 2.0 1.0 1.0
#            write(*,*) "1  - Time to start trajectory     - iStart"
#            write(*,*) "2  - Number of steps to skip      - iSkip"
#            write(*,*) "3  - Time to stop trajectory      - iStop"
#            write(*,*) "4  - Offset residue numbering     - Offset"
#            write(*,*) "5  - Number of time series        - Natoms"
#            write(*,*) "6  - Dimenionality of the problem - Ndim"
#            write(*,*) "7  - Set debuging flag value      - Debug"				
#            write(*,*) "8  - Set Bechmark flag value      - Benchmark"			
#
#
# mpirun -np 2 embdrun.exe 10000 2 1 1 1 1 10 1
#         write(*,*) "1  - Number of Time Frames            - Nframes"
#         write(*,*) "2  - Number of time series            - Natoms"
#         write(*,*) "3  - Dimenionality of the problem     - Ndim"
#         write(*,*) "4  - Time lag minimum value           - T1"
#         write(*,*) "5  - Time lag maximum value           - T2"
#         write(*,*) "6  - Embedded dimension minimum value - M1"
#         write(*,*) "7  - Embedded dimension maximum value - M2"
#         write(*,*) "8  - Set debuging flag value          - Debug"
#
#
# mpirun -np 2 symbrun.exe 10000 2 1 1 1 1000 
#          write(*,*) "1  - Number of Time Frames                                                - Nframes"
#          write(*,*) "2  - Number of time series                                                - Natoms"
#          write(*,*) "3  - Dimenionality of the problem                                         - Ndim"
#          write(*,*) "4  - Set debuging flag value                                              - Debug"			
#          write(*,*) "5-   Symbolic method flag (1: 0s, 1s; 2: 0s, 1s, 2s, etc; 3: Monte Carlo) - qSymbolic"			
#          write(*,*) "6  - Set Nr of Monte Carlo Steps                                          - Nmc"		
#
# mpirun -np 2 lterun.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#
#          write(*,*) "1  - Length of time series                                                                   - Nframes"
#          write(*,*) "2  - Number of time series                                                                   - Natoms"
#          write(*,*) "3  - Dimenionality of the problem                                                            - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete (m=1); 2: Discrete (m>1); 3: Schreiber)  - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization                     - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling)                 - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                                      - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                                             - Rcut"
#          write(*,*) "9  - Confidence level for averages                                                           - statP"
#          write(*,*) "10 - Set debuging flag value                                                                 - Debug"
#
# mpirun -np 2 lte.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#          write(*,*) "1  - Length of time series                                                  - Nframes"
#          write(*,*) "2  - Number of time series                                                  - Natoms"
#          write(*,*) "3  - Dimenionality of the problem                                           - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete; 2: Schreiber)          - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization    - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                     - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                            - Rcut"
#          write(*,*) "9  - Confidence level for averages                                          - statP"
#          write(*,*) "10 - Set debuging flag value                                                - Debug"
#
# 
# mpirun -np 2 mirun.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#          write(*,*) "1  - Length of time series                                                  - Nframes"
#          write(*,*) "2  - Number of time series                                                  - Natoms"
#          write(*,*) "3  - Dimenionality of the problem                                           - Ndim"
#          write(*,*) "4  - Flag for Method of MI calculation (1: Discrete; 2: Schreiber)          - qMIMethod"
#          write(*,*) "5  - Flag for Method of Shuffling of MI (1: Permutation; 2: Block shuffling)- qMIShuffle"
#          write(*,*) "6  - Number of Shuffling of time series                                     - Nshuffles"
#          write(*,*) "7  - Cutoff for Mutual Information minimum value                            - Rcut"
#          write(*,*) "8  - Confidence level for averages                                          - statP"
#          write(*,*) "9  - Set debuging flag value                                                - Debug"
#

############################################################################################################################# 
