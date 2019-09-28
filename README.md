# sifm
Symbolic Information Flow Measurement

# sifm
Symbolic Information Flow Measurement

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Id: Readme.md,v 1.0 02-07-2019, IBU
#
#                This source code is part of
#
#   Symbolic Information Flow Measure Code for studying the
#   information flow in dynamical systems
#
#                        VERSION 1.0
#
# Written by Hiqmet Kamberaj.
# Copyright (C) 2019 Hiqmet Kamberaj.
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
#
## Embedded Parameters computations:
# Makefile_embd --- make file for compiling and creating the executable embdrun.exe (MPI version);
# Makefile_embd_ser --- make file for compiling and creating the executable embdrun_ser.exe (serial version);
# This class contains the following files:
#
# ---------------------------------------------------------------------------------------------------------------------------------
# demo_embd.f90:
# This file is initializing the EMBD_CLASS, which is the module used for calculations of the embedded parameters. It contains the subroutines/functions 
# used to read the input parameters for the computations of the embedded parameters, such as:
#         write(*,*) "1  - Number of Time Frames            - Nframes"
#         write(*,*) "2  - Number of time series            - Natoms"
#         write(*,*) "3  - Dimensionality of the problem    - Ndim"
#         write(*,*) "4  - Time lag minimum value           - T1"
#         write(*,*) "5  - Time lag maximum value           - T2"
#         write(*,*) "6  - Embedded dimension minimum value - M1"
#         write(*,*) "7  - Embedded dimension maximum value - M2"
#         write(*,*) "8  - Set debugging flag value         - Debug"
#
# This subroutine calls the function 'EMBD_DRIVER', which is part of the module EMBD_CLASS found in these fortran files:
# 1) embdmodule.f90       -- Parallel (using MPI protocols) of the EMBD_CLASS
# 2) embdmodule_ser.f90   -- Serial version of the EMBD_CLASS
# EMBD_CLASS contains the following functions:
# subroutine EMBD_DRIVER(Nf, Nd, Na, M1, M2, Tau1, tau2, db) -- it is the driver subroutine of the EMBD_CLASS
# subroutine MPI_Broadcast() -- it used to broad coast to all processors all attributes of the class EMBD_CLASS
# subroutine read_xyz() -- it is used to read the input time series representing the dynamical variables
# subroutine print_embdparam() -- It prints out the embedded dimension parameters
# subroutine Allocate_embeddedDim() -- It is used to allocate memory for the global dynamical variables of the EMBD_CLASS
# subroutine deAllocate_EmbeddedDim() -- It is used to free the memory allocated for the global dynamical variables of the EMBD_CLASS
# subroutine CompEmbeddedDimensions_fast() -- Computes the time shift and state vector dimension of time series - fast routine of MPI
# subroutine getTimeLag_fast(Idof, tau_arg) -- Calculates the time Lag using the fast routine - MPI version
# subroutine getOptEmbDim_fast(idof, Mopt_arg) -- Calculates the embedded dimension using the fast routine - MPI version
# subroutine CompEmbeddedDimensions_slow() -- Computes the time shift and state vector dimension of time series - slow routine of MPI
# subroutine getTimeLag_slow(N, x, tau_min, tau_max, tau_arg, MI) -- Calculates the time Lag using the slow routine - MPI version
# subroutine getOptEmbDim_slow(n, X, m1, m2, tau, Mopt, count_fnn) -- Calculates the embedded dimension using the slow routine - MPI version
# subroutine CompEmbeddedDimensions(gtau, gmopt) -- Computes the time shift and state vector dimension of time series - serial routine
# subroutine getTimeLag(N, x, tau_min, tau_max, tau_arg, MI) -- Calculates the time Lag using the serial routine
# subroutine getOptEmbDim(n, X, m1, m2, tau, Mopt, count_fnn)  -- Calculates the embedded dimension using the serial routine 
# -------------------------------------------------------------------------------------------------------------------------------------------------
#
## Symbolic analysis using the method by Kamberaj & van der Vaart :
# Makefile_symb --- make file for compiling and creating the executable symbrun.exe (MPI version);
# Makefile_symb_ser --- make file for compiling and creating the executable symbrun_ser.exe (serial version);
# This class contains the following files:
#----------------------------------------------------------------------------------------------------------------------------------------------------
# demo_symb.f90:
# This file is initializing the SYMB_CLASS, which is the module used for symbolization of the dynamical variables. It contains the subroutines/functions 
# used to read the input parameters for symbolization of the time series, such as:
#          write(*,*) "1  - Number of Time Frames                                                - Nframes"
#          write(*,*) "2  - Number of time series                                                - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                        - Ndim"
#          write(*,*) "4  - Set debugging flag value                                             - Debug"			
#          write(*,*) "5-   Symbolic method flag (1: 0s, 1s; 2: 0s, 1s, 2s, etc; 3: Monte Carlo) - qSymbolic"			
#          write(*,*) "6  - Set Nr of Monte Carlo Steps                                          - Nmc"		
# This subroutine calls the function 'SYMB_DRIVER', which is part of the module SYMB_CLASS found in these fortran subroutines/functions:
# 1) symbmodule.f90       -- Parallel (using MPI protocols) of the SYMB_CLASS
# 2) symbmodule_ser.f90   -- Serial version of the SYMB_CLASS
# SYMB_CLASS contains the following subroutines/functions:
# subroutine read_embdparam() -- it reads the embedded parameters calculated from the EMBD_CLASS functions
# subroutine write_xyzs() -- It prints out the symbolic Trajectories of dynamical variables
# subroutine Allocate_Symbolics() -- It is used to allocate memory for the global dynamical variables of the SYMB_CLASS
# subroutine deAllocate_Symbolics() -- It is used to free the memory allocated for the global dynamical variables of the SYMB_CLASS
# subroutine symbolize_trajectoryMC(Ndof, nframes, X, global_XS, tau, mopt, Nmc) -- Symbolizes the time series using Monte Carlo Method 1 (0,1,2,...)
#                                                                                -- Parallel and Serial versions 
# subroutine symbolize_trajectory(Natoms, Ndim, nframes, X, XS, topt, mopt) -- Symbolizes the time series using Method 2 (0,1,2, ...) (only serial routine) 
# subroutine symbolize(Natoms, Ndim, nframes, X, XS) -- Symbolizes the time series using Method 3 (0 and 1) (only serial routine)
# subroutine symbolic_entropy1D(ndata,xs,m,tau,H) -- Computes the Shannon entropy of one dimensional symbolic time series
# --------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Symbolic Transfer Entropy:
# Makefile_te --- make file for compiling and creating the executable te.exe (MPI version);
# Makefile_te_ser --- make file for compiling and creating the executable te_ser.exe (serial version);
# This class contains the following files:
#----------------------------------------------------------------------------------------------------------------------------------------------------
# demo_te.f90
# This file is initializing the TE_CLASS, which is the module used for calculation of transfer entropy. It contains the subroutines/functions 
# used to read the input parameters for computation of transfer entropy between the time series, such as:
#          write(*,'("Enter:")')
#          write(*,*) "1  - Length of time series                                                                   - Nframes"
#          write(*,*) "2  - Number of time series                                                                   - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                                           - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete (m=1); 2: Discrete (m>1); 3: Schreiber)  - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization                     - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling)                 - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                                      - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                                             - Rcut"
#          write(*,*) "9  - Confidence level for averages                                                           - statP"
#          write(*,*) "10 - Set debugging flag value                                                                - Debug"
# This subroutine calls one of the following functions (which are part of the TE_CLASS in mte.f90 (MPI version) and mte_ser.f90 (serial version):
#     IF (qTEMethod == 1) THEN
#          CALL STE1_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
#      ELSEIF (qTEMethod == 2) THEN
#          CALL STE2_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
#      ELSEIF (qTEMethod == 3) THEN
#          CALL STE3_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
#      ELSE
#          CALL TE_DRIVER(Nframes, Ndim, Natoms, qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)      
#      ENDIF
# TE_CLASS contains the following subroutines/functions:
# subroutine write_TE(qTeShuffle, model) -- It prints out the transfer entropy
# subroutine Allocate_TransferEntropy(qteShuffle, model) -- It allocates memory for global dynamical variables of TE_CLASS
# subroutine deAllocate_TransferEntropy() -- It frees allocated memory for global dynamical variables of TE_CLASS
# subroutine MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, model) -- It broad cast to all processors global variables of TE_CLASS (MPI version)
# subroutine getCrossSTransferEntropyNDIM_MPIDOF(qteNorm, debug) -- It computes symbolic TE between two dynamical variables X(i,1,1:d), X(i,2,i:d), where
#                                                                -- i = number of pairs; d = dimensionality of the problem - MPI version
#                                                                -- No shuffling applies
#                                                                -- using method by Kamberaj & Van der Vaart
# subroutine getSTransferEntropyNDIM_MPIDOF(qteNorm, debug) -- It computes symbolic TE between dynamical variables presented as matrix X(1:N,1:d) 
#                                                           -- using method by Kamberaj & Van der Vaart
#                                                           -- MPI version
#                                                           -- No shuffling applies
# subroutine getSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP) -- It computes symbolic TE between dynamical variables presented as matrix X(1:N,1:d) 
#                                                           -- using method by Kamberaj & Van der Vaart
#                                                           -- MPI version
#                                                           -- Shuffling applies
# subroutine getTransferEntropyNDIM_MPIDOF(qteNorm, debug) -- It computes symbolic TE between dynamical variables presented as matrix X(1:N,1:d) 
#                                                           -- using Schreiber Method
#                                                           -- MPI version
#                                                           -- No shuffling applies
# subroutine getTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP) -- It computes symbolic TE between dynamical variables presented as matrix X(1:N,1:d) 
#                                                           -- using Schreiber Method
#                                                           -- MPI version
#                                                           -- Shuffling applies
#
# Some auxiliary functions / subroutines:
# symbolic_TE_entropies1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy) -- It computes symbolic transfer entropies between two symbolic strings time series x and y 
#                                                                            -- using method by Kamberaj & Van der Vaart
# subroutine symbolic_TE_entropiesND(ndata,ndim,xs,ys,m1,m2,tau1,tau2, Txy,Tyx,Hxx,Hyy,hx,hy) --- The same as above but for Ndim symbolic time series 
#                                                                                             --- using method by Kamberaj & Van der Vaart
# subroutine symbolic_TE_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy) -- It computes only Txy symbolic transfer entropy between two symbolic strings time series x and y 
#                                                                   -- using method by Kamberaj & Van der Vaart
# subroutine symbolic_TE_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy) - The same as above but for Ndim symbolic time series using method by Kamberaj & Van der Vaart
# subroutine symbolic_TE2_entropies1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy) -- It computes symbolic transfer entropies between two symbolic strings time series x and y 
#                                                                                        -- using Schreiber method
# subroutine symbolic_TE2_entropiesND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy) - The same as above but for Ndim real variable time series
# subroutine symbolic_TE2_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy) -- It computes only Txy symbolic transfer entropy between two symbolic strings time series x and y 
#                                                                    -- using Schreiber method
# subroutine symbolic_TE2_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy) - The same as above but for Ndim real variable time series
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Local Symbolic Transfer Entropy:
# Makefile_lte --- make file for compiling and creating the executable lte.exe (MPI version);
# Makefile_lte_ser --- make file for compiling and creating the executable lte_ser.exe (serial version);
# This class contains the following files:
#----------------------------------------------------------------------------------------------------------------------------------------------------
# demo_lte.f90
# This file is initializing the LTE_CLASS, which is the module used for calculation of local transfer entropy. It contains the subroutines/functions 
# used to read the input parameters for computation of local transfer entropy between the time series, such as:
#          write(*,*) "1  - Length of time series                                                                   - Nframes"
#          write(*,*) "2  - Number of time series                                                                   - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                                           - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete (m=1); 2: Discrete (m>1); 3: Schreiber)  - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization                     - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling)                 - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                                      - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                                             - Rcut"
#          write(*,*) "9  - Confidence level for averages                                                           - statP"
#          write(*,*) "10 - Set debugging flag value                                                                 - Debug"
# This subroutine calls one of the following functions (which are part of the LTE_CLASS in localmte.f90 (MPI version) and localmte_ser.f90 (serial version):
#     IF (qTEMethod == 1) THEN
#          CALL LSTE_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP) -- driver for symbolic local transfer entropy
#      ELSE
#          CALL LTE_DRIVER(Nframes, Ndim, Natoms, qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP) -- driver for symbolic local transfer entropy using Schreiber method     
#      ENDIF
# subroutine WRITE_LTE(qTEShuffle) -- Writes on output the computed local transfer entropies
# subroutine Allocate_LocalTransferEntropy(qTeShuffle, model) -- It allocates memory for global dynamical variables of the LTE_CLASS
# subroutine deAllocate_LocalTransferEntropy() -- It frees the allocated memory for the global dynamical variables of the LTE_CLASS
# subroutine MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, model) -- It broad casts the global variables of LTE_CLASS to all processors (MPI version)
# subroutine getLocalSTransferEntropyNDIM_MPIDOF(qteNorm, debug) -- It computes local symbolic transfer entropy without shuffling (MPI version) using method by Kamberaj & Van der Vaart
# subroutine getLocalSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP) -- It computes local symbolic transfer entropy with shuffling (MPI version) 
#                                                                                                  -- using method by Kamberaj & Van der Vaart
# subroutine getLocalTransferEntropyNDIM_MPIDOF(qteNorm, debug) -- It computes local symbolic transfer entropy without shuffling (MPI version) using Schreiber method
# subroutine getLocalTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP) -- It computes local symbolic transfer entropy with shuffling (MPI version) 
#                                                                                                 -- using Schreiber method
# 
# Auxiliary functions/subroutines:
# subroutine symbolic_LTE2_entropies(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy,LTyx) -- It computes local symbolic transfer entropies (LTxy, LTyx) between two symbolic time series xs, ys 
#                                                                                -- using Schreiber method
# subroutine symbolic_LTE2_entropy(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy) -- It computes local symbolic transfer entropy (LTxy) between two symbolic time series xs, ys 
#                                                                         -- using Schreiber method
# subroutine symbolic_LTE_entropies(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy,LTyx) -- It computes local symbolic transfer entropies (LTxy, LTyx) between two symbolic time series xs, ys 
#                                                                               -- using method by Kamberaj & Van der Vaart
# subroutine symbolic_LTE_entropy(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy) -- It computes local symbolic transfer entropy (LTxy) between two symbolic time series xs, ys 
#                                                                        -- using method by Kamberaj & Van der Vaart
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Symbolic Mutual Information:
# Makefile_mi --- make file for compiling and creating the executable mirun.exe (MPI version);
# Makefile_mi_ser --- make file for compiling and creating the executable mirun_ser.exe (serial version);
# This class contains the following files:
#----------------------------------------------------------------------------------------------------------------------------------------------------
# demo_mi.f90
# This file is initializing the MI_CLASS, which is the module used for calculation of symbolic mutual information between time series. It contains the subroutines/functions 
# used to read the input parameters for computation of the symbolic mutual information, such as:
#          write(*,*) "1  - Length of time series                                                  - Nframes"
#          write(*,*) "2  - Number of time series                                                  - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                          - Ndim"
#          write(*,*) "4  - Flag for Method of MI calculation (1: Discrete; 2: Schreiber)          - qMIMethod"
#          write(*,*) "5  - Flag for Method of Shuffling of MI (1: Permutation; 2: Block shuffling)- qMIShuffle"
#          write(*,*) "6  - Number of Shuffling of time series                                     - Nshuffles"
#          write(*,*) "7  - Cutoff for Mutual Information minimum value                            - Rcut"
#          write(*,*) "8  - Confidence level for averages                                          - statP"
#          write(*,*) "9  - Set debugging flag value                                               - Debug"
# This subroutine calls one of the following functions (which are part of the MI_CLASS in minfo.f90 (MPI version) and minfo_ser.f90 (serial version):
#      IF (qMIMethod == 1) THEN
#          CALL SMI_DRIVER(Nframes, Ndim, Natoms, qMIShuffle, debug, Nshuffles, Rcut, StatP) -- It computes symbolic MI using method by Kamberaj & Van der Vaart
#      ELSE
#          CALL  MI_DRIVER(Nframes, Ndim, Natoms, qMIShuffle, debug, Nshuffles, Rcut, StatP) -- It computes symbolic MI using method by Schreiber    
#      ENDIF
# subroutine write_MI(qMIShuffle) -- It writes the symbolic MI in an output file
# subroutine Allocate_MI(qmimethod, qMIShuffle) -- It allocates memory for global dynamical variables of the MI_CLASS
# subroutine deAllocate_MI() -- It frees memory allocated for global dynamical variables of the MI_CLASS
# subroutine MPI_BroadCast(qMIShuffle, debug, Nshuffles, r0, statP, model) -- It broad casts the global variables of the MI_CLASS to all processors (MPI version)
# subroutine getSMINDIM_MPIDOF(debug) -- It computes symbolic MI using the method by Kamberaj & van der Vaart (MPI version) without shuffling
# subroutine getSMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statP) -- It computes symbolic MI using the method by Kamberaj & van der Vaart (MPI version) with shuffling
# subroutine getMINDIM_MPIDOF(debug) -- It computes symbolic MI using the method by Schreiber (MPI version) without shuffling
# subroutine getMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statP) -- It computes symbolic MI using the method by Schreiber (MPI version) with shuffling
#
# Auxiliary functions / subroutines:
# subroutine symbolic_MI_ND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,mi) -- It computes symbolic MI between two symbolic time series xs, ys using method by Kamberaj & van der Vaart
# subroutine symbolic_MI2_ND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,MI) -- It computes symbolic MI between two symbolic time series xs, ys using method by Schreiber
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## System configuration setups:
# Makefile_setup --- make file for compiling and creating the executable sysrun.exe;
# This class contains the following files:
#----------------------------------------------------------------------------------------------------------------------------------------------------
# demo_setup.f90
# This file is initializing the SYSTEM_CLASS, which is the module used for reading/creating time series for the benchmark cases. It contains the subroutines/functions 
# used to read/create dynamical variables for the case studies, such as:
#            write(*,*) "1  - Time to start trajectory     - iStart"
#            write(*,*) "2  - Number of steps to skip      - iSkip"
#            write(*,*) "3  - Time to stop trajectory      - iStop"
#            write(*,*) "4  - Offset residue numbering     - Offset"
#            write(*,*) "5  - Number of time series        - Natoms"
#            write(*,*) "6  - Dimensionality of the problem - Ndim"
#            write(*,*) "7  - Set debugging flag value      - Debug"				
#            write(*,*) "8  - Set Benchmark flag value      - Benchmark"
#
# This subroutine calls one of the following functions (which are part of the SYSTEM_CLASS in syssetup.f90:
#      IF (benchmark == 1) then
#          CALL SETUP_DRIVER(Nframes, Ndim, Natoms, Offset, iStart, iSkip, iStop, debug, benchmark, Ax, Cxy, Sigma_x, Sigma_y)
#      else
#          CALL SETUP_DRIVER(Nframes, Ndim, Natoms, Offset, iStart, iSkip, iStop, debug, benchmark)
#      endif
# subroutine mdTrajEnergy(Istart, iSkip, iStop) --- It reads the energy from a Molecular dynamics simulations for each frame
# subroutine mdTrajPDBVmd(Istart, iSkip, iStop) --- It reads the molecular dynamics trajectory in PDB format as created using VMD software, each residue is represented by 
#                            --- CA atom coordinates, or
#                            --- Collective coordinates representing essential degrees of freedom for each amino acid computed using Machine Learning Approach 
# subroutine Benchmark1(Ax,Cxy,Sigmax,Sigmay,iStart,iSkip,iStop) --- Benchmark 1
!! The Multivariate Gaussian Process
!!
# subroutine Benchmark2(Istart, Iskip, Istop):
!!  Mutual information in this Experiment:
!!  Y = X+Z
!!  X uniformly distributed in [-1/2:1/2]
!!  Z uniformly distributed in [-a/2:a/2]
!!  f_y(y) = convolution(X,Z) 
!!  I(X;Y) = h(Y) - h(Z)
!!  In Nats: I(X;Y) = a/2  - ln(a) if a <= 1; otherwise I(X;Y)=1/(2*a) if a>=1.
!!  In Bits: I(X;Y) = [a/2 - ln(a)]/ln(2) if a <= 1; otherwise I(X;Y)=[1/(2*a)]/ln(2) if a>=1.
!!  Example: a=0.5, then I(X;Y) = 1.36067 (bits); h(Y)=0.49075 (bits)
!!
# subroutine Benchmark3(Istart, Iskip, Istop)
!!  Mutual information in the this Experiment:
!!  Y = X+Z
!!  X  uniformly distributed in [-1/2:1/2]
!!  Z1 uniformly distributed in [-a/2:a/2]
!!  Z2 uniformly distributed in [-a/4:a/4]
!!  f_y(y) = convolution(X,Z1,Z2)
!!
# subroutine Benchmark4(Istart, Iskip, Istop)
!!
!!  Mutual information in the this Experiment:
!!  Throwing a fair die:
!!  X -> Variable representing the value of the face closeer to us
!!  Y -> Variable representing the value of the top of the die
!!  Possible values are 1, 2, 3, 4, 5, 6
!!  Theoretical value:
!!  I(T,B) = H(B) - H(B|T) = log(6) - log(4) = log(3) - 1
!!  H(B) = 2.584962501; H(B|T) = 2.0
!!
# subroutine Benchmark5(Istart, Iskip, Istop)
!!
!!  Mutual information in the this Experiment:
!!  Flipping a fair Coin:
!!  X -> Variable representing the value of the bottom of the coin
!!  Y -> Variable representing the value of the top of the coin
!!  Head = 1 and Tail = 0
!!  Theoretical value:
!!  I(T,B) = H(B) - H(B|T) = log(2) - 0 = 1
!!
# subroutine write_xyz() --- It writes in output the real variable trajectory coordinates
# subroutine write_xyzs() --- It writes in output the symbolic trajectory coordinates
# subroutine read_xyz -- It reads trajectory coordinates
# subroutine Allocate_System(Benchmark) --- It allocates memory for global variables of the SYSTEM_CLASS
# subroutine deAllocate_System() --- It frees allocated memory for global variables of the SYSTEM_CLASS
#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Auxiliary Program functions/subroutines of the class TEUTILS_CLASS
# This class contains the following routines:
#
# function new_unit()  result (result) -- Returns a new unit for a new file to be open
# subroutine iswap(x,y) -- Swaps the values of two integers
# subroutine iswap_vec(D,x,y) -- Swaps the values of two integer arrays of size D
# subroutine rswap(x,y) -- Swaps the values between two real variables
# subroutine rswap_vec(x,y) -- Swaps the values between two real arrays of size D
# Interface Swap --- is an interface for iswap, iswap_vec, rswap, rswap_vec
# subroutine rperm(N0, N, p, IS) --- performs permutation of the array elements p from N0 to N
# subroutine BlockShuffle1D(N,x,m,tau,IS) --- performs block shuffling of the one-dimensional integer array X
# subroutine BlockShuffleND(N,ndim,x,m,tau,IS) ---  performs block shuffling of the N-dimensional integer array X
# subroutine BlockShuffleND2(N,ndim,x,m,tau,IS) --- performs block shuffling of the N-dimensional real array X
# subroutine BlockShuffle1D2(N,x,m,tau,IS) --- performs block shuffling of the one-dimensional real array X
# subroutine shuffleND(n,ndim,x,is) --- performs simple shuffling of the N-dimensional integer array X
# subroutine shuffle1D(n,x,is) --- performs simple shuffling of the one-dimensional integer array X
# subroutine shuffleND2(n,ndim,x,is) --- performs simple shuffling of the N-dimensional real array X
# subroutine shuffle1D2(n,x,is) --- performs simple shuffling of the one-dimensional real array X
# FUNCTION psi(xx) RESULT(fn_val) -- computes the derivative of the gamma function
# subroutine addStrings_slow(a,b) -- Slow (depending on the machine) routine of adding two strings
# subroutine addStrings(a,b) --- Fast (depending on the machine) routine of adding two strings
# SUBROUTINE addStrings2(ST,STMAX,STLEN,ADST,ADLEN) --- Fast (depending on the machine) routine of adding two strings
# function toString(a) result(b) - it converts integer to string
# function toInteger(a) result(b) -- it converts string to integer
# subroutine UPCASE (STRING) --- it converts a string to upper case letters
# subroutine DOWNCASE (STRING) --- it converts a string to lower case letters
# function mean(X,n) result(r) --- computer the mean value of the elements of a vector
# function rmin1d(X) --- Finds the minimum of an array of real numbers
# function rmax1d(X) --- Finds the maximum of an array of real  numbers
# function imax1d(X) result(r) --- Finds the maximum of an array of integer numbers
# function imin1d(X) result(r) --- Finds the minimum of an array of integer numbers
# subroutine sort(n,arr,indx) --- It sorts elements of an array (machine depending)
# subroutine sort2(n,r0,indx0,indx) --- It sorts elements of an array (machine depending)
# function vcopy(Xin,istart,M,istep,Sign) result(Xou) --- it copies elements of an array to another
# function getCovariance(Nframes,X, Y) result(Cov) --- it computes covariance matrix
# function find(i,list) result(r) --- it find the element i in the list of elements
# Subroutine statTest(N, r0, val_in, val_out, stat_prob) -- performs a statistical test
# FUNCTION erfcc_s(x); FUNCTION erfcc_v(x) --- it computes the error function value (with interface erfcc)
# function square(x) result(r) -- computes the square of a number
# Subroutine TransferFunctionInn_2d(z, r, ifunc) --- transfer function of Machine Learning (two-dimensional array)
# subroutine TransferFunctionInn_1d(Z, R, ifunc) --- transfer function of Machine Learning (one-dimensional array)
# function TransferFunctionInn_0d(z, ifunc) result(r) --- definition of the transfer function for Machine Learning
# Subroutine TransferFunctionInn_derivative_2d(z, r, ifunc) -- derivative of transfer function of Machine Learning (two-dimensional array)
# subroutine TransferFunctionInn_derivative_1d(Z, R, ifunc) --- derivative of transfer function of Machine Learning (one-dimensional array)
# function TransferFunctionInn_derivative_0d(z, ifunc) result(r) --- definition of the derivative of transfer function for Machine Learning
# subroutine getTimeLag(A_n, Tau) --- compute the time lag from the correlation function
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Auxiliary Program functions/subroutines of the class TELinkList_class
# This class contains the following routines:
#
# subroutine LinkList(head,tail,temp) --- it builds up the link list
# function compareStrings(a,b) result(r) --- a function to compare two strings if they are equal or not
# subroutine LocalLinkList(head,tail,temp,visited_state) --- builds up a local link list for each visited state
# subroutine freeList(head) --- deletes the link list and its heap memory storage
# DeleteList(head) -- it just deletes the link list 
# function SymbolicShannonEntropy(head) result(Entropy) --- it computes the Shannon entropy of symbolic states 
# subroutine SymbolicJointProbabilityTraj(head, total, prob, frequency) -- It computes the joint prob distribution of symbolic trajectory of states
# function getTotalStates(head) result(total) -- gets the total number of visited states
# function CountStates(head) result(Nstates) -- the same as 'getTotalStates'
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Program Control Variables in the class SIFM_KINDS
# This class contains the following definitions in the file 'defs.f90':
! Kind constants for system-independent precision specification.
! H. Kamberaj
!    integer, parameter :: sifm_real4 = kind(0.)   !-- real single precision
!    integer, parameter :: sifm_real8 = kind(0.d0) !-- real double precision
!
!-- Kinds for specified integer ranges.
!    integer, parameter :: sifm_int1 = selected_int_kind(2)  !-- 99 max
!    integer, parameter :: sifm_int2 = selected_int_kind(4)  !-- 9,999 max
!    integer, parameter :: sifm_int4 = selected_int_kind(8)  !-- 99,999,999 max
!    integer, parameter :: sifm_int8 = selected_int_kind(16) !-- 9,999,999,999,999,999 max
!
!-- Kind for working real precision.
!    integer, parameter :: sifm_real = sifm_real8
!-- Big constants. These are less than huge so that we have some
!-- room for comparisons and other arithmetic without overflowing.
!    real(sifm_real), parameter :: rzero = 0.0_sifm_real
!    real(sifm_real), parameter :: one  = 1.0_sifm_real
!    real(sifm_real), parameter :: toBits  = one/log(2.0_sifm_real)
!    real(sifm_real), parameter :: half = 0.5_sifm_real
!    real(sifm_real), parameter :: M_PI = 3.141592653589793238462643383279502884197_sifm_real
!    real(sifm_real), parameter :: M_2PI = 6.283185307179586476925286766559005768394_sifm_real
!    real(sifm_real), parameter :: ToDeg = 180.0_sifm_real/M_PI
!    real(sifm_real), parameter :: LOG_ZERO=-30.0_sifm_real
!    real(sifm_real), parameter :: max_exp = maxexponent(rzero)
!    real(sifm_real), parameter :: min_exp = minexponent(rzero)
!    real(sifm_real), parameter :: kB = 0.00198614_sifm_real
!    REAL(sifm_real), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sifm_real
!    real(sifm_real), parameter :: PRTINY = tiny(rzero) * 16.0_sifm_real
!    real(sifm_real), parameter :: PRBIG  = huge(sifm_real) / 16.0_sifm_real
!    integer, parameter         :: PIBIG  = huge(sifm_int8)/16
!    character(len=1), dimension(26) :: Symbol = (/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',  &
!                                                  'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',  &
!                                                  's', 't', 'u', 'v', 'w', 'x', 'y', 'z' /)
#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Program Control Variables in the class TEMPI_CLASS for the MPI version
# This class contains the following definitions in the file 'defs_mpi.f90':
!type val_t
!    real(sifm_real) :: value
!    integer :: rank
!end type val_t
!      integer, save :: NumProcs = 1
!      integer, save :: MyID = 0
!      integer, parameter :: master = 0
!      integer , parameter :: msgtag2 = 12
!      integer, save, dimension ( MPI_STATUS_SIZE ) :: status
!contains
# subroutine MPI_start() --- Start MPI communications
# subroutine MPI_finish() --- Finalizes the MPI communications
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random number generator routines in the class TERANDOM_CLASS
# This class contains the following definitions in the file 'random.f90':
# Contains:
# function RANDF(ig) result(random_rtn) --- Returns a double precision uniformly distributed random number in (0,1) from the ig-th stream
# subroutine  SetSeed (ig,s) -- sets seed value s for stream ig
# FUNCTION GAUSS_rand(g) result(RND) --- generates a Gaussian random deviate of 0.0 mean and standard deviation 1 for the stream g\\
# subroutine rndnum_iniall(nrnd, rngc, mynod) --- it initialize all needed for random number generators:
#                                             --- whether a sequence of Gaussians or u(0,1)
#                                             --- number of streams
#                                             --- and the processor ID (MPI version) 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
## BSANN folder contains python codes for performing Machine Learning to encode the information from the real dynamical variables.
# The following files are included:
# driver_bsann_reducedim.py -- it contains the functions to test benchmark data:
# 1) myTest1("data/benchmark1/lorenzmapode45.dat",3) -- accepts as input (X,Y,Z) dynamical variables and encodes this information into the output Z 
#                                                       Z represents reduced dimension space of the Lorenz map system.
# 2) myTest2("data/benchmark2/coords.pdb",3) -- represents the case of the fragment C2 of the protein G generated using Molecular Dynamics Simulations;
#                                            Inputs  -- atomic coordinates in PDB format
#                                            Outputs -- encoded variable (one-dimension array) for each amino acid.
# nnclass.py -- python code for performing BSANN optimization
# trainingData_mdsim.py --- python code for manipulating input molecular dynamics simulations data.
#
# Auxiliary Libraries:
# Numpy is needed to run the python codes.
#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set up system configuration file
# sysrun.exe 10001 1 20000 0 2 1 1 1 0.0 2.0 1.0 1.0
#            write(*,*) "1  - Time to start trajectory      - iStart"
#            write(*,*) "2  - Number of steps to skip       - iSkip"
#            write(*,*) "3  - Time to stop trajectory       - iStop"
#            write(*,*) "4  - Offset residue numbering      - Offset"
#            write(*,*) "5  - Number of time series         - Natoms"
#            write(*,*) "6  - Dimensionality of the problem - Ndim"
#            write(*,*) "7  - Set debugging flag value      - Debug"				
#            write(*,*) "8  - Set Benchmark flag value      - Benchmark"			
#
#
# mpirun -np 2 embdrun.exe 10000 2 1 1 1 1 10 1
#         write(*,*) "1  - Number of Time Frames            - Nframes"
#         write(*,*) "2  - Number of time series            - Natoms"
#         write(*,*) "3  - Dimensionality of the problem    - Ndim"
#         write(*,*) "4  - Time lag minimum value           - T1"
#         write(*,*) "5  - Time lag maximum value           - T2"
#         write(*,*) "6  - Embedded dimension minimum value - M1"
#         write(*,*) "7  - Embedded dimension maximum value - M2"
#         write(*,*) "8  - Set debugging flag value         - Debug"
#
#
# mpirun -np 2 symbrun.exe 10000 2 1 1 1 1000 
#          write(*,*) "1  - Number of Time Frames                                                - Nframes"
#          write(*,*) "2  - Number of time series                                                - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                        - Ndim"
#          write(*,*) "4  - Set debugging flag value                                             - Debug"			
#          write(*,*) "5-   Symbolic method flag (1: 0s, 1s; 2: 0s, 1s, 2s, etc; 3: Monte Carlo) - qSymbolic"			
#          write(*,*) "6  - Set Nr of Monte Carlo Steps                                          - Nmc"		
#
# mpirun -np 2 lte.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#          write(*,*) "1  - Length of time series                                                                   - Nframes"
#          write(*,*) "2  - Number of time series                                                                   - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                                           - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete (m=1); 2: Discrete (m>1); 3: Schreiber)  - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization                     - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling)                 - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                                      - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                                             - Rcut"
#          write(*,*) "9  - Confidence level for averages                                                           - statP"
#          write(*,*) "10 - Set debugging flag value                                                                - Debug"
#
# mpirun -np 2 te.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#          write(*,*) "1  - Length of time series                                                  - Nframes"
#          write(*,*) "2  - Number of time series                                                  - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                          - Ndim"
#          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete; 2: Schreiber)          - qTEMethod"
#          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization    - qTENorm"
#          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling - qTEShuffle"
#          write(*,*) "7  - Number of Shuffling of time series                                     - Nshuffles"
#          write(*,*) "8  - Cutoff for Mutual Information minimum value                            - Rcut"
#          write(*,*) "9  - Confidence level for averages                                          - statP"
#          write(*,*) "10 - Set debugging flag value                                               - Debug"
#
# 
# mpirun -np 2 mirun.exe 10000 2 1 1 0 1 10 0.02 0.95 1
#          write(*,*) "1  - Length of time series                                                  - Nframes"
#          write(*,*) "2  - Number of time series                                                  - Natoms"
#          write(*,*) "3  - Dimensionality of the problem                                          - Ndim"
#          write(*,*) "4  - Flag for Method of MI calculation (1: Discrete; 2: Schreiber)          - qMIMethod"
#          write(*,*) "5  - Flag for Method of Shuffling of MI (1: Permutation; 2: Block shuffling)- qMIShuffle"
#          write(*,*) "6  - Number of Shuffling of time series                                     - Nshuffles"
#          write(*,*) "7  - Cutoff for Mutual Information minimum value                            - Rcut"
#          write(*,*) "8  - Confidence level for averages                                          - statP"
#          write(*,*) "9  - Set debugging flag value                                               - Debug"
#

############################################################################################################################# 
######################################################################################################################### 
