%------------------------------------------------------------------------------%
%> @file    TransportPack.m
%> @brief   Main page documentation.
%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%
%> @mainpage MATLAB TransportPack
%>
%>
%> @section intro Introduction
%>
%> This package of MATLAB classes and functions provides a basic neutron
%> transport
%> capability in an easy-to-use environment with which many engineering students
%> are already familiar.  The purpose of this package is two-fold.  First, it is
%> to provide illustrative implementations of transport algorithms for use
%> in the classroom.  Second, it is to provide a simple framework where new
%> algorithms may be quickly implemented and tested.
%>
%>
%> @section problems Problems Solved
%> 
%> TransportPack solves a variety of neutron transport problems.  These include
%> both fixed source and eigenvalue problems in one, two, or three dimensions.
%> Currently, only the discrete ordinates method is implemented, but other
%> approaches may be added in the future (including other deterministic
%> approaches like the Method of Characteristics and stochastic (i.e. Monte 
%> Carlo) methods.
%>
%> @note One and three dimensional routines are under development.
%>
%> @section solver Solvers
%> 
%> For the deterministic problems, a main goal was to illustrate use of advanced
%> solution techniques.  These include Krylov solvers for within-group
%> iterations and various acceleration techniques for both source iteration and
%> eigenvalue (i.e. power) iteration.  In particular, for the within-group
%> ("inner") problem, the solvers available are
%> - source iteration
%> - %Livolant acceleration
%> - GMRES iteration
%> 
%> Preconditioners such as diffusion synthetic acceleration (DSA), coarse
%> mesh finite difference (CMFD) acceleration, and coarse mesh rebalance
%> (CMR) are possible future additions applicable to all inner iteration
%> solvers.
%>
%> For eigenvalue problems, currently only an unaccelerated power method 
%> is employed.  Future additions will be use of the Arnoldi method via
%> MATLAB's own "eigs" function, CMFD acceleration, and possible other
%> "low order" schemes.
%> 
%> @section other Structure and Documentation
%> 
%> This package is written completely in an object-oriented fashion, making
%> extensive use of MATLAB classes.  While this approach may be foreign to new 
%> (and old) users of MATLAB, we hope users will quickly see the benefits of
%> such an approach.  Chief examples include 
%> - encapsulation (Classes and methods can be used without worrying about their
%>   implementation)
%> - better data management (Many objects are created once and only once, and so
%>   data is not needlessly copied over and over again)
%>
%> Worth noting is the documenation, which is created automatically by
%> Doxygen.  However, Doxygen doesn't have support for MATLAB directly, and
%> so a tool is used that makes the .m files appear to be in c++ syntax.
%>
%-------------------------------------------------------------------------%
