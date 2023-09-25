                       Prolate_swf

  Profcn is available as both a subroutine version provided as
  the module Prolate_swf and a stand alone version profcn. It was
  originally developed by arnie lee van buren and jeffrey boisvert
  over 20 years ago and has been improved a number of times since
  then.

  Table of Contents
  1. Purpose
  2. Introduction
  3. Input and Output
  4. Accuracy of results
  5. Obtaining the expansion d coefficients
  6. Obtaining the eigenvalues


 1. Purpose

  To calculate the first and second kind prolate radial functions r1
  and r2 and their first derivatives r1d and r2d for a given order m,
  a range of degrees l beginning at m and for a specific size parameter c
  and shape parameter x. To calculate the first kind prolate angular
  functions and their first derivatives with respect to eta for the same
  values of m, l and c and for a set of values of the angular coordinate
  eta. The subroutine version prolate_swf calculates values for a single
  input value of m . The stand alone version profcn calculates values for a
  range of values of m.

 2. Introduction

  Profcn is written in free format fortran. It is designed around the
  maximum number of decimal digits ndec and the maximum exponent nex
  available in real arithmetic. Procedures used in profcn allow for
  exponents much larger than nex since the resulting floating point
  radial function values are given as a characteristic and an integer
  exponent.

  Profcn can be run in either double precision or quadruple precision
  arithmetic. The choice is set in the module param provided in the github
  repository. If this is not available, then create param as follows:

    module param
    integer, parameter :: knd = selected_real_kind(8)
    logical, parameter :: debug = .true.
    logical, parameter :: warn = .true.
    logical, parameter :: output = .true.
    end module param

  Set the value of knd in the parenthesis to either 8 for double
  precision or 16 for quadruple precision arithmetic. Some compilers
  require that param be compiled prior to rather than after profcn.
  The logicals in param are described below in the discussion of the
  output files.

  Some computers may have more than 8 bytes for double precision
  data and more than 16 bytes for quadruple precision data or may use
  kind values that do not correspond to the number of bytes. In this
  case just use the appropriate integers for the kind parameters in
  module param. Also change the values of kindd and kindq set in
  statement 5 below below the comments section to the kind values for
  double precision data and quadruple precision data, respectively.

  The methods used in profcn are documented in the journal articles:
  A. L. Van Buren and J. E. Boisvert, "Accurate calculation of
  prolate spheroidal radial functions of the first kind and their
  first derivatives," Quart. Appl. Math. vol 60, pp. 589-599 (2002)

  A. L. Van Buren and J. E. Boisvert, "Improved calculation of prolate
  spheroidal radial functions of the second kind and their first
  derivatives," Quart. Appl. Math. vol 62, pp. 493-507 (2004)
  [available as a pdf file on the website].
  
  [Both articles available at github.com/mathieuandspheroidalwavefunctions]
  
  Expressions to calculate the radial functions of the first kind and its
  first derivatives for x = 1.0 when m = 0 are obtained from the limits of
  equations given in the first reference. The resulting function values
  provided by profcn are nearly fully accurate. The radial functions of
  the first kind and their first derivatives are equal to 0.0 when m is
  unequal to 0. [The corresponding radial functions of the second kind and
  their first derivatives are infinite for all m.]  

 3. Input and Output

  Following is a description of the input and output parameters in the
  call statement for the subroutine version. After that will be a
  description of the the input and output files associated with the
  stand alone version. Note that these output files, if desired, can
  also be obtained when running the subroutine version. See comments about
  this below.

  A sample input and resulting output from profcn is provided by the
  files profcndat (text version of the input file profcn.dat for the
  stand alone version), profort20 (text version of the output file
  fort.20 giving the resulting radial functions) and profort30 (text
  version of the output file fort.30 giving the resulting angular
  functions).

  Subroutine Version of profcn

    subroutine profcn(c,m,lnum,ioprad,x1,iopang,iopnorm,narg,arg, &
                      r1c,ir1e,r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr, &
                      s1c,is1e,s1dc,is1de,naccs)

        real(knd), intent (in)  ::  c, x1, arg(narg)
        integer, intent (in)    ::  m, lnum, ioprad, iopang, iopnorm, narg
        real(knd), intent (out) ::  r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                                    s1c(lnum, narg), s1dc(lnum, narg)
        integer, intent (out)   ::  ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                                    is1e(lnum, narg), is1de(lnum, narg), & 
                                    naccr(lnum), naccs(lnum, narg)

      Input and output parameters appearing in the subroutine call
      statement are defined below:

          c      : desired value of the size parameter (= kd/2, where
                   k = wavenumber and d = interfocal length) (either
                   real*8 or real*16)
          m      : desired value for the order m (integer)
          lnum   : number of values desired for the degree l equal
                   to m, m + 1, m + 2, ..., m + lnum - 1 (integer)
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed
          x1     : value of the radial coordinate x minus 1.0.
                   [real*8 or real*16]. This choice is made to avoid
                   subtraction errors in calculating quantities
                   containing x - 1 when x is close to unity.
                   (a nominal value can be entered for x1 if ioprad = 0).
                   x1 can range from very large to extremely small
                   values, including 0.0. Note that for x1 = 0.0,
                   ioprad must be set to 1 since the radial functions
                   of the second kind and their first derivatives are
                   infinite. 
          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed
          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      legendre function [i.e., we use the Meixner and
                      Schafke normalization scheme.] This norm
                      becomes very large as m becomes large. The
                      angular functions are computed below as
                      a characteristic and an exponent to avoid
                      overflow.
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.
                      This is very useful since it removes the
                      need to calculate a normalization factor
                      when using the angular function values given
                      here. It also eliminates any chance for
                      overflow when the characteristics and exponents
                      are combined to form the angular functions.
          narg   : number of values of the angular coordinate eta for
                   which angular functions are calculated (integer)
          arg:     vector containing the values of eta for which
                   angular functions are desired (real*8 or real*16)
          r1c   :  either real*8 or real*16 vectors of length lnum
          r1dc     containing the characteristics for the radial
                   radial functions of the first kind r1 and their
                   first derivatives
          ir1e   : integer vectors of length lnum containing the
          ir1de    exponents corresponding to r1c and r1dc
          r2c    : real*8 or real*16 vectors of length lnum containing
          r2dc     the characteristics for the radial functions of the
                   second kind r2 and their first derivatives
          ir2e   : integer vectors of length lnum containing the
          ir2de    exponents corresponding to r2c and r2dc
          naccr  : integer vector of length lnum containing the estimated
                   accuracy of the radial functions
          s1c,   : two-dimensional arrays s1c(lnum,narg) and
          s1dc     s1dc(lnum,narg) that contain narg calculated
                   characteristics for the angular functions and
                   their first derivatives for each of the lnum
                   values of l (real*8 or real*16)
                   For example, s1(10,1) is the characteristic
                   of the angular function for l = m +10 -1 and
                   for the first value of eta given by arg(1)
          is1e,    integer arrays is1e(lnum,narg) and is1de(lnum,narg)
          is1de    containing the exponents corresponding to s1c and
                   s1dc
          naccs  : two-dimensional array naccs(lnum,narg) that contains
                   narg estimated accuracy values for the angular functions
                   for each of the lnum values of l

  Stand alone Version of profcn

     Input Data

     Input parameters are read from unit 1 in the file profcn.dat
     assumed to be in the directory of profcn.f90. Profcn.dat
     contains the following lines of data:

       line 1:
          mmin   : minimum value for m. (integer)
          minc   : increment for m. (integer)
          mnum   : number of values of m. (integer)
          lnum   : number of values of l [l=m, l=m+1,
                   ..., l=m+lnum-1]. (integer)

       line 2:
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed

          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed

          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner-
                      Schafke normalization scheme.]
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.

       line 3:
          c      : value of the size parameter (= kd/2, where k =
                   wavenumber and d = interfocal length) (real(knd))
          x1     : value of the radial coordinate x minus one (real(knd))
                   (a nominal value of 10.0e0_knd can be entered for x1
                   if ioprad = 0). [When x1 = 0.0e0_knd, ie. x = 1.0e0_knd,
                   radial functions of the first kind and their first
                   derivatives are = 0.0e0_knd when m is not equal to 0.
                   Special formulas obtained as the limiting form of
                   equations given in the journal article cited above
                   for calculating radial functions of the first kind and
                   their first derivatives are used to obtain values
                   for m = 0. For all m, the radial functions of the second
                   kind and their first derivatives are infinite. Thus
                   ioprad must be set = 1 when x = 1.0e0_knd.]

       line 4:
          ioparg : (integer)
                 : =0 if both arg1 and darg are angles in degrees
                 : =1 if arg1 and darg are dimensionless values of eta

          arg1   : first value for the angle coordinate (in degrees
                   or dimensionless if eta) for which angular
                   functions are to be computed. (real(knd))

          darg   : increment used to calculate additional desired
                   arguments for angular functions. (real(knd))

          narg   : number of desired angle arguments. (integer)
                   (line 4 is not read when iopang = 0)

     Output files

     These output files are also available using the subroutine version
     of profcn. Generation of each of the files is controlled by a logical
     specified in the module param located before the program. False
     suppresses the output file and true enables it. The logical debug
     controls fort.30 and fort.40, the logical output controls fort.20
     and fort.30 and warn controls fort.60.

      fort.20

     This file contains values for all radial functions that have
     been calculated.
     The first line in the file contains the values for x, c, and
     m, formatted as follows (see statements 115 and 120 in subroutine
     main):

                x      : e23.14 in real*8; e40.30 in real*16
                c      : e23.14 in real*8; e40.30 in real*16
                m      : i5

     Each subsequent line in fort.20 contains radial functions
     for given values of l. The first line contains values for l = m,
     the next for l=m+1 and continuing to l=m+lnum-1. The radial
     functions are preceeded by the value of l and followed by the
     accuracy, equal to the estimated number of accurate decimal digits
     in the radial functions as measured either using the Wronskian
     or estimated based on subtraction errors and degree of convergence
     when the Wronskian is used to obtain r2 and r2d.

       The output and corresponding format for each line is as follows
       (see statements 690, 700, and 710 in subroutine main). [Only 15
       decimal digits are written to fort.20 using these formats. This
       is equivalent to the number of decimal digits in real*8
       arithmetic. If more digits are desired in real*16 arithmetic,
       the format statements can be modified. The exponents are formatted
       as I6, allowing for exponents from -99999 to +99999. This can be
       changed to accommodate a larger exponent range by increasing the
       format to I7 or larger.]

         for ioprad = 1 or 2:

               l      : value for l (i6)
               r1c    : characteristic of the prolate radial function
                        of first kind (f17.14)
               ir1e   : exponent of the prolate radial function of
                        first kind (i6)
               r1dc   : characteristic of the first derivative of the
                        prolate radial function of first kind (f17.14)
               ir1de  : exponent of the first derivative of the
                        prolate radial function of first kind (i6)

        for ioprad = 2, each line also includes;

               r2c    : characteristic of the prolate radial function
                        of second kind (f17.14)
               ir2e   : exponent of the prolate radial function of
                        second kind (i6). If the exponent for any
                        function is greater than 9999, the format
                        can be increased to i6 or higher. Note that
                        the procedures used in this program allow for
                        exponents much larger than those allowed
                        in real(knd) arithmetic on the users computer
                        since the floating point function values
                        provided are given as a characteristic and an
                        integer exponent. Use of ratios in calculating
                        the functions eliminates overflow during the
                        calculations.
               r2dc   : characteristic of the first derivative of the
                        prolate radial function of second kind (f17.14)
               ir2de  : exponent of the first derivative of the
                        prolate radial function of second kind (i6).
                        [See comment above for ir2e.]
               naccr  : accuracy: equal to the number of decimal digits
                        of agreement between the theoretical Wronskian
                        and the calculated Wronskian or estimated when
                        the Wronskian is used to obtain r2 and r2d
                        (i2). This is a measure of the accuracy of the
                        radial function of the second kind and its
                        first derivative. The algorithms for calculting
                        the radial functions of the first kind and
                        their first derivatives are both robust with
                        accurate results.

                        To indicate that the accuracy estimate naccr
                        has been obtained using the Wronskian, the
                        letter w follows naccr in the output. This
                        distinguishes it from (1) the case where the
                        theoretical Wronskian is used to obtain the
                        value of the denominator (series) appearing
                        in the variable eta expressions for the radial
                        function of the second kind and its
                        first derivative. Here naccr is instead
                        calculated using the degree of convergence and
                        associated subtraction errors of the two
                        numerators. The letter e follows naccr in the
                        output to designate this case. This version of
                        profcn may also use the Wronskian to compute the
                        leading coefficients when computing r2 and r2d
                        using the taditional Legendre function
                        expansion. This can sometimes provides accurate
                        values for x less than or equal to 1.01 when
                        all of the other methods fail and the
                        coefficients are inaccurate due to large
                        subtraction errors incurred in their
                        computation. This case is also indicated in the
                        output by the use of the letter e following the
                        accuracy value naccr.

     fort.30

     This file contains values for all angular functions
     that have been calculated. Its first line contains the values
     for c and m, formatted as follows (see statements 60 and 70 in
     subroutine main).

                c      : e23.14 (knd = 8) or e40.30 (knd = 16)
                m      : i5

     The second line in fort.30 contains the value for the first l (=m)
     formatted as follows (see statement 140 in subroutine main):

                l      : i6

     This is followed by a series of narg lines. Each line contains
     a desired value of angle (ioparg = 0) or angular coordinate eta
     (ioparg =1) followed by the corresponding angular functions and
     accuracy. Specific output and format for each line is as follows.
     [only 15 decimal digits are written to fort.30 using these
     formats. If more digits are desired in real*16 arithmetic, the
     formats can be modified.]

        for iopang = 1:

               arg    : for ioparg = 0, angle in degrees (f17.14; see
                        statement 750 in subroutine main)
            or barg   ; for ioparg = 1, angular coordinate eta
                        (f17.14; see statement 750 in subroutine main)
               s1c    : characteristic of the prolate angular function
                        of first kind (f17.14; see statement 750)
               is1e   : exponent of the prolate angular function of
                        first kind (i5; see statement 750)

        for iopang = 2, each line also includes:
               s1dc   : characteristic of the first derivative of the
                        prolate angular function of first kind (f17.14;
                        see statement 760 in subroutine main)
               is1de  : exponent of the first derivative of the
                        prolate angular function of first kind (i5;
                        see statement 760 in subroutine main)

        for iopang = 1 or 2:
               naccs  : accuracy: estimate of the number of decimal
                        digits of accuracy in both the angular function
                        (and its first derivative when iopang = 2). It
                        is a conservative estimate based on the
                        calculated subtraction error in the series
                        calculation of the angular function. When the
                        accuracy estimate is equal to 0, the
                        corresponding angular functions are set equal
                        to zero. (i2; see statements 750 and 760 in
                        subroutine main).

     fort.40 and fort.50

     These files are diagnostic files. They contain information about
     the specific techniques used and the numbers of terms required
     for the radial function and angular function calculations,
     respecively. They are annotated and should be self-explanatory.

     fort.60

     This file may be of interest to the user. It is recommended that
     the user utilize this file, especially when using profcn outside
     the ranges where the estimated accuracy is expected to be at least
     5 decimal digits. (see section 4 on accuracy below).
     Whenever the estimated accuracy for the radial functions falls below
     a designated integer value, the associated values of x, c, m, and l
     are written to fort.60. This integer is currently set equal to 6 in
     the write statement for fort.60 found just before the line numbered
     600 in subroutine main. The reason for choosing 6 is that it is expected
     that 6 accurate decimal digits are sufficient for most applications.
     The integer can be changed to any desired value in the write(60,*),
     If the eigenvalue routine in subroutine conver fails to converge to
     an eigenvalue between the two neigboring eigenvalues after 50 tries
     using progressively refined starting values, the m, l, and c values
     where this occurs will be written to fort.60. Note that this has never
     been observed over the many years that this routine has been used for
     prolate eigenvalues.

  4. Accuracy of Results

  Use of double precision arithmetic can provide accurate results up to
  values of c ranging from moderate to very large, depending on the value
  for x and m. Use of quadruple precision arithmetic provides good results
  over much wider parameter ranges. However, profcn runs significantly
  faster using double precision arithmetic. For both double and quadruple
  precision, the radial functions of the first kind r1 and their first
  derivatives r1d are nearly fully accurate unless extremely close to a root
  or when x is extremely close to unity where r1d can lose a few digits
  of accuracy when m = 0 and l is odd. Values of r1 and r1d for x = 1.0
  are nearly fully accurate.

  The radial functions of the second kind r2 and their first derivatives
  r2d are usually accurate to ten or more digits in double precision and
  25 or more digits in quadruple precision. But r2 and r2d can be much less
  accurate than this, especially when c is very large, x is between about
  1.001 and 1.1, and the values for m and l are intermediate.

  Profcn was tested extensively for both double and quadruple precision
  arithmetic. The compiler used provided a precision of 15 decimal digits in
  double precision (real*8) and 31 decimal digits in quadruple precision
  (real*16).

  A set of values for x was chosen that includes those regions where
  profcn has difficulty providing highly accurate values for r2 and
  r2d. Testing for real*8 arithmetic included all orders m up to 1000.
  Testing using real*16 arithmetic included values of m from 0 to 500
  in steps of 10 and 550 to 1000 in steps of 50. For each value of m,
  l went from m up to a sufficiently large value where r1 and r1d have
  magnitudes smaller than about 10 to the -300 power.

  The testing results are summarized in the table below. For each of
  the listed x values, the table gives the largest tested value of c
  for which all of the r2 and r2d values were found to have an
  estimated accuracy of at least 5 decimal digits, except in the case
  where either one is near a root of r2 or r2d or where the variable eta
  method is used to obtain r2 and r2d and where the Wronskian is used
  to obtain a sufficiently accurate value for the denominator term in
  this method. Examination of the test results showed that the function
  value near a root is smaller in magnitude than those values for
  neighboring degrees l by an order of magnitude no less than the
  difference between its estimated accuracy and the estimated accuracy
  of the neighboring function values. Since the neighboring values have
  an estimated accuracy of 5 or more digits, the effective accuracy of
  the function value near the root is at least 5 digits. This is
  because its contribution to the solution of problems involving these
  functions is reduced by an amount corresponding to its reduction in
  magnitude and accuracy.

  The estimated accuracy is usually found by comparing the theoretical
  value for the Wronskian to the value computed using the radial
  function values. Sometimes, the Wronskian is used to determine values
  for r2 and r2d and is unavailable to estimate the accuracy. Here a
  conservative estimate of accuracy is obtained using subtraction
  errors and the degree of convergence of the series involved in their
  calculation. The choice of 5 digits is based on the likelihood that
  5 digits of accuracy is sufficient for most applications of these
  functions. A rare 4 digit result should also be acceptable.

  Nearly all of the values obtained for r2 and r2d are much more
  accurate than 5 digits. Lower accuracy tends to occur at higher
  values of c and then only for a limited range of intermediate values
  of m and l. For a desired value of x, a conservative estimate of the
  limit on c is given by the smaller of the limits for the two values
  of x that it lies between.
  
  It is possible that a rare 4-digit result other than at a root will
  occur even when c is less than the appropriate table value. It is not
  expected that this will be a problem in use of these function values
  to solve problems involving them.
  
                        Approximate Upper limit for c

          x         real*8 arithmetic    real*16 arithmetic

     1.000005            5000
     1.00001             5000
     1.00005             3500                  5000
     1.0001              2500                  5000
     1.0005              1200                  3100
     1.001               1000                  2900
     1.002                800                  2100
     1.003                650                  1700
     1.004                600                  1600
     1.005                550                  1450
     1.006                500                  1350
     1.007                480                  1300
     1.008                450                  1200
     1.009                450                  1150
     1.01                 400                  1100
     1.02                 350                   900*
     1.03                 300                   800*
     1.04                 290                   800*
     1.05                 290                   900*
     1.06                 290                  1000*
     1.07                 290                  2000*
     1.08                 300                  5000
     1.09                 300
     1.10                4000**
     1.101               5000
     1.102               4000**
     1.103               2000**
     1.104               5000
     1.105               5000
     1.11                5000

   * For real*16 arithmetic and for x from about 1.02 to 1.07, much
     larger values of c than those give in the table above can provide
     at least 5 digits of accuracy for values of m up to a maximum
     value less than 1000 but possibly large enough for many
     applications. I summarize this by listing the largest value of m
     for specific values of x and c that provide 5 or more digits of
     accuracy up to the specified value of m. Following the value of x
     is a series of values of c, each followed in parenthesis by the
     maximum value of m for that value of c.

     x = 1.02: 1000(90); 1500(130); 2000(170); 3000(230); 4000(290);
               5000(350)
     x = 1.03: 1000(120); 1500(180); 2000(240); 3000(330); 4000(420);
               5000(490)
     x = 1.04: 1000(180); 1500(250); 2000(310); 3000(440); 4000(550);
               5000(650)
     x = 1.05: 1000(190); 1500(310); 2000(420); 3000(550); 4000(700);
               5000(850)
     x = 1.06: 1100(300); 1500(390); 2000(500); 3000(750); 4000(900);
               5000(1000)
     x = 1.07: 3000(950); 4000(1000); 5000(1000)

  ** There is a single 4 digit result that is not near a root for c
     = 5000 when x = 1.10 and when x = 1.102 and for c = 3000 and 5000
     when x = 1.103.

  Additional testing showed that profcn should provide highly accurate
  results for values of m above 1000 and for l well above those used
  for the table. Profcn should also provide accurate results for c up
  to at least 5000 when x is greater than 1.11 for real*8 arithmetic
  and when x is greater than or equal to 1.08 for real*16 arithmetic.
  It should provide accurate results for c up to at least 5000 for
  x from 1.00005 down to values at least as small as 1.000000000001
  when using either real*8 or real*16 arithmetic.

  Use of the Wronskian to estimate accuracy can sometimes overestimate
  the accuracy of either r2 or r2d. This happens when one of the two
  terms in the Wronskian, either r1*r2d or r2*r1d, is smaller in
  magnitude that the other term. Here its effect on the Wronskian
  is reduced according to the degree to which it is smaller. Thus the
  accuracy of either r2 or r2d, whichever is in the smaller term, can
  be less than the Wronksian test indicates. If the reason for the
  smaller term is that one of the radial functions is near a root, then
  its contribution to a calculation using these functions will be
  reduced accordingly. It is also possible to underestimate the
  accuracy if the two terms r1*r2d and r2*r1d are nearly equal. Profcn
  checks for this rare occurrence and adjusts the Wronskian accuracy
  estimate when it occurs.

  An integer called minacc is used in profcn to designate the
  minimum number of accurate digits desired for r2 and r2d. The value
  for minacc controls which methods are used to calculate r2 and r2d
  and in which order. Minacc is set equal to 10 for real*8 arithmetic.
  It is recommended that this not be changed. Minacc is set equal to
  15 for real*16 arithmetic. If more accuracy is desired, minacc can
  be increased. If greater speed is desired in difficult regions where
  real*8 arithmetic is insufficient, minacc can be reduced to a value
  as low as 8 digits. The value of minacc is set in a statement early
  in the program.

  The calculated angular functions are highly accurate except for lower
  values of l (less than about 2c/pi) when c is large. They tend to be
  less accurate the larger the value of c, the smaller the value of
  l - m and the closer eta is to unity (i.e., the closer theta is to 0
  degrees). However, the loss of accuracy (in decimal digits) is due to
  subtraction error and is accompanied by a proportional decrease in
  the magnitude of the angular function relative to its corresponding
  associated Legendre function. This decrease in magnitude almost
  always results in a corresponding reduction in the magnitude of their
  contribution to the solution of physical problems involving prolate
  spheroidal functions. Thus the lower accuracy in some of the angular
  functions almost always has insignificant impact on calculated
  solutions.

  5. Obtaining the d expansion coefficients

  The user may desire values for the d coefficients that appear in
  the expression for the angular functions as well as in many of the
  expressions used to calculate the radial functions. Ratios of
  successive d coefficients are stored in the vector enr where enr(k)
  = d(subscript 2k+ix) divided by d(subscript 2k-2+ix). The vector enr
  is calculated in the subroutine dnorm in statement 20 and passed to
  subroutine main. The number lim2 of d coefficients calculated for a
  given l is chosen to be sufficient to compute radial and angular
  functions for that l. The size of lim2 necessary to compute r1, r1d
  and s1, s1d and the normaliation factors ranges for low l from less
  than 100 for low c and somewhat less than int(c) for large c. Lim2
  increases with increasing l, eventually becoming less than l in size.
  The size of lim2 needed to compute r2 and r2d can be comparable to
  this unless they are computed using Neumann function expansions. Then
  lim2 can be much larger, especially for x very close to unity. Note
  that the vector enr returned by the subroutine conver contains scaled
  ratios. The scaling factors are removed in subroutine dnorm to obtain
  the desired d coefficient ratios.

  The d coefficients themselves can be obtained starting with the value
  for d with the subscript l - m. If iopnorm is set = 0, Oblfcn uses
  the Meixner-Schafke scheme for normalizing the angular functions.
  Here they have the same norm as the corresponding associated Legendre
  functions. Computation of this normalization is very accurate since
  the series involved has only positive terms. The subroutine dnorm
  computes d(subscript l-m) for this normalization and returns it as
  a characteristic dmlms and an exponent idmlmse to subroutine main.
  Use of an exponent avoids possible overflow of d(subscript l-m) for
  extremely large c and m. When the user sets iopnorm = 1 so that the
  angular functions have unit norm, the corresponding characteristic
  and exponent for d(subscript l-m) are calculated in subroutine s1
  and returned to subroutine main as dmlms1 and idmlms1e. Values for
  the characteristics and exponents of d(subscript l-m) for the  Morse-
  Feshbach and Flammer normalizations are computed in dnorm and r1bes,
  respectively, and returned to main as dmlmf, idmlmfe and dmlf,
  idmlfe. Calculation of the Morse and Feshbach normalization suffers
  subtraction errors for lower values of l-m and large c that increase
  as c increases. The value for dmlmf will have reduced accuracy in
  this case.

  6. Obtaining the eigenvalues

  The eigenvalues for the prolate functions are computed in subroutine
  conver and returned to main where they are stored in the vector eig(l+1).
  There is such a vector created for each value of m.