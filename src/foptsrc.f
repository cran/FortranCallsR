
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Library General Public
C License as published by the Free Software Foundation; either
C version 2 of the License, or (at your option) any later version.
C
C This library is distributed in the hope that it will be useful, 
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
C GNU Library General Public License for more details.
C
C You should have received a copy of the GNU Library General 
C Public License along with this library; if not, write to the 
C Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
C MA 02111-1307 USA


C This is a toy model for a box bounded one-dimensioanl optimizer which re-
C quires only a function call to the objective function fun(n, x, value).
C The optimzer itself is written in Fortran, the objective function will be
C passed from the R environment via a C wrapper.

C So this example can be considered as a starting point for similar tasks.
C These are tasks where we have Fortran Code which optimizes, finds zeroes,
C integrates or even is good for other things, and the objective function
C should be made availalble through the R environment.

C The toy model considered here is Brent's one dimensional optimzer (1973). 
C We have slightly modified the original code which we copied from GAMS
C 'golden oldies' section. The optimizer and objective function are Fortran 
C subroutines. 

C Our goal would be later to implement other Fortran Solvers (which have a 
C similar structure) into the R environment. This can be done in a straight-
C forward way since the toy model can eaily be extended to higher dimensions,
C since the gradient and/or Hessian can be included in the same way as the
C objective function, and evenmore the box bounds can easily be extended to
C handle more complex constraints. This will be done in a second step for
C the 'sequential quadratic programming' optimizer written in Fortran by 
C Ladislav Luksan.

C Diethelm Wuertz
C ETH Zurich
C www.rmetrics.org
  

c ------------------------------------------------------------------------------   
c PART I: Objective Function


      subroutine func(n, x, value)
      implicit double precision (a-h, o-z)
      double precision x(n), value(1)    
      call cfunc(n, x, value)
      return
      end      

      
c ------------------------------------------------------------------------------   
c PART II: Optimizer 


      subroutine fopt(par, nf, xlo, xup, value, mf, rpar, np)
      implicit double precision (a-h, o-z)
      double precision par(1), xlo(1), xup(1), value(1), rpar(1)
      
      x = par(1)
      tol = rpar(1)
      ax = xlo(1)
      bx = xup(1)
    
c     Input | Output:  
c     x(1) | output - an approximation to the minimum of f(x)
c     xlo(1), xup(1) | Input - lower and upper endpoint of interval
c     value(1) | Output - the value of f at the mimimum
c     rpar(1) | Input - tolerance parameter
c     nf, mf, np | Not used, here all integers are one
c
c ------------------------------------------------------------------------------     
c Copy from GAMS' golden oldies (with some minor modifications):
c License unknown ?
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result (.ge.0.)
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a
c        minimum
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
c  then fmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
c  c is the squared inverse of the golden ratio
      c=0.5d0*(3.0d0-dsqrt(5.0d0))
c
   10 continue
c  eps is approximately the square root of the relative machine precision.
c  To get d1mach, mail netlib send d1mach from core
c  DW: not required here for the toy model we use: eps = 1.0d-12
c  eps=d1mach(4) replaced by
      eps = 1.0d-8
      tol1=eps+1.0d0
      eps=dsqrt(eps)
c
      a=ax
      b=bx
      v=a+c*(b-a)
      w=v
      x=v
      e=0.0d0
      call func(1, x, fx)
      fv=fx
      fw=fx
      tol3=tol/3.0d0
c
c  main loop starts here
c
   20 xm=0.5d0*(a+b)
      tol1=eps*dabs(x)+tol3
      t2=2.0d0*tol1
c
c  check stopping criterion
c
      if (dabs(x-xm).le.(t2-0.5d0*(b-a))) go to 190
      p=0.0d0
      q=0.0d0
      r=0.0d0
      if (dabs(e).le.tol1) go to 50
c
c  fit parabola
c
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0d0*(q-r)
      if (q.le.0.0d0) go to 30
      p=-p
      go to 40
   30 q=-q
   40 r=e
      e=d
   50 if ((dabs(p).ge.dabs(0.5d0*q*r)).or.(p.le.q*(a-x))
     2          .or.(p.ge.q*(b-x))) go to 60
c
c  a parabolic-interpolation step
c
      d=p/q
      u=x+d
c
c  f must not be evaluated too close to ax or bx
c
      if (((u-a).ge.t2).and.((b-u).ge.t2)) go to 90
      d=tol1
      if (x.ge.xm) d=-d
      go to 90
c
c  a golden-section step
c
   60 if (x.ge.xm) go to 70
      e=b-x
      go to 80
   70 e=a-x
   80 d=c*e
c
c  f must not be evaluated too close to x
c
   90 if (dabs(d).lt.tol1) go to 100
      u=x+d
      go to 120
  100 if (d.le.0.0d0) go to 110
      u=x+tol1
      go to 120
  110 u=x-tol1
  120 continue
      call func(1, u, fu)
c
c  update  a, b, v, w, and x
c
      if (fx.gt.fu) go to 140
      if (u.ge.x) go to 130
      a=u
      go to 140
  130 b=u
  140 if (fu.gt.fx) go to 170
      if (u.ge.x) go to 150
      b=x
      go to 160
  150 a=x
  160 v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
      go to 20
  170 if ((fu.gt.fw).and.(w.ne.x)) go to 180
      v=w
      fv=fw
      w=u
      fw=fu
      go to 20
  180 if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) go to 20
      v=u
      fv=fu
      go to 20
c
c  end of main loop
c
  190 continue
      
      par(1) = x
      call func(1, x, fx)
      value(1) = fx
      
      return
      end
c ------------------------------------------------------------------------------
