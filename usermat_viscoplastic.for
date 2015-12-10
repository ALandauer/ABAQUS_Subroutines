      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
    !
      INCLUDE 'ABA_PARAM.INC'
    !     WARNING - the aba_param.inc file declares
    !        Implicit real*8(a-h,o-z)
    !     This means that, by default, any variables with
    !     first letter between a-h or o-z are double precision.
    !     The rest are integers.
    !     Note that this also means that if you type a variable
    !     name incorrectly, the compiler won't catch your typo.
    !
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


    !
    !      DDSDDE(NTENS,NTENS)
    !         Jacobian matrix of the constitutive model.
    !         DDSDDE(I,J) defines the change in the Ith stress component
    !         at the end of the time increment caused by an infinitesimal
    !         perturbation of the Jth component of the strain increment array.
    !         Unless you invoke the unsymmetric equation solution capability
    !         for the user-defined material, ABAQUS/Standard will use only
    !         the symmetric part of DDSDDE. The symmetric part of the matrix
    !         is calculated by taking one half the sum of the matrix and its transpose.


    !      STRESS(NTENS)
    !         This array is passed in as the stress tensor at the beginning
    !         of the increment and must be updated in this routine to be the
    !         stress tensor at the end of the increment. If you specified
    !         initial stresses (â€œInitial conditions,â€? Section 19.2.1), this
    !         array will contain the initial stresses at the start of the
    !         analysis. The size of this array depends on the value of NTENS
    !         as defined below. In finite-strain problems the stress tensor
    !         has already been rotated to account for rigid body motion in
    !         the increment before UMAT is called, so that only the corotational
    !         part of the stress integration should be done in UMAT. The
    !         measure of stress used is â€œtrueâ€? (Cauchy) stress.
    !
    !      STATEV(NSTATV)
    !         An array containing the solution-dependent state variables.
    !         These are passed in as the values at the beginning of the
    !         increment unless they are updated in user subroutines USDFLD
    !        (â€œUSDFLD,â€? Section 25.2.39) or UEXPAN (â€œUEXPAN,â€? Section 25.2.20),
    !        in which case the updated values are passed in. In all cases
    !         STATEV must be returned as the values at the end of the increment.
    !         The size of the array is defined as described in
    !        â€œAllocating spaceâ€? in â€œUser subroutines: overview,â€? Section 25.1.1.
    !
    !         In finite-strain problems any vector-valued or tensor-valued
    !         state variables must be rotated to account for rigid body
    !         motion of the material, in addition to any update in the
    !         values associated with constitutive behavior. The rotation
    !         increment matrix, DROT, is provided for this purpose.
    !
    !      SSE, SPD, SCD
    !         Specific elastic strain energy, plastic dissipation, and
    !         â€œcreepâ€? dissipation, respectively. These are passed in as
    !         the values at the start of the increment and should be
    !         updated to the corresponding specific energy values at
    !         the end of the increment. They have no effect on the solution,
    !         except that they are used for energy output.
    !
    !     Only in a fully coupled thermal-stress analysis
    !      RPL
    !         Volumetric heat generation per unit time at the end of the increment
    !         caused by mechanical working of the material.
    !
    !     DDSDDT(NTENS)
    !          Variation of the stress increments with respect to the temperature.
    !
    !     DRPLDE(NTENS)
    !           Variation of RPL with respect to the strain increments.
    !
    !     DRPLDT
    !           Variation of RPL with respect to the temperature.
    !
    !     Variables that can be updated
    !
    !     PNEWDT
    !        Ratio of suggested new time increment to the time increment being
    !        used (DTIME, see discussion later in this section). This variable
    !        allows you to provide input to the automatic time incrementation
    !        algorithms in ABAQUS/Standard (if automatic time incrementation is chosen).
    !        For a quasi-static procedure the automatic time stepping that ABAQUS/Standard
    !        uses, which is based on techniques for integrating standard creep laws
    !        (see â€œQuasi-static analysis,â€? Section 6.2.5), cannot be controlled from within
    !        the UMAT subroutine.
    !        PNEWDT is set to a large value before each call to UMAT.
    !        If PNEWDT is redefined to be less than 1.0, ABAQUS/Standard must abandon the
    !        time increment and attempt it again with a smaller time increment. The
    !        suggested new time increment provided to the automatic time integration
    !        algorithms is PNEWDT Ã— DTIME, where the PNEWDT used is the minimum value
    !        for all calls to user subroutines that allow redefinition of PNEWDT for this
    !        iteration.
    !        If PNEWDT is given a value that is greater than 1.0 for all calls to user
    !        subroutines for this iteration and the increment converges in this iteration,
    !        ABAQUS/Standard may increase the time increment. The suggested new time increment
    !        provided to the automatic time integration algorithms is PNEWDT Ã— DTIME, where
    !        the PNEWDT used is the minimum value for all calls to user subroutines for
    !        this iteration.
    !        If automatic time incrementation is not selected in the analysis procedure,
    !        values of PNEWDT that are greater than 1.0 will be ignored and values of
    !        PNEWDT that are less than 1.0 will cause the job to terminate.
    !
    !    Variables passed in for information
    !
    !     STRAN(NTENS)
    !         An array containing the total strains at the beginning of the increment.
    !         If thermal expansion is included in the same material definition, the
    !         strains passed into UMAT are the mechanical strains only (that is, the
    !         thermal strains computed based upon the thermal expansion coefficient have
    !         been subtracted from the total strains). These strains are available for output
    !         as the â€œelasticâ€? strains.
    !
    !         In finite-strain problems the strain components have been rotated to account for
    !         rigid body motion in the increment before UMAT is called and are approximations
    !         to logarithmic strain.

    !     DSTRAN(NTENS)
    !         Array of strain increments. If thermal expansion is included in the same
    !         material definition, these are the mechanical strain increments (the total
    !         strain increments minus the thermal strain increments).
    !
    !     TIME(1)
    !         Value of step time at the beginning of the current increment.
    !
    !     TIME(2)
    !          Value of total time at the beginning of the current increment.
    !
    !     DTIME
    !        Time increment.
    !
    !     TEMP
    !         Temperature at the start of the increment.
    !
    !     DTEMP
    !         Increment of temperature.
    !
    !     PREDEF
    !        Array of interpolated values of predefined field variables at this point
    !        at the start of the increment, based on the values read in at the nodes.
    !
    !      DPRED
    !        Array of increments of predefined field variables.
    !
    !      CMNAME
    !        User-defined material name, left justified. Some internal material models are given
    !		 names starting with the â€œABQ_â€? character string. To avoid conflict, you should not use â€œABQ_â€? as the leading string for CMNAME.
    !
    !      NDI
    !        Number of direct stress components at this point.
    !
    !      NSHR
    !        Number of engineering shear stress components at this point.
    !
    !      NTENS
    !        Size of the stress or strain component array (NDI + NSHR).
    !
    !      NSTATV
    !         Number of solution-dependent state variables that are associated with
    !         this material type (defined as described in â€œAllocating spaceâ€? in â€œUser
    !         subroutines: overview,â€? Section 25.1.1).
    !
    !      PROPS(NPROPS)
    !         User-specified array of material constants associated with this user material.
    !
    !      NPROPS
    !         User-defined number of material constants associated with this user material.
    !
    !      COORDS
    !         An array containing the coordinates of this point. These are the current
    !         coordinates if geometric nonlinearity is accounted for during the step
    !         (see â€œProcedures: overview,â€? Section 6.1.1); otherwise, the array contains
    !         the original coordinates of the point.
    !
    !     DROT(3,3)
    !          Rotation increment matrix. This matrix represents the increment of rigid
    !          body rotation of the basis system in which the components of stress
    !          (STRESS) and strain (STRAN) are stored. It is provided so that vector- or
    !          tensor-valued state variables can be rotated appropriately in this subroutine:
    !          stress and strain components are already rotated by this amount before UMAT
    !          is called. This matrix is passed in as a unit matrix for small-displacement
    !          analysis and for large-displacement analysis if the basis system for the
    !          material point rotates with the material (as in a shell element or when a
    !          local orientation is used).
    !
    !      CELENT
    !          Characteristic element length, which is a typical length of a line across
    !          an element for a first-order element; it is half of the same typical length
    !          for a second-order element. For beams and trusses it is a characteristic length
    !          along the element axis. For membranes and shells it is a characteristic length
    !          in the reference surface. For axisymmetric elements it is a characteristic length
    !          in the  plane only. For cohesive elements it is equal to the constitutive
    !          thickness.
    !
    !      DFGRD0(3,3)
    !          Array containing the deformation gradient at the beginning of the increment.
    !          See the discussion regarding the availability of the deformation gradient for
    !          various element types.
    !
    !     DFGRD1(3,3)
    !            Array containing the deformation gradient at the end of the increment.
    !           The components of this array are set to zero if nonlinear geometric effects
    !           are not included in the step definition associated with this increment. See
    !           the discussion regarding the availability of the deformation gradient for
    !           various element types.
    !
    !      NOEL
    !           Element number.
    !
    !      NPT
    !           Integration point number.
    !
    !      LAYER
    !          Layer number (for composite shells and layered solids).
    !
    !      KSPT
    !          Section point number within the current layer.
    !
    !      KSTEP
    !         Step number.
    !
    !     KINC
    !         Increment number.

    !      user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
    !      and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
    !
    !     Local variables
    
      double precision Emod,xnu,Y,eps_0,xn,epsdot_0,xm,tol,sequiv,error,
     1e,f,c,dfde,eplas,gamma0,del_ij,del_ik,del_jl,del_jk,del_il,del_kl,
     2enew,beat,factor,count0
      dimension Ss_i(NTENS),Sn_i(NTENS),de_i(NTENS),Ss_ij(NDI,NDI),
     1Sn_ij(NDI,NDI),de_ij(NDI,NDI),Cmat(NDI,NDI,NDI,NDI),
     2stress_ij(NDI,NDI)
      integer i,j,k,l
	
      ddsdde = 0.d0
      Cmat = 0.d0
      sequiv = 0.d0
      error = 0.d0
      c = 0.d0
      f = 0.d0
      count0 = 0.d0

      Emod = props(1)
      xnu = props(2)
      Y = props(3)
      eps0 = props(4)
      xn = props(5)
      epsdot0 = props(6)
      xm = props(7)

      eplas = statev(1)

      !Compute deviatoric strain increment and stress
      do i = 1,NTENS
        if (i .lt. NDI+1) then
            de_i(i) = dstran(i) - sum(dstran(1:NDI))/3.d0
            Sn_i(i) = stress(i) - sum(stress(1:NDI))/3.d0
        else
            de_i(i) = dstran(i)
            Sn_i(i) = stress(i)
        end if
		
      end do
C      write(6,*) Sn_i
	  !Compute the matrix representations, since the element locs are hard-coded
	  !use if-then structure
      if (NTENS == 6) then
      de_ij(1,1) = de_i(1)
      de_ij(2,2) = de_i(2)
      de_ij(3,3) = de_i(3)
      de_ij(1,2) = de_i(4)
      de_ij(1,3) = de_i(5)
      de_ij(2,3) = de_i(6)
      de_ij(2,1) = de_i(4)
      de_ij(3,1) = de_i(5)
      de_ij(3,2) = de_i(6)

      Sn_ij(1,1) = Sn_i(1)
      Sn_ij(2,2) = Sn_i(2)
      Sn_ij(3,3) = Sn_i(3)
      Sn_ij(1,2) = Sn_i(4)
      Sn_ij(1,3) = Sn_i(5)
      Sn_ij(2,3) = Sn_i(6)
      Sn_ij(2,1) = Sn_i(4)
      Sn_ij(3,1) = Sn_i(5)
      Sn_ij(3,2) = Sn_i(6)
	  
      else if (NTENS == 4) then
      de_ij(1,1) = de_i(1)
      de_ij(2,2) = de_i(2)
      de_ij(1,2) = de_i(3)
      de_ij(2,1) = de_i(3)
	  
      Sn_ij(1,1) = Sn_i(1)
      Sn_ij(2,2) = Sn_i(2)
      Sn_ij(1,2) = Sn_i(3)
      Sn_ij(2,1) = Sn_i(3)
      end if
	  
      !Compute the elastic predictors
      do i = 1,NDI
        do j = 1,NDI
            Ss_ij(i,j) = Sn_ij(i,j) + Emod/(1.d0+xnu)*de_ij(i,j)
            sequiv = sequiv + Ss_ij(i,j)*Ss_ij(i,j)
        end do
      end do
      sequiv = sqrt(1.5*sequiv)
      
      !N-R loop
      e = 10.d0**(-15.d0)
      error = Y
      tol = 10.d0**(-6.d0)*Y
      if (sequiv*sum(stran(1:NTENS)*2) == 0) then
        e = 0.d0
      else
        do while (error .gt. tol .and. count0 < 250)
      c=(1.d0+(eplas+e)/eps0)**(1.d0/xn)*(e/(dtime*epsdot0))**(1.d0/xm)
      f = sequiv/Y - 1.5*e*Emod/(Y*(1.d0+xnu))- c
      dfde=-1.5*Emod/(Y*(1.d0+xnu))-c*(1.d0/(xn*(eplas+e+e0))
     1 +1.d0/(xm*e))
      enew = e - f/dfde
      count0 = count0 + 1
      if (enew<0.d0) then
                !        e must be >0, so if new approx to e <0 the solution
                !        must lie between current approx to e and zero.
                e = e/10.d0
            else
                e = enew
            end if
            err = abs(f)
        end do
      end if
	     
      if (sequiv*e == 0.d0) then
       stress(1:NTENS) = 0.d0
	     
      else
C       Calculate stress (update values)
       do i = 1,NTENS
        if (i .lt. NDI+1) then
         stress(i)=(1-3.d0*Emod*e/(2.d0*(1.d0+xnu)*sequiv))*Ss_i(i)+
     1   (sum(stress(1:NDI))+Emod*sum(dstran(1:NDI))/(1.d0-2.d0*xnu))
     2   /3.d0
        else
         stress(i)=(1.d0-3.d0*Emod*e/(2.d0*(1.d0+xnu)*sequiv))*Ss_i(i)
        end if
      
       end do 
      end if
	  
      if (NTENS == 6) then
       stress_ij(1,1) = stress(1)
       stress_ij(2,2) = stress(2)
       stress_ij(3,3) = stress(3)
       stress_ij(1,2) = stress(4)
       stress_ij(1,3) = stress(5)
       stress_ij(2,3) = stress(6)
       stress_ij(2,1) = stress(4)
       stress_ij(3,1) = stress(5)
       stress_ij(3,2) = stress(6)
      else if (NTENS == 4) then
       stress_ij(1,1) = stress(1)
       stress_ij(2,2) = stress(2)
       stress_ij(1,2) = stress(3)
       stress_ij(2,1) = stress(3)
      end if	 
C      write(6,*) stress_ij
	  
      if (sequiv*e == 0) then
       beta = 1.d0
       factor = 0.d0
      else 
       beta = 1.d0/(1.d0+1.5*Emod*e/((1.d0+xnu)*sequiv))
       gamma0=beta*(1.5*Emod/((1.d0+xnu)*sequiv)+
     1  (1.d0/(xn*(eps0+eplas+e))+1.d0/(xm*e)))
       factor = 1.5*1.5*Emod*(e-1.d0/gamma0)/((1.d0+xnu)*sequiv**3.d0)
      end if
	 
    !find the tangent stiffness matrix
	 !Calculate C_ijkl
      do i = 1,NDI
        do j = 1,NDI
            do k = 1,NDI
                do l = 1,NDI

            del_ij = 0.d0
            del_ik = 0.d0
            del_jl = 0.d0
            del_jk = 0.d0
            del_il = 0.d0
            del_kl = 0.d0

            if (i==j) then del_ij = 1.d0
            if (i==k) then del_ik = 1.d0
            if (j==l) then del_jl = 1.d0
            if (j==k) then del_jk = 1.d0
            if (i==l) then del_il = 1.d0
            if (k==l) then del_kl = 1.d0

      Cmat(i,j,k,l)=beta*Emod/(1.d0+xnu)*((del_ik*del_jl+del_kj*del_il)
     1  /2.d0-del_ij*del_kl/3.d0+factor*stress_ij(i,j)*stress_ij(k,l))
     2  +Emod/(3.d0*(1.d0-2.d0*xnu))*del_ij*del_kel
     
C     1 Emod/(1.d0+xnu)*(1.d0-3.d0*Emod*e/(2.d0*
C     1 (1.d0+xnu)*sequiv))*(0.5*(del_ik*del_jl+del_jk*del_il)
C     2 -1.d0/3.d0*(del_ij*del_kl))+Emod/(1.d0+xnu)*(9.d0*Emod*
C     3 (e-1.d0/gamma0))/(4.d0*(1.d0+xnu)*sequiv)*(Ss_ij(i,j)/sequiv)
C     4 *Ss_ij(k,l)/sequiv+Emod*del_ij*del_kl/(3.d0*(1.d0-2.d0*xnu))

                end do
            end do
        end do
      end do
	  
      
	  
	  !Fill in the DDEDDE matrixs from the elements of the complete C_ijkl matrix
	  !Either for the 3D or 2D case, since the elements are hard-coded use if-then
      if (NTENS == 6) then
      DDSDDE(1,1) = Cmat(1,1,1,1)
      DDSDDE(1,2) = Cmat(1,1,2,2)
      DDSDDE(1,3) = Cmat(1,1,3,3)
      DDSDDE(1,4) = Cmat(1,1,1,2)
      DDSDDE(1,5) = Cmat(1,1,3,1)
      DDSDDE(1,6) = Cmat(1,1,2,3)
      DDSDDE(2,1) = Cmat(2,2,1,1)
      DDSDDE(2,2) = Cmat(2,2,2,2)
      DDSDDE(2,3) = Cmat(2,2,3,3)
      DDSDDE(2,4) = Cmat(2,2,1,2)
      DDSDDE(2,5) = Cmat(2,2,3,1)
      DDSDDE(2,6) = Cmat(2,2,2,3)
      DDSDDE(3,1) = Cmat(3,3,1,1)
      DDSDDE(3,2) = Cmat(3,3,2,2)
      DDSDDE(3,3) = Cmat(3,3,3,3)
      DDSDDE(3,4) = Cmat(3,3,1,2)
      DDSDDE(3,5) = Cmat(3,3,3,1)
      DDSDDE(3,6) = Cmat(3,3,2,3)
      DDSDDE(4,1) = Cmat(1,2,1,1)
      DDSDDE(4,2) = Cmat(1,2,2,2)
      DDSDDE(4,3) = Cmat(1,2,3,3)
      DDSDDE(4,4) = Cmat(1,2,1,2)
      DDSDDE(4,5) = Cmat(1,2,3,1)
      DDSDDE(4,6) = Cmat(1,2,2,3)
      DDSDDE(5,1) = Cmat(3,1,1,1)
      DDSDDE(5,2) = Cmat(3,1,2,2)
      DDSDDE(5,3) = Cmat(3,1,3,3)
      DDSDDE(5,4) = Cmat(3,1,1,2)
      DDSDDE(5,5) = Cmat(3,1,3,1)
      DDSDDE(5,6) = Cmat(3,1,2,3)
      DDSDDE(6,1) = Cmat(2,3,1,1)
      DDSDDE(6,2) = Cmat(2,3,2,2)
      DDSDDE(6,3) = Cmat(2,3,3,3)
      DDSDDE(6,4) = Cmat(2,3,1,2)
      DDSDDE(6,5) = Cmat(2,3,3,1) 
      DDSDDE(6,6) = Cmat(2,3,2,3)
      else if (NTENS == 4) then
      DDSDDE(1,1) = Cmat(1,1,1,1)
      DDSDDE(1,2) = Cmat(1,1,2,2)
      DDSDDE(1,3) = Cmat(1,2,2,1)
      DDSDDE(2,1) = Cmat(2,2,1,1)
      DDSDDE(2,2) = Cmat(2,2,2,2)
      DDSDDE(2,3) = Cmat(2,2,2,1)
      DDSDDE(3,1) = Cmat(2,1,1,1)
      DDSDDE(3,2) = Cmat(2,1,2,2)
      DDSDDE(3,3) = Cmat(2,1,2,1)
      end if
	  

      !update state vars
      if (e*sequiv == 0.d0) then
       statev(1) = 0.d0
      else
       statev(1) = eplas+e
      end if
    !    for debugging, you can use
    !      write(6,*) ' Hello '
    !    Output is then written to the .dat file
	

    !
    !     NOTE: ABAQUS uses engineering shear strains,
    !     i.e. stran(ndi+1) = 2*e_12, etc...

        RETURN
      END
	 
