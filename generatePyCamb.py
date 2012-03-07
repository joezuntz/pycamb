#!/usr/bin/env python

# Most of the camb parameters are members of a derived type called CambParams, 
# the main instance of which is called CP.  Parameters added here are by default
# assumed to be in that type.

# If the parameter name starts with @ that signifies it is *not* part of the camb params type but
# is defined elsewhere in CAMB. For example, the dark energy parameters w_lam and cs2_lam are defined
# in the LambdaGeneral module in equations.F90.  To be modified by this code the parameters
# need to be accessible in camb.f90

# There are various members of CambParams which are themselves derived types, such as 
# Init__Power and Reion.  To access these variables, use a double underscore to 
# indicate what would be a % symbol in fortran.
# e.g. reion__redshift

# Parameters can be given default values to take if they are not set.
# CAMB already provides defaults for may values in the subroutine
# CAMB_SetDefParams in camb.f90; if these are not suitable then you can set them here.
# Gloabl parameters (starting in @) must have default values.

# Parameters with nasty names, like InitPower__ScalarPowerAmp(1), can be given friendly names
# in the "alias" dictionary.

# These are the  numerical parameters that we will let the user pass to camb.
numericalParams=['omegab', 'omegac', 'omegav', 'omegan','H0','TCMB',
                'yhe','Num_Nu_massless','Num_Nu_massive','omegak',
                'reion__redshift','reion__optical_depth','reion__fraction',"reion__delta_redshift",
                'Scalar_initial_condition','scalar_index','scalar_amp',
                'scalar_running','tensor_index','tensor_ratio',
                '@AccuracyBoost','@lAccuracyBoost','@lSampleBoost','@w_lam','@cs2_lam',
]

# The boolean (logical) options that the user can pass to camb.
logicalParameters=['WantScalars', 'WantTensors','reion__reionization','reion__use_optical_depth','@w_perturb']

#Default values, a dictionary
defaultValues={'@AccuracyBoost':1.0,'@lAccuracyBoost':1.0,'@lSampleBoost':1.0,'scalar_amp':2.1e-9,'@w_lam':-1.0,'@w_perturb':False,'@cs2_lam':1.0}

#JZ These map the friendly names in the parameters above to the initial power vectors
alias={
'scalar_index': 'InitPower__an(1)',
'scalar_amp': 'InitPower__ScalarPowerAmp(1)',
'scalar_running':'InitPower__n_run(1)',
'tensor_index':'InitPower__ant(1)',
'tensor_ratio': 'InitPower__rat(1)'
}




#If you just want to add new parameters, don't change anything below here.





nparams=len(numericalParams)+len(logicalParameters)

endl="\n"
from os import system


#Generate the fortran subroutine that is called directly from python.
def makeFortranClsCaller(number_parameters):
    code= """
    subroutine getcls(paramVec,lmax,cls)
        use camb
        implicit none
        real, intent(in) :: paramVec(%d)
        integer, intent(in) :: lmax
        real, intent(out) ::  cls(2:lmax,4)
        type(CAMBparams) :: P
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        P%%max_l=lmax
        P%%Max_l_tensor=lmax
        P%%Max_eta_k=2*lmax
        P%%Max_eta_k_tensor=2*lmax
        call CAMB_GetResults(P)
        call CAMB_GetCls(cls, lmax, 1, .false.)
        cls=cls*output_cl_scale
    end subroutine getcls
    
    """ % number_parameters
    
    return code


def makeFortranAgeCaller(number_parameters):
    code= """
    subroutine getage(paramVec,age)
        use camb
        implicit none
        real, intent(in) :: paramVec(%d)
        double precision, intent(out) ::  age
        type(CAMBparams) :: P
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        age = CAMB_GetAge(P)
    end subroutine getage

    """ % number_parameters

    return code


def makeFortranTransfersCaller(number_parameters):
    code = """subroutine gentransfers(paramVec,lmax,nred,redshifts)
        use camb
        implicit none
        real, intent(in) :: paramVec(%d)
        integer, intent(in) :: nred, lmax
        double precision, intent(in), dimension(nred) :: redshifts
        type(CAMBparams) :: P
        integer :: nr, i
        nr = size(redshifts)
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        P%%WantTransfer = .true.
        P%%max_l=lmax
        P%%Max_l_tensor=lmax
        P%%Max_eta_k=2*lmax
        P%%Max_eta_k_tensor=2*lmax
        P%%transfer%%num_redshifts = nr
        do i=1,nr
            P%%transfer%%redshifts(i)=redshifts(i)
        enddo
        call CAMB_GetResults(P)
        allocate(transfers(Transfer_max,MT%%num_q_trans,nred))
        allocate(transfers_k(MT%%num_q_trans))
        allocate(transfers_sigma8(nred))
        transfers = MT%%TransferData
        transfers_k = MT%%q_trans
        transfers_sigma8 = MT%%sigma_8(:,1)
    end subroutine gentransfers

    subroutine freetransfers()
        deallocate(transfers)
        deallocate(transfers_k)
        deallocate(transfers_sigma8)
    end subroutine freetransfers

""" % number_parameters
    return code

def makeFortranHeader():
    return """module pycamb_mod
    double precision, dimension(:,:,:), allocatable :: transfers
    double precision, dimension(:), allocatable :: transfers_k,transfers_sigma8
    real, dimension(:,:,:), allocatable :: matter_power
    double precision, dimension(:,:), allocatable :: matter_power_kh
    double precision :: output_cl_scale=7.4311e12
    contains
    """



def makeFortranExtraFunctions(number_parameters):
    code="""
    function angularDiameter(paramVec,z)
        use ModelParams, only : CAMBparams, camb_angulardiameter => AngularDiameterDistance
        use camb, only : CAMB_SetDefParams, CAMBParams_Set
        implicit none
        double precision, intent(in) :: z
        double precision :: angularDiameter
        real, intent(in) :: paramVec(%d)
        integer error
        type(CAMBparams) :: P        
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        call CAMBParams_Set(P,error)
        angularDiameter = camb_angulardiameter(z)
    end function angularDiameter

    subroutine angularDiameterVector(paramVec,n,z,ang)
        !Should be in temporal order ie redshift decreasing
        use ModelParams, only : CAMBparams, DeltaTime, rofchi
        use camb, only : CAMB_SetDefParams, CAMBParams_Set, CP
        implicit none
        integer, intent(in) :: n
        integer :: nz
        double precision,dimension(n), intent(in) :: z
        double precision,dimension(n), intent(out) :: ang
        real, intent(in) :: paramVec(%d)
        integer error
        integer i,j
        type(CAMBparams) :: P        
        nz=size(z)
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        call CAMBParams_Set(P,error)
        ang(nz) = rofchi(DeltaTime(1/(1+z(nz)),1.0_8)/CP%%r)
        do i=1,nz-1
            j=nz-i
            ang(j) = rofchi(DeltaTime(1.0/(1.0+z(j)),1./(1.+z(j+1)))/CP%%r) + ang(j+1)
        enddo
        ang = ang * CP%%r/(1+z)
    end subroutine angularDiameterVector

    """ % (number_parameters,number_parameters)
    return code
    
    
def makeFortranMatterPowerCaller(number_parameters):
    code = """
    subroutine getpower(paramVec,maxk,dlogk,nred,redshifts)
        use camb
        implicit none
        real, intent(in) :: paramVec(%d)
        integer, intent(in) :: nred
        integer :: lmax
        double precision, intent(in), dimension(nred) :: redshifts
        type(CAMBparams) :: P
        integer :: nr, i
        real, intent(in) ::  maxk, dlogk
        real, parameter :: minkh = 1.0e-4
        integer in,itf, points, points_check
        nr = size(redshifts)
        call CAMB_SetDefParams(P)
        call makeParameters(paramVec,P)
        P%%WantTransfer = .true.
        lmax=10000
        P%%max_l=lmax
        P%%Max_l_tensor=lmax
        P%%Max_eta_k=2*lmax
        P%%Max_eta_k_tensor=2*lmax
        P%%transfer%%num_redshifts = nr
        
        P%%transfer%%kmax = maxk * (P%%h0/100._dl)
        P%%transfer%%k_per_logint = dlogk
        
        do i=1,nr
            P%%transfer%%redshifts(i)=redshifts(i)
        enddo
        call CAMB_GetResults(P)

        itf=1
         points = log(MT%%TransferData(Transfer_kh,MT%%num_q_trans,itf)/minkh)/dlogk+1
        allocate(matter_power(points,CP%%InitPower%%nn,nr))
        allocate(matter_power_kh(points,nr))
        
        do itf=1, CP%%Transfer%%num_redshifts
            points_check = log(MT%%TransferData(Transfer_kh,MT%%num_q_trans,itf)/minkh)/dlogk+1
             if (points_check .ne. points)  stop 'Problem with pycamb assumption on k with z'
             do in = 1, CP%%InitPower%%nn
              call Transfer_GetMatterPower(MT,matter_power(:,in,itf), itf, in, minkh,dlogk, points)
!              if (CP%%OutputNormalization == outCOBE) then
!                 if (allocated(COBE_scales)) then
!                  outpower(:,in) = outpower(:,in)*COBE_scales(in)
!                 else
!                  if (FeedbackLevel>0) write (*,*) 'Cannot COBE normalize - no Cls generated'
!                 end if
!             end if
             end do
             do i=1,points
              matter_power_kh(i,itf)=minkh*exp((i-1)*dlogk)
             end do
        enddo !End redshifts loop
    end subroutine getpower

    subroutine freepower()
        deallocate(matter_power)
        deallocate(matter_power_kh)
    end subroutine freepower
"""% number_parameters
    return code

#Generate the fortran subroutine called by the first one to set the CAMBparams instance's values from the vector we pass it. 
def makeFortranParamSetter(numericalParameters,logicalParameters,defaultValues):
    code="""
    subroutine makeParameters(paramVec,P)
        use camb
        implicit none
        real, intent(in) :: paramVec(%d)
        type(CAMBparams) P
    """ % (len(numericalParameters)+len(logicalParameters))
    
    nn=len(numericalParameters)
    nm=len(numericalParameters) + len(logicalParameters)
    
    for p,friendly_name in enumerate(numericalParameters):
        if alias.has_key(friendly_name):
            name=alias[friendly_name]
        else:
            name=friendly_name
        fortranName="P%%%s" % name.replace('__','%')
        if name.startswith("@"):
            fortranName="%s" % (name.lstrip("@").replace('__','%') )
            if name not in defaultValues:
                raise ValueError("Global parameter %s requires default argument (or it is retained between runs).  Edit the top of %s" % (name,__file__))

        if friendly_name in defaultValues:
            print "%s has default" % name
            code +="""
if (paramVec(%d) .ne. -1.6375e30) then
    %s = paramVec(%d)
else
    %s = %e
endif
""" % (p+1,fortranName,p+1,fortranName,defaultValues[friendly_name])
        else:
            print "%s has NO default" % name     
            code += "if (paramVec(%d) .ne. -1.6375e30) %s = paramVec(%d)\n" % (p+1,fortranName,p+1)



    for p,friendly_name in enumerate(logicalParameters):
        if alias.has_key(friendly_name):
            name=alias[friendly_name]
        else:
            name=friendly_name
        fortranName="P%%%s" % name.replace('__','%')
        if name.startswith("@"):
            fortranName="%s" % (name.replace('__','%').lstrip("@") )
            if name not in defaultValues:
                raise ValueError("Global parameter %s require default argument (or it is retained between runs).  Edit the top of %s" % (name,__file__))
        if friendly_name in defaultValues:
            if defaultValues[friendly_name]:
                defaultValue = '.true.'
            else:
                defaultValue = '.false.'    
            code += """ 
if (paramVec(%d) .ne. -1.6375e30) then 
    if(paramVec(%d) .ne. 0.0) then 
        %s=.true. 
    else 
        %s=.false. 
    endif
else
    %s=%s
endif
"""  % (nn+p+1,nn+p+1,fortranName,fortranName,fortranName,defaultValue) + endl
        else:  #No default value
            code += """
if (paramVec(%d) .ne. -1.6375e30) then 
    if(paramVec(%d) .ne. 0.0) then 
        %s=.true. 
    else 
        %s=.false. 
    endif
endif
""" % (nn+p+1,nn+p+1,fortranName,fortranName)
        
    code += "    end subroutine makeParameters" + endl +endl
    
    return code

#Generate the python code that will wrap the fortran.
def makePython(numericalParameters,logicalParameters,defaultValues):

    param_string=(", ".join([name.lstrip("@") for name in numericalParameters+logicalParameters])) + "\n"
    alias_string="\n    ".join(["%s ---> %s" % (name.strip("@"),alias[name]) for name in alias])
    defparam_string = "\n".join("%s (%e)" % (name.lstrip("@"),param) for (name,param) in defaultValues.items())
    header="""
__all__=["camb","age","angular_diameter","matter_power","transfers"]
from numpy import zeros,repeat, array, float64
import sys
import _pycamb
_getcls = _pycamb.pycamb_mod.getcls
_getage = _pycamb.pycamb_mod.getage
_gentransfers = _pycamb.pycamb_mod.gentransfers
_freetransfers = _pycamb.pycamb_mod.freetransfers
_getpower = _pycamb.pycamb_mod.getpower
_freepower = _pycamb.pycamb_mod.freepower
_angulardiameter = _pycamb.pycamb_mod.angulardiameter
_angulardiametervector = _pycamb.pycamb_mod.angulardiametervector

numericalParameters=%s
logicalParameters=%s
defaultValues=%s
np=len(numericalParameters)+len(logicalParameters)
nn=len(numericalParameters)
nm=len(numericalParameters)+len(logicalParameters)
""" % (str(numericalParameters),str(logicalParameters),str(defaultValues))

    pvec_builder ="""
def build_pvec(**parameters):
    pvec=repeat(-1.6375e30,np)
    input_params=parameters.keys()
    for n,p in enumerate(numericalParameters):
        pin=p
        if p.startswith('@'): pin=p.lstrip('@')
        if pin in input_params:
            pvec[n]=parameters[pin]
            input_params.remove(pin)
    for n,p in enumerate(logicalParameters):
        pin=p
        if p.startswith('@'): pin=p.lstrip('@')
        if pin in input_params:        
            if parameters[pin]:
                pvec[n+nn]=1.0
            else:
                pvec[n+nn]=0.0
            input_params.remove(pin)

    if input_params:
        print "WARNING: Unrecognized parameters:"
        for p in input_params:
            print p
    return pvec
"""

    camb_code = """    

def camb(lmax,**parameters):
    \"""
    Run camb up to the given lmax, with the given parameters and return the Cls.  Parameter names are case-insensitive.
    
    Any parameters that are not specified in the input are left with camb's default values,as given by the CAMB_SetDefParams subroutine in camb.f90.
    
    You can either specify parameters as keywords:
        cl = camb(1000,H0=72.0,omegab=0.04)
        
    Or using a dictionary:
        camb_params={"H0":72.0,"omegab":0.04}
        cl=camb(1000,**camb_params)
    
    The latter method gives you more flexibilty when writing codes.
    
    Valid parameters are all specified in the script that generated this code.  Most are elements of the CAMBparams derived type specified in modules.f90 .
    Parameters that are not in CAMBparams can be set by pre-pending an underscore.
    
    Parameters that are members of a sub-typewithin the CAMBparams type, (that is, all those accessed in fortran using CP%%(something)%%(varname)  ) such as the reionization parameters, should be given with the percent symbols replaced with a double underscore, __.
    
    For example, to specify the redshift of reionization, given by CP%%reion%%redshift in camb use:
        cl=camb(1000,reion__redshift=11.0)
    
    Boolean (logical) parameters can be specified as you would expect:
        cl=camb(1000,reionization=False)
        
    You can let more parameters be passed into camb by modifiying the top of, and running, generatePyCamb.py

    In this code, valid normal parameters are:
    %s
    
    And parameters with default values are:
%s
    
    
    Parameters which take a vector (like the primordial power spectrum parameters) are not yet properly implemented here, 
    so only one spectral index, amplitude, etc., at a time is possible.  Keyword parameters that are mapped to camb names are:
        %s
    
    \"""
    pvec=build_pvec(**parameters)
    cls=_getcls(pvec,lmax)
    return cls.transpose()
    
    """ % (param_string,defparam_string,alias_string)

    age_code = """

def age(**parameters):
    \"""
Get the age of the unverse with the given parameters.  See the docstring for pycamb.camb for more info on parameters.
    \"""
    pvec=build_pvec(**parameters)
    age=_getage(pvec)
    return age
    
"""
    trans_code = """
def transfers(redshifts=[0],**parameters):
    \"""
Generate transfer functions for the cosmology.  Specify the redshifts in DECREASING order (ie in temporal order).
The cosmological parameters that can be varied in the code are the same as for the 'camb' function - see help on that function for details
    \"""
    lmax=1000
    pvec=build_pvec(**parameters)
    redshifts = array(redshifts,dtype=float64)
    ordered_redshifts = redshifts.copy()
    ordered_redshifts.sort()
    ordered_redshifts=ordered_redshifts[::-1]
    if not (redshifts == ordered_redshifts).all(): sys.stderr.write("WARNING:  Re-ordered redshift vector to be in temporal order.  Ouput will be similarly re-ordered.\\n")
    if len(redshifts)>500: raise ValueError("At most 500 redshifts can be computed without changing the hardcoded camb value")
    
    _gentransfers(pvec,lmax,ordered_redshifts)
    T = _pycamb.pycamb_mod.transfers.copy()
    K = _pycamb.pycamb_mod.transfers_k.copy()
    S = _pycamb.pycamb_mod.transfers_sigma8.copy()
    _freetransfers()
    return K,T,S

"""
    power_code = """
def matter_power(redshifts=[0],maxk=1.,logk_spacing=0.02,**parameters):
    \"""
    Generate matter power spectra for the cosmology.  Specify the redshifts in DECREASING order (ie in temporal order).
    You can also set the maximum k value and the spacing in log(k)
    The cosmological parameters that can be varied in the code are the same as for the 'camb' function - see help on that function for details
    \"""
    pvec=build_pvec(**parameters)
    redshifts = array(redshifts,dtype=float64)
    ordered_redshifts = redshifts.copy()
    ordered_redshifts.sort()
    ordered_redshifts=ordered_redshifts[::-1]
    if not (redshifts == ordered_redshifts).all(): sys.stderr.write("WARNING:  Re-ordered redshift vector to be in temporal order.  Ouput will be similarly re-ordered.\\n")
    if len(redshifts)>500: raise ValueError("At most 500 redshifts can be computed without changing the hardcoded camb value")
    _getpower(pvec,maxk,logk_spacing,ordered_redshifts)
    power=_pycamb.pycamb_mod.matter_power.copy()
    kh=_pycamb.pycamb_mod.matter_power_kh.copy()
    _freepower()
    return kh.squeeze(),power.squeeze()
"""
    extra_code = """
def angular_diameter(z,**parameters):
    \"""
    Compute the angular diameter distance to the given redshift for the specified cosmology.
    You can also set the maximum k value and the spacing in log(k)
    The cosmological parameters that can be varied in the code are the same as for the 'camb' function - see help on that function for details
        \"""    
    pvec=build_pvec(**parameters)
    if isinstance(z,float) or len(z)==1:
        return _angulardiameter(pvec,z)
    redshifts = array(z,dtype=float64)
    ordered_redshifts = redshifts.copy()
    ordered_redshifts.sort()
    ordered_redshifts=ordered_redshifts[::-1]
    if not (redshifts == ordered_redshifts).all(): sys.stderr.write("WARNING:  Re-ordered redshift vector to be in temporal order.  Ouput will be similarly re-ordered.\\n")
    return _angulardiametervector(pvec,ordered_redshifts)
        
    
    """


    return header + pvec_builder + camb_code + age_code + trans_code + power_code + extra_code
    

def makeFortranFooter():
    return "end module pycamb_mod"
    
def main():
    fortran_file=open("py_camb_wrap.f90","w")
    fortran_file.write(makeFortranHeader())
    number_parameters = len(numericalParams) + len(logicalParameters)
    fortran_file.write(makeFortranClsCaller(number_parameters))
    fortran_file.write(makeFortranAgeCaller(number_parameters))
    fortran_file.write(makeFortranTransfersCaller(number_parameters))
    fortran_file.write(makeFortranMatterPowerCaller(number_parameters))
    fortran_file.write(makeFortranExtraFunctions(number_parameters))
    fortran_file.write(makeFortranParamSetter(numericalParams,logicalParameters,defaultValues))
    fortran_file.write(makeFortranFooter())
#JZ This is where you would add calls to other routines that make more fortran code, for example to wrap other camb functions.
    fortran_file.close()
    
    python_file=open("pycamb.py","w")
    python_file.write(makePython(numericalParams,logicalParameters,defaultValues))
    python_file.close()
    compile_command="f2py -c -m _pycamb -L. -lcamb -lgomp py_camb_wrap.f90 skip: makeparameters :"
    print compile_command
    system(compile_command)

    
if __name__=="__main__":
    main()
