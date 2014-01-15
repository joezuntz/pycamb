#!/usr/bin/env python
import numpy.f2py
import re
import warnings

# All the numerical parameters that we will let the user pass to camb.
# If the parameter name starts with @ that signifies it is *not* part of the camb params type but
# is defined elsewhere. For example, the dark energy parameters w_lam and cs2_lam are defined in the
# LambdaGeneral module in equations.F90
numericalParams=['omegab', 'omegac', 'omegav', 'omegan','H0','TCMB',
                'yhe','Num_Nu_massless','Num_Nu_massive','omegak',
                'reion__redshift','reion__optical_depth','reion__fraction',"reion__delta_redshift",
                'Scalar_initial_condition','scalar_index','scalar_amp',
                'scalar_running','tensor_index','tensor_ratio','nonlinear',
                '@lAccuracyBoost','@lSampleBoost','@w_lam','@cs2_lam',
                '@AccuracyBoost'             , 
]



defaultValues={'@AccuracyBoost':1.0,
'@lAccuracyBoost':1.0,
'@lSampleBoost':1.0,
'scalar_amp':2.1e-9,
'@w_lam':-1.0,
'@w_perturb':False,
'@cs2_lam':1.0,
'nonlinear':0,
}


# The boolean options that the user can pass to camb.
logicalParameters=['WantScalars', 'WantTensors','reion__reionization','reion__use_optical_depth','@w_perturb','DoLensing']


# These map the friendly names in the parameters above to the initial power vectors
alias={
'scalar_index': 'InitPower__an(1)',
'scalar_amp': 'InitPower__ScalarPowerAmp(1)',
'scalar_running':'InitPower__n_run(1)',
'tensor_index':'InitPower__ant(1)',
'tensor_ratio': 'InitPower__rat(1)'
}



nparams=len(numericalParams)+len(logicalParameters)

endl="\n"
from os import system


def template_file(infile, outfile, parameters):
    text = open(infile).read()
    for param,value in parameters.items():
        pattern = r"\$%s\$"%param
        count=0
        for i in xrange(1000):
            text,c = re.subn(pattern, str(value), text)
            count += c
            if c==0: break
        else:
            raise ValueError("Got into a long loop on template_file.  Did your template values contain $$ symbols?")
        if count==0:
            warnings.warn("Unused template parameter: %s" % param)
    if "$" in text: raise ValueError("Template contains unfilled parameters")
    open(outfile,"w").write(text)




#Generate the fortran subroutine called by the first one to set the CAMBparams instance's values from the vector we pass it. 
def makeFortranParamSetter(numericalParameters,logicalParameters,defaultValues):
    code=""
    
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
            # print "%s has default" % name
            code +="""
if (paramVec(%d) .ne. -1.6375e30) then
    %s = paramVec(%d)
else
    %s = %e
endif
""" % (p+1,fortranName,p+1,fortranName,defaultValues[friendly_name])
        else:
            # print "%s has NO default" % name     
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
        
    
    return code



import pprint
#Generate the python code that will wrap the fortran.
def makePython(numericalParameters,logicalParameters,defaultValues):
    param_string=(", ".join([name.lstrip("@") for name in numericalParameters+logicalParameters])) + "\n"
    alias_string="\n\t".join(["%s ---> %s" % (name.strip("@"),alias[name]) for name in alias])
    defparam_string = "\n".join("%s (%f)" % (name.lstrip("@"),param) for (name,param) in defaultValues.items())

    params = {"numericalParameters":numericalParameters, 
              "logicalParameters":logicalParameters, 
              "defaultValues":defaultValues,
              "param_string":param_string,
              "alias_string":alias_string,
              "defparam_string":defparam_string,
              }

    template_file("templates/pycamb.py.template", "src/__init__.py", params)

def makeFortran(numericalParams,logicalParameters,defaultValues):
    number_parameters = len(numericalParams) + len(logicalParameters)
    param_caller_function = makeFortranParamSetter(numericalParams,logicalParameters,defaultValues)
    params = {"number_parameters":number_parameters, "param_caller_function":param_caller_function}
    template_file("templates/py_camb_wrap.f90.template", "src/py_camb_wrap.f90",params)


    
def main():
    makeFortran(numericalParams, logicalParameters, defaultValues)
    makePython(numericalParams, logicalParameters, defaultValues)

if __name__=="__main__":
    main()
