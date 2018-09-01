from sympy import diff
from sympy import lambdify
from sympy.core.add import Add
from sympy.abc import x
from sympy.abc import t
from numpy import array as np_array
from numpy import linalg as la
import timeit

from builtins import isinstance
 
def wrapper(func, a_list):
    def wrapped():
        if (len(a_list)==5):
            return func(a_list[0],a_list[1],a_list[2],a_list[3],a_list[4])
        else:
            return func()
    return wrapped


def newton_ralphson(expr,initial_guess,strategy,tolerance = 1e-5, max_iteration = 100):
    if not (isinstance(initial_guess,list)):
        print("initial guess mush be a list")
        return []
    r = lambdify((x,t),expr)
    exp_der_u = diff(expr,x)
    exp_der_l = diff(expr,t)
    u_k = initial_guess[0] ## k is 0
    lambda_k = initial_guess[1] ## k is 0
    drdu = lambdify(x,exp_der_u)
    drdl = lambdify(t,exp_der_l)
    i = 0
    print(r)
    if(strategy == "displacement control" or strategy == "d"):
        assert(len(initial_guess)==2 and isinstance(expr, Add))
        while(tolerance < abs(r(u_k,lambda_k))):
            if(i == max_iteration):
                break
            lambda_k = lambda_k + (-r(u_k,lambda_k)/drdl(lambda_k))
            i+=1
        print("r error is : ",r(u_k,lambda_k))
    elif(strategy == "arc length control" or strategy == "a"):
        ##[u_guess,l_guess,u_eq,l_eq,arc_leng ]
        assert(len(initial_guess)==5 and isinstance(expr, Add))
        arc_length = initial_guess[4]
        u_eq = initial_guess[2]
        l_eq = initial_guess[3]
        c_exp = (x - u_eq)**2 + (t - l_eq)**2 - (arc_length*arc_length)
        c_der_u_exp = diff(c_exp,x)
        c_der_l_exp = diff(c_exp,t)
        c = lambdify((x,t),c_exp)
        dcdu = lambdify(x,c_der_u_exp)
        dcdl = lambdify(t,c_der_l_exp)
        while(tolerance < abs(r(u_k,lambda_k))):
            if(i == max_iteration):
                break
            lhs_1_1 = float(drdu(u_k))
            lhs_1_2 = float(drdl(lambda_k))
            lhs_2_1 = float(dcdu(u_k))
            lhs_2_2 = float(dcdl(lambda_k))
            rhs_0 = float(-r(u_k,lambda_k))
            
            lhs = np_array([[lhs_1_1,lhs_1_2],[lhs_2_1,lhs_2_2]])
            
            rhs = np_array([[rhs_0],[0]])
            
            sol = la.solve(lhs, rhs)
            
            u_k += sol[0]
            lambda_k += sol[1]
            
            i+=1
        print("r error is : ",r(u_k,lambda_k),"\t, c error is : ",c(u_k,lambda_k))
    elif(strategy =="load control" or strategy == "l"):
        assert(len(initial_guess)==2 and isinstance(expr, Add))
        while(tolerance < abs(r(u_k,lambda_k))):
            if(i == max_iteration):
                break
            u_k += (-r(u_k,lambda_k)/drdu(u_k))
            i+=1
        print("r error is : ",r(u_k,lambda_k))
    else:
        print("you entered wrong parameters, strategy a is for arc-length and d is for disp f.e. ")
        return []
         
         
         
         
         
    return [float(u_k),float(lambda_k)]
                 

# u = Symbol('u')
# l = Symbol('l')

exp = (2/125) * (16*x - (6*x*x) + 0.5*x*x*x) + 0.5*x*x + 1.3*t*t
# print(type(exp))
# print(isinstance(exp, Add))

c_exp = (x - -0.15667)**2 + (t - 0.15247)**2 - (0.1*0.1)
c = lambdify((x,t),c_exp)
r = lambdify((x,t),exp)
# exp_der_u = diff(exp,x)
# exp_der_l = diff(exp,t)
# drdu = lambdify(x,exp_der_u)
# drdl = lambdify(t,exp_der_l)
# print(r(-0.25450538,0.17346463),c(-0.25450538,0.17346463))


# print(drdu(-0.156667))
# print(drdl(0.23))
# print(newton_ralphson(exp, [-0.255121,0.17,-0.15667,0.15247,0.1],"a",1e-8,40))

short_list1= [exp, [-0.15667,0.15247],"d",1e-8,40]
short_list2= [exp, [-0.05667,0.12],"l",1e-8,40]

wrapped1 = wrapper(newton_ralphson,short_list1)
wrapped2 = wrapper(newton_ralphson,short_list2)


print(timeit.timeit(wrapped1,number =1))
print(timeit.timeit(wrapped2,number =1))

# print(newton_ralphson(exp, [-0.15667,0.15247],"d",1e-8,40))
# 
# print(newton_ralphson(exp, [-0.05667,0.12],"l",1e-8,40))



