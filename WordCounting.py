textfile = open("C:/Artyom/NEWTXT/Chakrabarty/Shakespeare_Byron.txt", 'r') #Uncomment this to read directly from file
######### dim_q should be even! Calculations take about dim_q**2/1000 min !!!
dim_q =40
text = textfile.read() #Uncomment this to store all text in text variable
#text = ("  ") # Put text in .txt file and read from there
# Cleaning text and lower casing all words
for char in '.,:;«»“”!?—"': # symbols to be deleted
    text=text.replace(char,' ')
text = text.lower()
text=text.replace(" - ", " ")
#print(text)
#print(len(text))


from collections import Counter


copy_text = list(text.split())
#print(copy_text)

import math
n = len(copy_text) 
n_by_2 = math.floor(n/2)

#print(n, n_by_2)

# Not very important this but helps sometimes
# from nltk import FreqDist
# word_dist = FreqDist(copy_text)
# print(dict(word_dist))


R_list = [0]
diff_list = []
for i in range(len(copy_text)): 
    if copy_text[i] not in diff_list: # If word has not occurred before 
        diff_list.append(copy_text[i]) # gets appended to different word list
        R_list.append(len(diff_list)) # Number of different words is updated
        #print(i+1, len(R_list))
    elif copy_text[i] in diff_list: # If word has occurred before
        R_list.append(len(diff_list))# number of different words stays the same and gets appended to list 
        #print(i+1, len(R_list))
#print(R_list)


R_n = R_list[-1]
#print(R_n)
R_n_by_2 = R_list[n_by_2+1]
if n_by_2 == n/2:
    R_n_by_2= (R_list[n_by_2 + 1]+ R_list[n_by_2])/2
    
#print(R_n_by_2)
r_exp = R_n / R_n_by_2
theta_hat = math.log(r_exp, 2)
#print(theta_hat)
theta_hat=min(theta_hat,0.95)
#print(theta_hat)
theta=theta_hat
alpha=1/theta

import numpy as np
from numpy import sqrt, sin, cos, pi, exp, log

def myzeta(alpha,q):
    x=0
    for j in range(1,101):
        x=x+exp(log(j+q-1)*(-alpha))
        
    x=x+exp(log(100.5+q-1)*(1-alpha))/(alpha-1)
    return  x




from scipy.special import factorial, binom, gamma, loggamma, zeta, gammainc
q=19.55
q_left=-0.9
q_right=40
for i in range(1,21):
    pp=[0]
    r_list=[0]
    C=1/zeta(alpha,q+1)
    ss=0
    M=1000
    N=M+0.5+q
    for j in range(1,M+1):
        p=C*(j+q)**(-alpha)
        pp.append(p)
    ER=0
    for j in range(1,M+1):
        ER=ER+1-(1-pp[j])**n
    ER=ER+(n*C)**theta*gammainc(1-theta,n*C*N**(-alpha))*gamma(1-theta)-N*(1-exp(-n*C*N**(-alpha)))
 #   print(R_list[-1])
  #  print(ER)
   # print(q)
    #print()
    
    if ER<R_list[-1]:
        q_left=q
        q=(q_left+q_right)/2
        
    if ER>R_list[-1]:
        q_right=q
        q=(q_left+q_right)/2

for k in range(1,n+1):
    ER=0
    for j in range(1,M+1):
        ER=ER+1-(1-pp[j])**k
    ER=ER+(k*C)**theta*gammainc(1-theta,k*C*N**(-alpha))*gamma(1-theta)-N*(1-exp(-k*C*N**(-alpha)))
    r_list.append(ER)
    

Z_list=[0]
for i in range(1,n+1):
    zz=(R_list[i]-r_list[i])/(R_list[n]**0.5)
    Z_list.append(zz)


    

omega_2=0
for i in range(1,n+1):
    print(copy_text[i-1],R_list[i],r_list[i],Z_list[i])
for i in range(1,n):
    omega_2=omega_2+Z_list[i]*(2*Z_list[i]+Z_list[i+1])/3/n

print()    
#print('n, R_n, theta_hat, q_hat, omega^2')    
print('n  ',n)
print('R_n ',R_n)
print('theta_hat ',theta_hat)
print('q_hat ',q)    
print('omega^2 ',omega_2)

#theta_hat=1

theta=theta_hat

def kernel(tau,t):
    return (tau+t)**theta_hat-(max(tau,t))**theta_hat
#print(kernel(0.5,0.5)) 

def kernel_0(s,t):
    return kernel(s,t)-s**theta*kernel(1,t)-t**theta*kernel(1,s)+s**theta*t**theta*kernel(1,1)

#print(kernel_0(0.5,0.5))



def kernel_hat(s,t):
    x=kernel_0(s,t)- t**theta*log(t)*(kernel(s,1)-2**theta*kernel(s, 1/2))/log (2)
    x=x-s**theta*log(s)*(kernel(t,1)-2**theta*kernel(t, 1/2))/log (2)
    x=x+s**theta*t**theta*(log(s)+log(t))*(kernel(1,1)-2**theta*kernel(1, 1/2))/log (2)
    x=x+s**theta*t**theta*log(s)*log(t)*(kernel(1,1)-2**(theta+1)*kernel(1, 1/2)+2**(2*theta)*kernel(1/2, 1/2))/log (2)**2
    return  x

#print(kernel_hat(0.5,0.5))



import scipy.integrate as integrate
import scipy.special as special
#result = integrate.quad(lambda x: kernel_0(1,x), 0, 1)
#print(result)
#print(result[0])



q_ij = [[0] * dim_q for i in range(dim_q)]

#from scipy.integrate import dblquad

print('matrix')


import mpmath
from mpmath import *
mp.dps = 25; mp.pretty = True

def A_i(theta,i):
    return pi*i*hyp1f2(theta/2+1, 3/2, theta/2+2, -i**2*pi**2/4)/(theta+2)

def B_i(theta,i):
    return pi*i*2**(theta+2)*hyp1f2(theta/2+1, 3/2, theta/2+2, -i**2*pi**2)/(theta+2) - A_i(theta,i)

def C_i(theta,i):
    return hyp1f2(theta/2+1, 1/2, theta/2+2, -i**2*pi**2/4)/(theta+2)

def D_i(theta,i):
    return 2**(theta+2)*hyp1f2(theta/2+1/2, 1/2, theta/2+3/2, -i**2*pi**2)/(theta+1)-2**(theta+2)*hyp1f2(theta/2+1, 1/2, theta/2+2, -i**2*pi**2)/(theta+2)-2*C_i(theta-1,i)+C_i(theta,i)

def E_ij(theta,i,j):
    return A_i(theta,i+j)/2+A_i(theta,i-j)/2

def F_i(theta,i):
    return (-1)**i* B_i(theta,i) + ((-1)**i-1)/pi/i

def G_i(theta,i):
    return - pi* i* hyp2f3(theta/2+1, theta/2+1, 3/2, theta/2+2, theta/2+2, -i**2*pi**2/4)/(theta+2)**2

def H_i(theta,i):
    x= pi* i *cos (pi* i /2) * (3/2)**(theta+2)* hyp1f2(theta/2+1, 3/2, theta/2+2, -9*i**2 *pi**2/16 )/(theta+2)
    x= x - pi* i *cos (pi* i /2) * 2**(-theta-2)* hyp1f2(theta/2+1, 3/2, theta/2+2, -i**2 *pi**2/16)/(theta+2)
    x=x-sin(pi*i/2)*(3/2)**(theta+1)* hyp1f2(theta/2+1/2, 1/2, theta/2+3/2, -9*i**2 *pi**2/16 )/(theta+1)
    x=x+sin(pi*i/2)*2**(-theta-1)* hyp1f2(theta/2+1/2, 1/2, theta/2+3/2, -i**2*pi**2/16)/(theta+1)
    x=x+(cos (pi* i /2)-1)*2**(-theta)/(pi*i) -  A_i(theta,i)
    x=x+pi*i*2**(-theta-2)* hyp1f2(theta/2+1, 3/2, theta/2+2, -i**2 *pi**2/16 )/(theta+2)
    return x

def J_ij(theta,i,j):
    x=0
    if i!=j:
        x=(i*A_i(theta,j)-j*A_i(theta,i)-(-1)**(i+j)*(i*B_i(theta,j)-j*B_i(theta,i)))/pi/(i**2-j**2)
    else:
        x=(A_i(theta,i)-B_i(theta,i))/2/pi/i-(C_i(theta,i)+D_i(theta,i))/2
    return x       
        

#for i in range(10):
#    print(J_ij(0,1,1))


q_ij_new = [[0] * (dim_q+2) for i in range(dim_q+2)]



for i in range(1,dim_q+1):
    print(i, end=' ')
    for j in range(1,i+1):
        x=J_ij(theta,i,j)-  (A_i(theta,i)-E_ij(theta,i,j))/pi/j  - (A_i(theta,j)-E_ij(theta,j,i))/pi/i
        x=x- A_i(theta,i)*F_i(theta,j)-A_i(theta,j)*F_i(theta,i) + kernel(1,1)*A_i(theta,i)*A_i(theta,j)
        x=x-G_i(theta,j)*(F_i(theta,i)-2**theta*H_i(theta,i))/log(2)
        x=x-G_i(theta,i)*(F_i(theta,j)-2**theta*H_i(theta,j))/log(2)
        x=x+(A_i(theta,i)*G_i(theta,j)+A_i(theta,j)*G_i(theta,i))*(kernel(1,1)-2**theta*kernel(1,1/2))/log(2)
        x=x+G_i(theta,i)*G_i(theta,j)*(kernel(1,1)-2**(theta+1)*kernel(1,1/2)+2**(2*theta)*kernel(1/2,1/2))/log(2)**2
        q_ij_new[i][j]=x
        q_ij_new[j][i]=x
        
#        print(q_ij_new[i][j], end=' ')
print()

q_ij_true = [[0.] * dim_q for i in range(dim_q)]


import math

for i in range(dim_q):
    for j in range(dim_q):
        q_ij_true[i][j]=math.floor(q_ij_new[i+1][j+1]*(10)**(18))*(10)**(-18)
#        print(q_ij_true[i][j], end=' ')
#    print()


#eighenvalues    

import numpy as np


from numpy import sqrt, sin, cos, pi



from numpy import linalg as LA


eiv, v = LA.eigh(q_ij_true)


#print('eigenvalues')
#print(eiv)


lambda_list=[]
for i in range(dim_q):
    ll=1/eiv[dim_q-i-1]
    lambda_list.append(ll)

#print('lambda_list')

#print(lambda_list)

#print(omega_2)


#Smirnov's formula


def Determ(lam,lambda_list):
    kk=math.floor(len(lambda_list)/2)
    D=1
    for k in range(1,2*kk+1):
        D=D*(1-lam/lambda_list[k-1])
    return D    


def Smirnov(x,lambda_list):
    F=1
    kk=math.floor(len(lambda_list)/2)
    for k in range(1,kk+1):
        F=F+(-1)**k/pi*integrate.quad(lambda la: exp(-la*x/2)*(-Determ(la,lambda_list))**(-0.5)/la, lambda_list[2*k-2], lambda_list[2*k-1])[0]
    return F




#p-value
p_value=1-Smirnov(omega_2,lambda_list)        
print('p-value ',p_value)

