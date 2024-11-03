#!/usr/bin/env python
# coding: utf-8

# Our task is to find the roots of a given function in the complex plane via the Newton-Raphson method.
# The iteration scheme is the same as for functions in real space:
# \begin{equation}
#     z_{n+1}=z_n-\frac{f(z)}{f'(z)},\quad\text{where}\quad z\in\mathbb{C}.
# \end{equation}
# 
# Firstly, we consider the following complex-valued function:
# \begin{equation}
#     f(z)=z^3-1.
# \end{equation}
# 
# Thus, we are asked to solve the equation: $z^3-1=0$.
# 
# The iteration procedure in our case reads:
# \begin{equation}
#     z_{n+1}=z_n-\frac{z^3_n-1}{3z^2_n}.
# \end{equation}

# In[2]:


import cmath   # to deal with complex numbers
import numpy as np
from matplotlib import pyplot as plt 


# In the next section we define de following functions:
# - $\textbf{f}$ evaluates the complex function.
# - $\textbf{Newton\_Raphson}$ evaluates the recurrence relation.
# - $\textbf{solve\_cnewton}$ performs the Newton_Raphson method for a given initial point and returns the coordinates of the starting point, the root (or the last guess if it does not converge), the function evaluated at the root, the number of steps $n$, and $\log_{10}(n)$.

# In[4]:


def f(z):
    return z**3 - 1


def Newton_Raphson(z):
    z_new = z - (z**3 - 1)/(3*z**2)
    return z_new


def solve_cnewton(duplet):
    counter = 0
    
    x_duplet, y_duplet = duplet
    z = complex(x_duplet, y_duplet)   # make the duplet a point of the complex plane
    #print("z=",z)
    
    for i in range(M):
        z_new = Newton_Raphson(z)
        counter = counter + 1
        
        f_new = f(z_new)
    
        if f_new == 0 or abs(z_new - z) < tolerance:   # we are done
            return x_duplet, y_duplet, z_new, f_new, counter, np.log10(counter)
        elif i == M - 1:
            return x_duplet, y_duplet, z_new, f_new, counter, np.log10(counter)

        z = z_new


# In the following section we construct the grid starting from setting the upper left point and the lower right point. In particular we set (-2,2) and (2,-2) as such points.
# We find the root for each point of the grid and construct the scatter plot with the starting points. We mark the starting points with three different colors (black, green and red) on the basis of the kind of root that we find.

# In[6]:


tolerance = 1e-9
tolerance2 = 1e-6
M = 200   # maximal number of steps
N = 200   # dimension of the grid (NxN square matrix)

grid_point1 = (-2, 2)   # upper left
grid_point2 = (2, -2)   # lower right
x1, y1 = grid_point1
x2, y2 = grid_point2

x_values = np.linspace(x1, x2, N)   # create the two "axes" of the grid
y_values = np.linspace(y1, y2, N)

x, y = np.meshgrid(x_values, y_values)   # create coordinates of the grid points
#print("x=",x)
#print("y=",y)
grid = np.dstack((x, y))   # construct the grid, each point of it corresponds to a duplet (x, y)
#print(grid)
#print(grid[N-1][N-1])

x_duplet, y_duplet, root, f_root, counter, log_counter = solve_cnewton(grid[2][1])   # test
print(x_duplet, y_duplet, root, f_root, counter, log_counter)

log_n_matrix = np.zeros((N,N))   # initialize the matrix for storing log10(n)

plt.figure()
for i in range(N):
    for j in range(N):
        x_duplet, y_duplet, root, f_root, counter, log_counter = solve_cnewton(grid[i][j])

        log_n_matrix[i][j] = log_counter 
        
        if abs(root.real - 1) < tolerance2:   # root 1
            plt.scatter(x_duplet, y_duplet, c="k", marker=".")
        elif abs(root.imag - 0.8660254040) < tolerance2:   # root (0.5, 0.8660254040)
            plt.scatter(x_duplet, y_duplet, c="g", marker=".")
        else:   # root (0.5, -0.8660254040)
            plt.scatter(x_duplet, y_duplet, c="r", marker=".")
        '''
        print(f"x0={x_duplet}, y0={y_duplet}, root={root}, f(root)={f_root}, counter={counter}, log={log_counter}")
        print("")
        '''
plt.xlabel("x")
plt.ylabel("y")
plt.title("Black: 1; green: (0.5, 0.86602); red: (0.5, -0.86602)")


# We now plot $\log_{10}(n)$ in the ($x$,$y$)-plane.

# In[8]:


plt.figure()
plt.pcolormesh(x, y, log_n_matrix, shading="auto", cmap="viridis")
plt.colorbar(label="log10(n)")
plt.xlabel("x")
plt.ylabel("y")


# In this section we repeat the entire procedure zooming in around a specific zone in order to highlight the fractal structure.
# In particular we set (0.355,-1.010) and (0.371,-1.045) as the new upper left and lower right points respectively.

# In[10]:


grid_point1 = (0.355, -1.010)   # upper left
grid_point2 = (0.371, -1.045)   # lower right
x1, y1 = grid_point1
x2, y2 = grid_point2

x_values = np.linspace(x1, x2, N)   # create the two "axes" of the grid
y_values = np.linspace(y1, y2, N)

x, y = np.meshgrid(x_values, y_values)   # create coordinates of the grid points

grid = np.dstack((x, y))   # construct the grid, each point of it corresponds to a duplet (x, y)

log_n_matrix = np.zeros((N,N))   # initialize the matrix for storing log10(n)

plt.figure()
for i in range(N):
    for j in range(N):
        x_duplet, y_duplet, root, f_root, counter, log_counter = solve_cnewton(grid[i][j])
        
        if abs(root.real - 1) < tolerance2:   # root 1
            plt.scatter(x_duplet, y_duplet, c="k", marker=".")
        elif abs(root.imag - 0.8660254040) < tolerance2:   # root (0.5, 0.8660254040)
            plt.scatter(x_duplet, y_duplet, c="g", marker=".")
        else:   # root (0.5, -0.8660254040)
            plt.scatter(x_duplet, y_duplet, c="r", marker=".")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Black: 1; green: (0.5, 0.86602); red: (0.5, -0.86602)")


# In[11]:


def f2(z):
    return 35 * z**9 - 180 * z**7 + 378 * z**5 - 420 * z**3 + 315 * z


def Newton_Raphson2(z):
    z_new = z - (35 * z**9 - 180 * z**7 + 378 * z**5 - 420 * z**3 + 315 * z)/(35 * 9 * z**8 - 180 * 7 * z**6 + 378 * 5 * z**4 - 420 * 3 * z**2 + 315)
    return z_new


def solve_cnewton2(duplet):
    counter = 0
    
    x_duplet, y_duplet = duplet
    z = complex(x_duplet, y_duplet)   # make the duplet a point of the complex plane
    #print("z=",z)
    
    for i in range(M):
        z_new = Newton_Raphson2(z)
        counter = counter + 1
        
        f_new = f2(z_new)
    
        if f_new == 0 or abs(z_new - z) < tolerance:   # we are done
            return x_duplet, y_duplet, z_new, f_new, counter, np.log10(counter)
        elif i == M - 1:
            return x_duplet, y_duplet, z_new, f_new, counter, np.log10(counter)

        z = z_new


# In[12]:


tolerance = 1e-9
tolerance2 = 1e-6
M = 200   # maximal number of steps
N = 200   # dimension of the grid (NxN square matrix)

grid_point1 = (-2, 2)   # upper left
grid_point2 = (2, -2)   # lower right
x1, y1 = grid_point1
x2, y2 = grid_point2

x_values = np.linspace(x1, x2, N)   # create the two "axes" of the grid
y_values = np.linspace(y1, y2, N)

x, y = np.meshgrid(x_values, y_values)   # create coordinates of the grid points
#print("x=",x)
#print("y=",y)
grid = np.dstack((x, y))   # construct the grid, each point of it corresponds to a duplet (x, y)
#print(grid)
#print(grid[N-1][N-1])

x_duplet, y_duplet, root, f_root, counter, log_counter = solve_cnewton2(grid[2][1])   # test
print(x_duplet, y_duplet, root, f_root, counter, log_counter)

log_n_matrix = np.zeros((N,N))   # initialize the matrix for storing log10(n)

#plt.figure()
for i in range(N):
    for j in range(N):
        x_duplet, y_duplet, root, f_root, counter, log_counter = solve_cnewton2(grid[i][j])

        log_n_matrix[i][j] = log_counter 
        
        if abs(root.real - 0) < tolerance2:   # root 0
            plt.scatter(x_duplet, y_duplet, c="k", marker=".")
        elif abs(root.real - 0.93774544) < tolerance2 and abs(root.imag - 0.65437520) < tolerance2:   # root (0.93774544, 0.65437520)
            plt.scatter(x_duplet, y_duplet, c="g", marker=".")
        elif abs(root.real - 0.93774544) < tolerance2 and abs(root.imag + 0.65437520) < tolerance2:   # root (0.93774544, -0.65437520)
            plt.scatter(x_duplet, y_duplet, c="r", marker=".")
        elif abs(root.real + 0.93774544) < tolerance2 and abs(root.imag - 0.65437520) < tolerance2:   # root (-0.93774544, 0.65437520)
            plt.scatter(x_duplet, y_duplet, c="darkorchid", marker=".")
        elif abs(root.real + 0.93774544) < tolerance2 and abs(root.imag + 0.65437520) < tolerance2:   # root (-0.93774544, -0.65437520)
            plt.scatter(x_duplet, y_duplet, c="orange", marker=".")
        elif abs(root.real + 1.48569) < tolerance2 and abs(root.imag - 0.295006) < tolerance2:   # root (-1.48569, 0.295006)
            plt.scatter(x_duplet, y_duplet, c="yellow", marker=".")
        elif abs(root.real + 1.48569) < tolerance2 and abs(root.imag + 0.295006) < tolerance2:   # root (-1.48569, -0.295006)
            plt.scatter(x_duplet, y_duplet, c="cyan", marker=".")
        elif abs(root.real - 1.48569) < tolerance2 and abs(root.imag - 0.295006) < tolerance2:   # root (1.48569, 0.295006)
            plt.scatter(x_duplet, y_duplet, c="w", marker=".")
        else:   # root (1.48569, -0.295006)
            plt.scatter(x_duplet, y_duplet, c="pink", marker=".")
        '''
        print(f"x0={x_duplet}, y0={y_duplet}, root={root}, f(root)={f_root}, counter={counter}, log={log_counter}")
        print("")
        '''


# In[ ]:


plt.figure()
plt.pcolormesh(x, y, log_n_matrix, shading="auto", cmap="viridis")
plt.colorbar(label="log10(n)")
plt.xlabel("x")
plt.ylabel("y")

