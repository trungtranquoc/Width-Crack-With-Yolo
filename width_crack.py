import numpy as np
import matplotlib.pyplot as plt
import functools
import math
from visibility import visibility_polygon

import cv2
from skimage.morphology import skeletonize
import sknw
dummy_ske = np.zeros([4, 4, 3], dtype=np.uint8)
sknw.build_sknw(dummy_ske)

# from lib_computing_crack import *

largest_range = 1000000000
delta = 0.0000005

#Solve Equation form x^3 + a.x^2 + b.x + c = 0
def solve(a, b, c, d):

    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)
        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j

        return np.array([x1, x2, x3])           # Returning One Real Root and two Complex Roots as numpy array.


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0


# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0


# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)

#End solve cubic-equation
     

def onSegment(p1, p2, p):
    """
    Check if the point p lies in the axis-aligned segment with endpoints p1, p2
    """
    a = (p1[0] - p[0], p1[1] - p[1])
    b = (p2[0] - p[0], p2[1] - p[1])    
    x1 = a[0] * b[1] 
    x2 = a[1] * b[0] 
    return (x1 >= x2 - delta and x1 <= x2 + delta) and p[0] >= min(p1[0] - delta,p2[0] - delta) and p[0] <= max(p1[0] + delta,p2[0] + delta) and p[1] >= min(p1[1] - delta,p2[1] - delta) and p[1] <= max(p1[1] + delta,p2[1] + delta)

#function calculate intersection
def line_intersection(p1, p2, q1, q2):
    """
    Find the intersection between the line containing p1, p2 and the line containing
    q1, q2
    """
    xdiff = (float(p1[0] - p2[0]), float(q1[0] - q2[0]))          #directed vector
    ydiff = (float(p1[1] - p2[1]), float(q1[1] - q2[1]))
    
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)

    if(div == 0):
        return (largest_range, largest_range)

    d = (det(p1, p2), det(q1, q2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return (x, y)

def angle_between_2_lines(p1,p2,q1,q2):
    xdiff = (p1[0] - p2[0], p1[1] - p2[1])          
    ydiff = (q1[0] - q2[0], q1[1] - q2[1])

    angle_diff = (xdiff[0] * ydiff[1] - xdiff[1] * ydiff[0])/math.sqrt((xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1]) * (ydiff[0] * ydiff[0] + ydiff[1] * ydiff[1]))
    if abs(angle_diff) > 1:
        return math.pi/2
    else:
        return abs(math.asin(angle_diff))



#function check whether there is segments go through p with 2 endpoints in 2 lines. (Lemma 1)
def check_visible(p1, p2, q1, q2, p):
    """
    Check whether 4 line going through p to p1,p2,q1,q2 passed opposite edge
    """

    # Check whether c lies between a and b with condition a,b,c are in straight line.
    def lies_between(a,b,c):
        return c[0] >= min(a[0] - delta,b[0] - delta) and c[0] <= max(a[0] + delta,b[0] + delta) and c[1] >= min(a[1] - delta,b[1] - delta) and c[1] <= max(a[1] + delta,b[1] + delta)

    if p2[0] == q1[0] and p2[1] == q1[1]:
        alpha1 = angle_between_2_lines(p,p1,p,p2)
        alpha2 = angle_between_2_lines(p,q1,p,q2)
        return (alpha1 + alpha2 >= math.pi) 

    a1 = line_intersection(p,p1,q1,q2)
    if a1[0] != largest_range and a1[1] != largest_range and lies_between(q1,q2,a1) and lies_between(a1,p1,p):
        return True
    a2 = line_intersection(p,p2,q1,q2)
    if a2[0] != largest_range and a2[1] != largest_range and lies_between(q1,q2,a2) == True and lies_between(p2,a2,p):
        return True
    b1 = line_intersection(p,q1,p1,p2)
    if(b1[0] != largest_range and b1[1] != largest_range and lies_between(p1,p2,b1) == True and lies_between(q1,b1,p)):
        return True
    b2 = line_intersection(p,q2,p1,p2)
    if(b2[0] != largest_range and b2[1] != largest_range and lies_between(p1,p2,b2) == True and lies_between(q2,b2,p)):
        return True

    return False

#Distance between 2 points:
def distance(a,b):
    xdiff = a[0] - b[0]
    ydiff = a[1] - b[1]
    return abs(math.sqrt(xdiff * xdiff + ydiff * ydiff))



#return shortest line segment that go through p and joint 2 edges
def shortest_edge(p1,p2,q1,q2,p):
    #Determine whether intersection of line cd and ab is lies on ab ? 
    def check_intersection(a,b,c,d):
        intersect = line_intersection(a,b,c,d)
        if(intersect[0] == largest_range):
            return False
        else:
            return (onSegment(a,b,intersect))

    #Step check: 
    if(check_visible(p1,p2,q1,q2,p) == False):
        return np.array([(-largest_range, -largest_range),(largest_range,largest_range)])

    def projection_point(x,y,p):
        xdiff = x[0] - y[0]
        ydiff = x[1] - y[1]
        q = (p[0]+ydiff, p[1]-xdiff)
        return line_intersection(x,y,p,q) 
    
    if ((p1[0]-p2[0])*(q1[1]-q2[1]) == (p1[1]-p2[1])*(q1[0]-q2[0])):
        #Calculate point i and j:
        i = projection_point(p1,p2,p)
        j = projection_point(q1,q2,p)
    else:
        #Intersection of line ab and line D through p and parallel cd:
        def line_intersection_of_parallel(a,b,c,d,p):
            xdiff = c[0] - d[0]
            ydiff = c[1] - d[1]
            p_new = (p[0]+xdiff, p[1]+ydiff)
            return line_intersection(a,b,p,p_new) 

        o = line_intersection(p1,p2,q1,q2)
        p_p = line_intersection_of_parallel(p1,p2,q1,q2,p)
        p_q = line_intersection_of_parallel(q1,q2,p1,p2,p)
        alpha = angle_between_2_lines(o,p_p,o,p_q)
        p_2 = distance(p, p_p) 
        p_1 = distance(p, p_q) 

        a2 = -(p_2 * math.cos(alpha))
        a3 = p_1 * p_2 * math.cos(alpha)
        a4 = -(p_1 * p_2 * p_2)
        array = solve(1, a2, a3, a4)

        for t in range (0,3):
            if(array[t].imag == 0 and array[t].real > 0):
                break
        t = array[t].real 

        #Calculate point i and j:
        ratio = (p_1 + t) / distance(p_p,o)
        i = (o[0] + ratio*(p_p[0] - o[0]), o[1] + ratio*(p_p[1] - o[1]))
        j = line_intersection(p,i,q1,q2) 

        # return np.array([p_p,p_q])
        # return np.array([i,j])
        a = (largest_range, largest_range)
    if i[0] <= max(p1[0] + delta,p2[0] + delta) and i[0] >= min(p1[0] - delta,p2[0] - delta) and j[0] <= max(q1[0] + delta,q2[0] + delta) and j[0] >= min(q1[0] - delta,q2[0] - delta) and i[1] <= max(p1[1] + delta,p2[1] + delta) and i[1] >= min(p1[1] - delta,p2[1] - delta) and j[1] <= max(q1[1] + delta,q2[1] + delta) and j[1] >= min(q1[1] - delta,q2[1] - delta) :
        return np.array([i,j])

    min_value = largest_range

    if(check_intersection(p1,p2,q1,p)):
        c_img = line_intersection(p1,p2,q1,p)
        if(distance(q1,c_img) < min_value and onSegment(q1,c_img,p)):
            min_value = distance(q1,c_img)
            i = c_img
            j = q1

    if(check_intersection(p1,p2,q2,p)):
        d_img = line_intersection(p1,p2,q2,p)
        if(distance(q2,d_img) < min_value and onSegment(q2,d_img,p)):
            min_value = distance(q2,d_img)
            i = d_img
            j = q2

    if(check_intersection(q1,q2,p1,p)):
        a_img = line_intersection(q1,q2,p1,p)
        if(distance(p1,a_img) < min_value and onSegment(p1,a_img,p)):
            min_value = distance(p1,a_img)
            i = p1
            j = a_img

    if(check_intersection(q1,q2,p2,p)):
        b_img = line_intersection(q1,q2,p2,p)
        if(distance(p2,b_img) < min_value and onSegment(p2,b_img,p)):
            min_value = distance(p2,b_img)
            i = p2
            j = b_img

    return np.array([i,j])

#Main algorithm
def width_crack(p, pol):
    # Check if line that go through a and b pass p
    def check_straight(a,b,p):
        return (onSegment(a,p,b) or onSegment(b,p,a))

    m = len(pol)-1
    count = 0
    k = 0
    t = []

     # Check if p is in the boundary
    while count <= m:
        if onSegment(pol[count], pol[(count + 1) % m], p):
            return (0,pol[count], pol[(count + 1) % m], p)
        if (check_straight(pol[count], pol[(count+1) % m], p) == False):
            t.append([pol[count], pol[(count+1) % m]])
            k = k+1 
        count += 1 

    i_o = 1
    while check_visible(t[0][0], t[0][1], t[i_o][0], t[i_o][1], p) == False:
        i_o = i_o + 1 

    min_value = largest_range

    i = 0
    j = i_o
    # test_set = (t[i][0], t[i][1], t[j][0], t[j][1])
    while i<=i_o and j<k:
        array = shortest_edge(t[i][0], t[i][1], t[j][0], t[j][1], p)
        if(distance(array[0], array[1]) < min_value):
            min_value = distance(array[0], array[1])
            min_x = array[0]
            min_y = array[1]
            test_set = (t[i][0], t[i][1], t[j][0], t[j][1])
        if(j < k-1 and check_visible(t[i][0], t[i][1], t[j+1][0], t[j+1][1], p) == True):
            j = j+1
        else:
            i = i+1

    return min_value, min_x, min_y
    # return [min_x, min_y]




