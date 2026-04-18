import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from PBE_diffraction_module import get_reciprocal_vectors, calculate_I, get_crystal, get_form_factor_array, get_j0_coeffs, get_j2_coeffs

name = sys.argv[1]
wavelength = float(sys.argv[2])

file_type = "f"

diffraction_type = ""
if "-xz" in sys.argv:
    diffraction_type = "xz"
elif "-xc" in sys.argv:
    diffraction_type = "xc"
elif "-n" in sys.argv:
    diffraction_type = "n"
elif "-nm" in sys.argv:
    diffraction_type = "nm"
elif "-m" in sys.argv:
    diffraction_type = "m"
else:
    raise ValueError("Need valid diffraction type")

p_on = False
t_on = False
l_on = False
partial_occupancy = False
occupancies = {}
if "-p" in sys.argv:
    p_on = True
if "-l" in sys.argv:
    l_on = True
if "-lp" in sys.argv:
    l_on = True
    p_on = True
if "-t" in sys.argv:
    t_on = True


crystal = get_crystal(name)
form_factor_array = get_form_factor_array(crystal, diffraction_type, partial_occupancy, occupancies)
rvectors = get_reciprocal_vectors(crystal)
j0_coeffs = 0
j2_coeffs = 0
if diffraction_type == "nm" or diffraction_type == "m":
    j0_coeffs = get_j0_coeffs(crystal, partial_occupancy, occupancies)
    j2_coeffs = get_j2_coeffs(crystal, partial_occupancy, occupancies)
    

two_theta = []
I_G = []
labels = []
I_G_g = []
two_theta_g = []
to_write = "h,k,l,|G|,2theta,d,I\n"

#diffraction calculations
lower = -9
upper = 9

for h in range(lower,upper+1):
    for k in range(lower,upper+1):
            for l in range(lower,upper+1):
                if (h == 0 and l == 0 and k == 0):
                    continue

                G = h*rvectors[0]+k*rvectors[1]+l*rvectors[2]
                hkl = np.array([h, k, l])

                mag_k = 2*math.pi/wavelength
                mag_G = math.sqrt(np.dot(G, G))

                OH = 0.5*math.sqrt(np.dot(G, G))/mag_k
                #test if 2theta is in range
                if (OH < -1 or OH > 1):
                    continue
                else:
                    theta_r = math.asin(OH)

                theta = math.degrees(theta_r)
                L = 1
                T = 1
                P = 1
                if l_on:
                    L = 1/((math.sin(theta_r)**2)*4*math.cos(theta_r))
                if p_on:
                    P = 0.5*(1+math.cos(2*theta_r)**2)
                if t_on:
                    B = 0.8
                    T = math.exp(-B*(math.sin(theta_r)/wavelength)**2)
                

                I = calculate_I(diffraction_type, form_factor_array, crystal, hkl, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies, L, P, T)
                #makes data a bit cleaner by removing rounding errors
                if (I < 0.000000000000001):
                    I = 0

                #hkl
                to_write+=str(h) + "," + str(k) + "," + str(l) + ","
                #|G|
                to_write+=str(mag_G) + ","
                #2theta
                to_write+=str(2*theta) + ","
                #interplanar spacing
                to_write+=str(2*math.pi/np.dot(G, G)) + ","
                #intensity
                to_write+=str(I) + "\n"
                

                #for plotting sharp peaks
                if (2*theta in two_theta):
                    I_G[two_theta.index(2*theta)] += I
                elif 2*theta > 0:
                    two_theta.append(2*theta)
                    I_G.append(I)
                
                #for plotting with gaussian
                mu = mag_G
                sigma = 0.01
                x = np.linspace(mu - 5*sigma, mu + 5*sigma, 101)
                for i in range(len(x)):
                    if (0.5*x[i]/mag_k > 1 or 0.5*x[i]/mag_k<-1):
                        continue
                    theta_2 = math.degrees(math.asin(0.5*x[i]/mag_k))
                    if 2*theta_2 in two_theta_g:
                        I_G_g[two_theta_g.index(2*theta_2)] += I * math.exp(-(mag_G-x[i])**2/(2*sigma**2))
                    elif 2*theta_2 > 0:
                        two_theta_g.append(2*theta_2)
                        I_G_g.append(I * math.exp(-(mag_G-x[i])**2/(2*(sigma**2)))) 
                
#make a function of just adding gaussians times data points centered at the peaks from 0-120
#in a different file

max = 0
for i in range(len(I_G_g)):
    if I_G_g[i] > max and two_theta_g[i] > 0 and two_theta_g[i] < 120:
        max = I_G_g[i]
print(max)
for i in range(len(I_G_g)):
    I_G_g[i] = I_G_g[i]/max

with open("./" + name + "/"+ diffraction_type+"_" + name + "_data.csv", "w") as f:
    f.write(to_write)

with open("./"+ name + "/" + diffraction_type+ "_" + name +"_peaks.csv", "w") as f:
    to_write2 = "2theta, I\n"
    for i in range(len(two_theta)):
        if I_G[i] > 0:
            to_write2 += str(two_theta[i]) + "," + str(I_G[i])+"\n"
    f.write(to_write2)
plt.title("My code")
plt.plot(two_theta_g, I_G_g)
plt.xlim(0, 120)
plt.ylim(0, 1.1)
plt.xlabel(r"2${\Theta}$ [deg]")
plt.ylabel(r"Intensity")
plt.show()
