import math
import cmath
import numpy as np
import ase.io.vasp
import ase.io.cif
import csv

def get_atomic_symbol(num):
    with open("params/Z_list.txt", "r") as f:
        lines = f.readlines()
        symbol = lines[num].split(" ")[1][:-1]
        return symbol

def get_atomic_num(symbol):
    with open("params/Z_list.txt", "r") as f:
        lines = f.readlines()
        num = 1
        current = lines[1].split(" ")[1][:-1]
        while (current != symbol):
            num+=1
            current = lines[num].split(" ")[1][:-1]
        return num

def get_cromer_f(form_factors, symbol, wavelength, theta_r):
    #get cromer form factors as a function of sine(theta)/wavelength
    nums = form_factors[symbol]
    form_factor = 0
    s = math.sin(theta_r)/wavelength
    for i in range(4):
        ai = float(nums[2*i])
        bi = float(nums[2*i+1])
        form_factor += ai * math.exp(-bi*(s**2))
    form_factor += float(nums[8])
    return form_factor

def get_magnetic_f(symbol, wavelength, theta_r, j0_coeffs, j2_coeffs):
    #get magnetic form factors as a function of sine(theta)/wavelength
    s = math.sin(theta_r)/wavelength
    j0 = 0
    j2 = 0

    if symbol in j0_coeffs:
        j0_atom = j0_coeffs[symbol]
        for i in range(3):
            j0 += j0_atom[2*i] * math.exp(-1*j0_atom[2*i+1]*s*s)
        j0 += float(j0_atom[6])

    if symbol in j2_coeffs:
        j2_atom = j2_coeffs[symbol]
        for i in range(3):
            j2 += float(j2_atom[2*i]) * math.exp(-1*float(j2_atom[2*i+1])*s*s)
        j2 += float(j2_atom[6])
        j2 = j2*s*s
    
    #j2 multiplied by 0 for transition metals!
    f = j0 + 0*j2/2

    return f

def get_reciprocal_vectors(crystal):
    #reciprocal lattice vectors
    reciprocal_vectors = []
    cell = crystal["vectors"]
    volume = np.dot(cell[0], np.cross(cell[1], cell[2]))
    reciprocal_vectors.append(2*math.pi*np.cross(cell[1],cell[2])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[2],cell[0])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[0],cell[1])/volume)
    return reciprocal_vectors

def calculate_I_xz(crystal, v, partial_occupancy, occupancies):
    #calculate intensity for scattering w/o angular dependence

    SG = 0
    positions = crystal["positions"]
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        f = get_atomic_num(crystal["symbols"][i])
        SG += f*e*crystal["occupancies"][i]
    return abs(SG)**2

def calculate_I_xc(form_factors, crystal, v, wavelength, theta_r, partial_occupancy, occupancies):
    #calculate intensity for scatting w cromer form factors
    SG = 0
    positions = crystal["positions"]
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        f = float(get_cromer_f(form_factors,crystal["symbols"][i],wavelength,theta_r))

        SG += f*e*crystal["occupancies"][i]
    return abs(SG)**2

def calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies):
    #nuclear neutron scattering
    SG = 0
    positions = crystal["positions"]
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        f = float(form_factors[crystal["symbols"][i]])

        SG += f*e*crystal["occupancies"][i]
    return abs(SG)**2

def calculate_I_nm(form_factors, crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies):
    #nuclear and magnetic neutron scattering
    # 45.3484518099
    
    I = calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies)
    I += (1/0.066745253805)*calculate_I_m(crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies)
    return I

def calculate_I_m(crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies):
    F_m = np.array([0,0,0])
    positions = crystal["positions"]
    for i in range(len(positions)):
        
        f = 0
        m = 0
        f = float(get_magnetic_f(crystal["symbols"][i],wavelength,theta_r, j0_coeffs, j2_coeffs))
        m = crystal["moments"][i]
        
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        F_m = np.add(F_m, f*e*m*crystal["occupancies"][i])
    k_hat = v / math.sqrt(np.dot(v, v))
    I = np.linalg.norm(np.cross(k_hat, np.cross(F_m, k_hat)))**2

    return I

def calculate_I(type, form_factors,crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies, L, P, T):
    if type == "xz":
        return L*P*T*calculate_I_xz(crystal, v, partial_occupancy, occupancies)
    elif type == "xc":
        return L*P*T*calculate_I_xc(form_factors, crystal, v, wavelength, theta_r, partial_occupancy, occupancies)
    elif type == "n":
        return L*P*T*calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies)
    elif type == "nm":
        return L*P*T*calculate_I_nm(form_factors, crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies)
    elif type == "m":
        return L*P*T*calculate_I_m(crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, partial_occupancy, occupancies)


def get_form_factor_array(crystal, diffraction_type, partial_occupancy, occupancies):
    #get form factors depending on type of diffraction
    #only for nuclear neutron or cromer form factors
    symbols = crystal["symbols"]

    form_factors = []
    if diffraction_type == "xc":
        with open("params/cromer_factors.csv", newline='') as f:
            cromer_form_factors = {}
            lines = f.readlines()
            for i in range(0,len(lines),2):
                line = lines[i][:-1]
                if line in symbols and line not in cromer_form_factors:
                    nums = lines[i+1][2:-1].split(" ")
                    cromer_form_factors[line]=[float(x) for x in nums]
            form_factors = cromer_form_factors

    elif diffraction_type == "n" or diffraction_type == "nm":
        with open("params/neutron_diffraction_lengths.csv", 'r') as file:
            all_neutron_form_factors = []
            csvreader = csv.reader(file)
            all_neutron_form_factors = next(csvreader)
            neutron_form_factors = {}

            for i in range(len(all_neutron_form_factors)):
                if (get_atomic_symbol((i)) in symbols):
                    neutron_form_factors[get_atomic_symbol((i))] = float(all_neutron_form_factors[i-1])
            form_factors = neutron_form_factors
    elif diffraction_type == "xz":
        return
    elif diffraction_type == "m":
        return
    else:
        raise ValueError('Not a valid type of diffraction for this code')
    
    return form_factors

def get_j0_coeffs(crystal, partial_occupancy, occupancies):
    #used to calculate nuclear scattering form factor
    symbols = crystal["symbols"]

    j0 = {}
    with open("params/j0.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            #for two letter element names 
            if line[0:2] in symbols and line[0:2] not in j0:
                nums = lines[i][4:-1].split("\t")
                j0[line[0:2]]=[float(x) for x in nums]
            #for single-letter element name
            elif line[0:1] in symbols and line[0:1] not in j0:
                nums = lines[i][:-1].split("\t")[1:-1]
                j0[line[0:1]]=[float(x) for x in nums]
    return j0

def get_j2_coeffs(crystal, partial_occupancy, occupancies):
    #used to calculate nuclear scattering form factor
    symbols = crystal["symbols"]

    j2 = {}
    with open("params/j2.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in symbols and line[0:2] not in j2:
                nums = lines[i][:-1].split("\t")[1:-1]
                j2[line[0:2]]=[float(x) for x in nums]
    return j2

def get_crystal(name):
    crystal = {}
    with open("crystals/"+name+".csv") as f:
        lines = f.readlines()

        #get cell vectors, assumes they're at 90 degree angles
        v1 = np.array([float(lines[0].split(" ")[1]), 0, 0])
        v2 = np.array([0, float(lines[1].split(" ")[1]), 0])
        v3 = np.array([0, 0, float(lines[2].split(" ")[1])])
        crystal["vectors"] = [v1, v2, v3]

        #get everything else
        crystal["symbols"] = []
        crystal["positions"] = []
        crystal["occupancies"] = []
        crystal["moments"] = []
        for i in range(4, len(lines)):
            line = lines[i].replace("   ", " ").replace("  ", " ").split(" ")
            if len(line) == 1:
                continue
            
            crystal["symbols"].append(line[1])

            position = np.array([float(line[4]), float(line[5]), float(line[6])])
            crystal["positions"].append(position)
            crystal["occupancies"].append(float(line[7]))

            moment = np.array([float(line[8]), float(line[9]), float(line[10])])
            crystal["moments"].append(moment)

    return crystal