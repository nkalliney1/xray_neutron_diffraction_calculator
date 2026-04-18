import math
import cmath
import numpy as np
import ase.io.vasp
import ase.io.cif
import csv

def get_atomic_symbol(num):
    with open("./params/Z_list.txt", "r") as f:
        lines = f.readlines()
        symbol = lines[num].split(" ")[1][:-1]
        return symbol

def get_atomic_num(symbol):
    with open("./params/Z_list.txt", "r") as f:
        lines = f.readlines()
        num = 1
        current = lines[1].split(" ")[1][:-1]
        while (current != symbol):
            num+=1
            current = lines[num].split(" ")[1][:-1]
        return num

def get_cromer_f(form_factors, symbol, wavelength, theta_r):
    #get cromer form factors as a function of sine(theta)/wavelength
    nums = form_factors["c"][symbol]
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
    f = j0 + 0*j2
    return f

def get_magnetic_moment(symbol, moments):
    if symbol in moments:
        return moments[symbol]*np.array([1,0,0])
    else:
        return 0
    
def get_moments(name):
    #get all magnetic moments
    moments = {}
    with open("./crystals/" + name + "-moments.csv") as f:
        lines = f.readlines()
        for line in lines:
            line_split = line.split(" ")
            moments[line_split[0]] = float(line_split[1])
    return moments

def get_reciprocal_vectors(crystal):
    #reciprocal lattice vectors
    reciprocal_vectors = []
    cell = crystal.cell
    volume = np.dot(cell[0], np.cross(cell[1], cell[2]))
    reciprocal_vectors.append(2*math.pi*np.cross(cell[1],cell[2])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[2],cell[0])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[0],cell[1])/volume)
    return reciprocal_vectors

def calculate_I_z(crystal, v, partial_occupancy, occupancies):
    #calculate intensity for scattering w/o angular dependence

    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*get_atomic_num(atom[0])
        else:
            f = crystal.numbers[i]
        SG += f*e
    return abs(SG)**2

def calculate_I_c(form_factors, crystal, v, wavelength, theta_r, partial_occupancy, occupancies):
    #calculate intensity for scatting w cromer form factors
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*float(get_cromer_f(form_factors,atom[0],wavelength,theta_r))
        else:
            f = float(get_cromer_f(form_factors,crystal.symbols[i],wavelength,theta_r))

        SG += f*e
    return abs(SG)**2

def calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies):
    #nuclear neutron scattering
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*float(form_factors["n"][atom[0]])
        else:
            f = float(form_factors["n"][crystal.symbols[i]])

        SG += f*e
    return abs(SG)**2


def calculate_I_m(crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies):
    F_m = np.array([0,0,0])
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        
        f = 0
        m = 0
        if partial_occupancy:
            for atom in occupancies[crystal.symbols[i]]:
                if (v[0] == 1 and v[1] == 1 and v[2] == 0):
                    print(atom)
                    print(float(get_magnetic_f(atom[0],wavelength,theta_r, j0_coeffs, j2_coeffs)))
                    print(get_magnetic_moment(atom[0], moments))
                
                f2 = atom[1]*float(get_magnetic_f(atom[0],wavelength,theta_r, j0_coeffs, j2_coeffs))
                m2 = atom[1]*get_magnetic_moment(atom[0], moments)

                f += f2
                m += m2
        else:
            f = float(get_magnetic_f(crystal.symbols[i],wavelength,theta_r, j0_coeffs, j2_coeffs))
            m = get_magnetic_moment(crystal.symbols[i], moments)
        
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        if (v[0] == 1 and v[1] == 1 and v[2] == 0):
            print(f*e*m)
        F_m = np.add(F_m, f*e*m)
    if (v[0] == 1 and v[1] == 1 and v[2] == 0):
       print(F_m)
    k_hat = v / math.sqrt(np.dot(v, v))
    I = np.linalg.norm(np.cross(k_hat, np.cross(F_m, k_hat)))**2

    return I

def calculate_I(type, form_factors,crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies, L, P, T):
    I = 0

    #intensities for x ray must be multiplied by r_e^2 bc f is in terms of r_e normally
    if "z" in type:
        I += (2.8e-15)**2*calculate_I_z(crystal, v, partial_occupancy, occupancies)
    if "c" in type:
        I += (2.8e-15)**2*calculate_I_c(form_factors, crystal, v, wavelength, theta_r, partial_occupancy, occupancies)
    if "n" in type:
        I += calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies)
    #pre factor is normalization into barns
    #I_m has a fixed constant (yr/2) that we are dividing out (and changing units into barns)
    #y = neutron gyromagnetic ratio
    #r = classical electron radius
    if "m" in type:
        I += (1/0.066745253805)*calculate_I_m(crystal, v, wavelength, theta_r, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies)
    return L*P*T*I

def get_crystal(name, file_type):
    if file_type == "v":
        crystal = ase.io.vasp.read_vasp("./crystals/"+name+".vasp")
    elif file_type == "f":
        crystal = read_crystal_file(name)
    else:
        crystal = ase.io.cif.read_cif("./crystals/"+name+".cif")
    return crystal

def get_form_factor_array(crystal, diffraction_type, partial_occupancy, occupancies):
    #get form factors depending on type of diffraction
    #only for nuclear neutron or cromer form factors
    
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    form_factors = {}
    if "c" in diffraction_type:
        with open("./params/cromer_factors.csv", newline='') as f:
            cromer_form_factors = {}
            lines = f.readlines()
            for i in range(0,len(lines),2):
                line = lines[i][:-1]
                if line in symbols and line not in cromer_form_factors:
                    nums = lines[i+1][2:-1].split(" ")
                    cromer_form_factors[line]=[float(x) for x in nums]
            form_factors["c"] = cromer_form_factors

    if "n" in diffraction_type:
        with open("./params/neutron_diffraction_lengths.csv", 'r') as file:
            all_neutron_form_factors = []
            csvreader = csv.reader(file)
            all_neutron_form_factors = next(csvreader)
            neutron_form_factors = {}

            for i in range(len(all_neutron_form_factors)):
                if (get_atomic_symbol((i)) in symbols):
                    neutron_form_factors[get_atomic_symbol((i))] = float(all_neutron_form_factors[i-1])
            form_factors["n"] = neutron_form_factors
    
    return form_factors

def get_j0_coeffs(crystal, partial_occupancy, occupancies):
    #used to calculate nuclear scattering form factor
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    j0 = {}
    with open("./params/j0.csv", newline='') as f:
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
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    j2 = {}
    with open("./params/j2.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in symbols and line[0:2] not in j2:
                nums = lines[i][:-1].split("\t")[1:-1]
                j2[line[0:2]]=[float(x) for x in nums]
    return j2

def get_occupancies(name):
    occupancies = {}

    with open("./crystals/"+name+"-header.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            nums = lines[i][3:-1].split(" ")
            occupancies[line[0:2]] = []
            for j in range(0,len(nums),2):
                occupancies[line[0:2]].append([nums[j], float(nums[j+1])])
    return occupancies
