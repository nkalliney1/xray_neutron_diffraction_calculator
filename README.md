# Diffraction Plotter

## How to use the code

### Factors under consideration
### Types of diffraction
This code can create diffraction plots for four types of diffraction. The type of diffraction must be specified by the user with a flag.

(a) x-ray diffraction without angular dependence (flag: `-xz`)

(b) x-ray diffraction with angular dependence (cromer form factors) (flag: `-xc`)

(c) nuclear scattering (flag: `-n`)

(d) magnetic nuclear scattering (flag: `-nm`)


### File types
This code can read crystal structures in `.cif` and `.vasp` files. The default is `.cif` files, but `.vasp` files can be used by adding the flag `-v`.

### Extra factors
Diffraction is impacted by a temperature factor. By default this is not included, but can be added to calculations if desired by including the flag `-t`.

Unpolarized x-ray diffraction is also impacted by an angular dependence called the Lorentz-Polarization factor. By default this is not included, but can be added to calculations if desired by including the flag `-lp`. This can also be added to other types of diffraction but is not physically meaningful.

### Wavelength
The wavelength of the scattering x-ray or neutron must be specified by the user.

### Partial occupancy
If the crystal has partial occupancy, add the flag `-po`.

### Command structure
When running the code, the terminal command follows this structure:
```
python3 diffraction.py crystal_name wavelength diffraction_type optional_flags
```
The crystal name should be given without the file extension, and wavelength should be in Angstroms.

Example command:
This command produces the magnetic neutron diffraction plot for a crystal titled BCC with partial occupancy. The wavelength of the incident neutrons is 1.54 A.
```
python3 diffraction.py BCC 1.54 -nm -po
```

## Files required
All crystals must be described using a `.cif` or `.vasp` file.

### Partial occupancy
Crystals that are described with partial occupancy should be structured like so:
In the `.cif` file, each location should be filled with a placeholder element, which does not necessarily have to be one of the elements in the crystal, but must be unique.
The user should create a second file, titled `crystalname-header.csv` that should be structured like so:
```
placeholder element_1 occupancy_frac1 element_2 occupancy_frac2 [etc]
```
for as many lines as necessary. Lines can be in any order.

For example, to describe AlCrVTi as a BCC lattice where the center position is occupied by V and Ti at equal fractions, and the vertex occupied by Al and Cr at equal fractions, and using Al as the placeholder for the vertices and Ti as the placeholder for the center:

End of `.cif` file:
```
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_fract_symmform
Al1 Al  1 a 0.00000 0.00000 0.00000 1.00000 0,0,0 
Ti1 Ti  1 b 0.50000 0.50000 0.50000 1.00000 0,0,0 
```

Header file:
```
Al Al 0.5 Cr 0.5
Ti Ti 0.5 V 0.5
```

### Magnetic neutron diffraction
For magnetic neutron diffraction, all atoms with magnetic moments should have their moments specified in any order n a file named `crystalname-moments.csv`. Moments should be in terms of the Bohr magneton. 

For example, for a crystal where Fe is the only atom with a magnetic moment, `Fe-moments.csv` is simply:
```
Fe 2
```

## Theoretical backing

### X-ray and nuclear neutron diffraction
The Bragg condition requires that $$\vec{\Delta k} = \vec{G} = h\mathbf{b_1} + k\mathbf{b_2} + l\mathbf{b_3}$$, where $$\mathbf{b_i}$$ are the recriprocal lattice vectors.

The intensity of diffraction then is proportional to the squared magnitude of the structure factor. 
$$I \propto |S_G(h, k, l)|^2$$

The structure factor is given by summing over the atoms in the unit cell.

$$S_G(\vec{\Delta k}) = \sum_{j} f_j exp(-2\pi i (hx_j + ky_j + lz_j))$$.

For nuclear neutron scattering, $$f_j$$ does not have an angular dependence, and is simply the atomic neutron scattering length.

For x-ray diffraction, the form factor is given by $$f_j(s) = \sum_{i=1}^{4}[a_ie^{-b_i s^2}] + c$$.

$$s = sin(\theta)/\lambda$$ where $$\theta$$ is the incident angle of the beam and $$\lambda$$ is its wavelength. The $$a_i$$ and $$b_i$$ are the Cromer factors and are element-specific.

### Elastic magnetic neutron scattering
Intensity of magnetic neutron scatting is given by $$I = I_n + I_m$$, where $$I_n$$ is nuclear neutron scattering (see above) and $$I_m$$ is magnetic neutron scattering.

The structure factor for magnetic neutron scattering is a vector because the relative orientation of the neutron and atom's magnetic moments affects the scattering. 
Recall that $$\vec{\Delta k} = h\mathbf{b_1} + k\mathbf{b_2} + l\mathbf{b_3}$$

The structure factor is given by is given by $$\vec{F_m}(\vec{\Delta k})= \sum_{j} \vec{\mu_j} f_j(s) exp(-2\pi i (hx_j + ky_j + lz_j)) $$.

The $$f_j(s)$$ are given by $$f_j(s) = <j_0> + \frac{g-2}{2}<j_2>$$ where $$g$$ is the neutron's g-factor and $$s$$ is defined as above.

$$<j_0>= A_0e^{-a_0s^2} + B_0e^{-b_0s^2} + C_0e^{-c_0s^2} + D_0$$ where $$A_0,a_0,B_0,b_0,C_0,c_0,D_0$$ are element specific factors.


$$<j_2>= [A_2e^{-a_2s^2} + B_2e^{-b_2s^2} + C_2e^{-c_2s^2} + D_2]s^2$$ where $$A_2,a_2,B_2,b_2,C_2,c_2,D_2$$ are element specific factors.

$$\vec{\mu_j}$$ is the atom's magnetic moment in terms of the Bohr magneton.

The intensity is given by $$I_m(\Delta\vec{k}) = |\hat{\Delta k} \times \vec{F_m}(\vec{\Delta k}) \times \hat{\Delta k}|^2$$
## Sources for coefficients
### Cromer factors
[Database](https://sources.debian.org/src/libxray-scattering-perl/3.0.1-2/data/cromann.dat/)

## Magnetic neutron scattering factors
[Database](https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html)

## Neutron scattering lengths
Used bound coherent scattering length.

[Database](https://www.ncnr.nist.gov/resources/n-lengths/list.html)
