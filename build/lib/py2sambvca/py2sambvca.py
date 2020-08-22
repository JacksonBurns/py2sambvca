import subprocess, re, glob, os


class py2sambvca():
    """
    Wrapper class for py2sambvca functions.

    Call this class to instantiate a py2sambvca object, which has methods to write input, call SambVca,
    and retreieve output.

    Parameters:
    xyz_filepath (str): Location of .xyz molecular coordinates file for writing input data
    sphere_center_atom_ids (list): ID of atoms defining the sphere center 
    z_ax_atom_ids (list): ID of atoms for z-axis 
    xz_plane_atoms_ids (list): ID of atoms for xz-plane
    atoms_to_delete_ids (list): ID of atoms to be deleted (default None)
    sphere_radius (float): Sphere radius in Angstrom (default 3.5)
    displacement (float): Displacement of oriented molecule from sphere center in Angstrom (default 0.0)
    mesh_size (float): Mesh size for numerical integration (default 0.10)
    remove_H (int): 0/1 Do not remove/remove H atoms from Vbur calculation (default 1) 
    orient_z (int): 0/1 Molecule oriented along negative/positive Z-axis (default 1)
    write_surf_files (int): 0/1 Do not write/write files for top and bottom surfaces (default 1) 
    path_to_sambvcax (str): Path to the SambVca executable. Only needed to use py2sambvca.calc()( default "/path/to/executable/sambvca.x")
    """
    def __init__(self,
                xyz_filepath,
                sphere_center_atom_ids,
                z_ax_atom_ids,
                xz_plane_atoms_ids,
                atoms_to_delete_ids=None,
                sphere_radius=3.5,
                displacement=0.0,
                mesh_size=0.10,
                remove_H=1,
                orient_z=1,
                write_surf_files=1,
                path_to_sambvcax="/path/to/executable/sambvca.x"):
        """
        See docstring for py2sambvca
        """
        # if atoms are requested to be deleted, assign them and the number of them
        if atoms_to_delete_ids is not None:
            self.n_atoms_to_delete = len(atoms_to_delete_ids)
            self.atoms_to_delete_ids = atoms_to_delete_ids
        else:  # otherwise, set to none to avoid bad writes in the future
            self.n_atoms_to_delete = None
            self.atoms_to_delete_ids = None
        # various other parameters
        self.sphere_center_atom_ids = sphere_center_atom_ids
        self.n_sphere_center_atoms = len(sphere_center_atom_ids)
        self.z_ax_atom_ids = z_ax_atom_ids
        self.n_z_atoms = len(z_ax_atom_ids)
        self.xz_plane_atoms_ids = xz_plane_atoms_ids
        self.n_xz_plane_atoms = len(xz_plane_atoms_ids)
        self.sphere_radius = sphere_radius
        self.displacement = displacement
        self.mesh_size = mesh_size
        self.remove_H = remove_H
        self.orient_z = orient_z
        self.write_surf_files = write_surf_files
        # open the xyz file, read the data
        with open(xyz_filepath, "r") as file:
            self.xyz_data = file.readlines()
        # assign the path to the calculator
        self.path_to_sambvcax = path_to_sambvcax

    def write_input(self):
        """
        Write input for the Sambvca buried-volume Fortran calculator based on the data entered
        when object was initialized.

        """
        # make file in the same cwd, which is where sambvca will look
        with open("py2sambvca_input.inp", "w") as file:
            # write atoms to be deleted, if there are any
            if self.atoms_to_delete_ids is not None:
                file.writelines([
                    str(self.n_atoms_to_delete) + "\n",
                    str(self.atoms_to_delete_ids).replace(",","").replace("[", "").replace("]", "") + "\n"
                    ])
            else:
                file.write("0\n")
            # write user settings
            file.writelines([
                str(self.n_sphere_center_atoms) + "\n",
                str(self.sphere_center_atom_ids).replace(",", "").replace("[", "").replace("]", "") + "\n",
                str(self.n_z_atoms) + "\n",
                str(self.z_ax_atom_ids).replace(",", "").replace("[", "").replace("]", "") + "\n",
                str(self.n_xz_plane_atoms) + "\n",
                str(self.xz_plane_atoms_ids).replace(",", "").replace("[", "").replace("]", "") + "\n",
                str(self.sphere_radius) + "\n",
                str(self.displacement) + "\n",
                str(self.mesh_size) + "\n",
                str(self.remove_H) + "\n",
                str(self.orient_z) + "\n",
                str(self.write_surf_files) + "\n",
                "103\n"
                ])
            # write radii
            file.writelines(radii_table)
            # write the atom coordinates
            file.writelines(self.xyz_data)

    def calc(self):
        """
        Call SambVca based on the executable path given on initiliazation of py2sambvca.

        Be sure to write_input() before calling this function.

        """
        try:
            subprocess.run(
                [self.path_to_sambvcax, "py2sambvca_input"],
                stderr=subprocess.DEVNULL
                )
            return True
        except Exception as e:
            print(e)
            return False

    def get_buried_vol(self):
        """
        Retrieves the buried volume from a SambVca output file in the current working directory
        or False if it cannot find it.

        """
        # open the file, read data
        with open("py2sambvca_input.out", 'r') as file:
            file_data = file.readlines()
        pattern = re.compile(r"    The %V Bur of the molecule is:     (\d*\.\d*)")
        for line in file_data:
            m = pattern.search(line)
            if m:
                return m[1]
                
        return False

    def clean_files(self):
        """
        Remove all input and output files associated with py2sambvca.

        """
        [os.remove(i) for i in glob.glob("py2sambvca_input*")]


radii_table = [
    'H       1.28\n',
    'HE      1.64\n',
    'LI      2.13\n',
    'BE      1.79\n',
    'B       2.25\n',
    'C       1.99\n',
    'N       1.81\n',
    'O       1.78\n',
    'F       1.72\n',
    'NE      1.80\n',
    'NA      2.66\n',
    'MG      2.02\n',
    'AL      2.15\n',
    'SI      2.46\n',
    'P       2.11\n',
    'S       2.11\n',
    'CL      2.05\n',
    'AR      2.20\n',
    'K       3.22\n',
    'CA      2.70\n',
    'SC      2.52\n',
    'TI      2.47\n',
    'V       2.42\n',
    'CR      2.41\n',
    'MN      2.40\n',
    'FE      2.39\n',
    'CO      2.34\n',
    'NI      1.91\n',
    'CU      1.64\n',
    'ZN      1.63\n',
    'GA      2.19\n',
    'GE      2.47\n',
    'AS      2.16\n',
    'SE      2.22\n',
    'BR      2.16\n',
    'KR      2.36\n',
    'RB      3.55\n',
    'SR      2.91\n',
    'Y       2.71\n',
    'ZR      2.61\n',
    'NB      2.55\n',
    'MO      2.54\n',
    'TC      2.53\n',
    'RU      2.49\n',
    'RH      2.46\n',
    'PD      1.91\n',
    'AG      2.01\n',
    'CD      1.85\n',
    'IN      2.26\n',
    'SN      2.54\n',
    'SB      2.41\n',
    'TE      2.41\n',
    'I       2.32\n',
    'XE      2.53\n',
    'CS      4.01\n',
    'BA      3.14\n',
    'LA      2.84\n',
    'CE      2.83\n',
    'PR      2.81\n',
    'ND      2.80\n',
    'PM      2.78\n',
    'SM      2.76\n',
    'EU      2.75\n',
    'GD      2.74\n',
    'TB      2.73\n',
    'DY      2.70\n',
    'HO      2.69\n',
    'ER      2.68\n',
    'TM      2.66\n',
    'YB      2.64\n',
    'LU      2.62\n',
    'HF      2.61\n',
    'TA      2.60\n',
    'W       2.55\n',
    'RE      2.53\n',
    'OS      2.53\n',
    'IR      2.49\n',
    'PT      2.01\n',
    'AU      1.94\n',
    'HG      1.81\n',
    'TL      2.29\n',
    'PB      2.36\n',
    'BI      2.42\n',
    'PO      2.30\n',
    'AT      2.36\n',
    'RN      2.57\n',
    'FR      4.07\n',
    'RA      3.31\n',
    'AC      2.89\n',
    'TH      2.87\n',
    'PA      2.84\n',
    'U       2.18\n',
    'NP      2.80\n',
    'PU      2.84\n',
    'AM      2.85\n',
    'CM      2.87\n',
    'BK      2.85\n',
    'CF      2.87\n',
    'ES      2.87\n',
    'FM      2.87\n',
    'M       2.88\n',
    'NO      2.88\n',
    'LR      2.88\n'
]