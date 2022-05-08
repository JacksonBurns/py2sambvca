import subprocess
import re
import glob
import os
from tempfile import mkdtemp
import shutil
import sys
import numpy as np
import pandas as pd
from collections import defaultdict


class py2sambvca:
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
    sphere_radius (float): Sphere radius in Angstrom (default 3.5). Only used for single point calculation.
    displacement (float): Displacement of oriented molecule from sphere center in Angstrom (default 0.0)
    mesh_size (float): Mesh size for numerical integration (default 0.10)
    remove_H (int): 0/1 Do not remove/remove H atoms from Vbur calculation (default 1)
    orient_z (int): 0/1 Molecule oriented along negative/positive Z-axis (default 1)
    write_surf_files (int): 0/1 Do not write/write files for top and bottom surfaces (default 1)
    path_to_sambvcax (str): Path to the SambVca executable. Only needed to use py2sambvca.calc()( default "/path/to/executable/sambvca.x")
    verbose (int): 0 for no output, 1 for some output, 2 for the most output
    """

    def __init__(
        self,
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
        path_to_sambvcax=None,
        verbose=1,
    ):
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
        verbose (int): 0 for no output, 1 for some output, 2 for the most output
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
        self.remove_H = int(remove_H)
        self.orient_z = int(orient_z)
        self.write_surf_files = int(write_surf_files)
        # open the xyz file, read the data
        with open(xyz_filepath, "r") as file:
            self.xyz_data = file.readlines()

        # assign the path to the calculator
        if path_to_sambvcax is None:
            if sys.platform == "win32":
                self.path_to_sambvcax = os.path.join(
                    "..", "executables", "sambvca21.exe"
                )
            else:
                self.path_to_sambvcax = os.path.join("..", "executables", "sambvca21.x")
        else:
            self.path_to_sambvcax = path_to_sambvcax
        self.verbose = verbose

        # make results accesible from object directly
        self.total_results = None
        self.quadrant_results = None
        self.octant_results = None

        # make a temporary folder to work in.
        self.tmp_dir = mkdtemp()
        self.input_file = os.path.join(self.tmp_dir, "py2sambvca_input.inp")

    def write_input(self):
        """
        Write input for the Sambvca buried-volume Fortran calculator based on the data entered
        when object was initialized.

        """
        # make file in the same cwd, which is where sambvca will look
        with open(self.input_file, "w") as file:
            # write atoms to be deleted, if there are any
            if self.atoms_to_delete_ids is not None:
                file.writelines(
                    [
                        str(self.n_atoms_to_delete) + "\n",
                        str(self.atoms_to_delete_ids)
                        .replace(",", "")
                        .replace("[", "")
                        .replace("]", "")
                        + "\n",
                    ]
                )
            else:
                file.write("0\n")
            # write user settings
            file.writelines(
                [
                    str(self.n_sphere_center_atoms) + "\n",
                    str(self.sphere_center_atom_ids)
                    .replace(",", "")
                    .replace("[", "")
                    .replace("]", "")
                    + "\n",
                    str(self.n_z_atoms) + "\n",
                    str(self.z_ax_atom_ids)
                    .replace(",", "")
                    .replace("[", "")
                    .replace("]", "")
                    + "\n",
                    str(self.n_xz_plane_atoms) + "\n",
                    str(self.xz_plane_atoms_ids)
                    .replace(",", "")
                    .replace("[", "")
                    .replace("]", "")
                    + "\n",
                    str(self.sphere_radius) + "\n",
                    str(self.displacement) + "\n",
                    str(self.mesh_size) + "\n",
                    str(self.remove_H) + "\n",
                    str(self.orient_z) + "\n",
                    str(self.write_surf_files) + "\n",
                    "103\n",
                ]
            )
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
            result = subprocess.run(
                [self.path_to_sambvcax, os.path.splitext(self.input_file)[0]],
                stderr=subprocess.DEVNULL,
            )
            result.check_returncode()
            return True
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(e)
            return False

    def get_buried_vol(self):
        """
        Retrieves the buried volume from a SambVca output file in the current working directory
        or False if it cannot find it.
        """
        m = self.get_regex(r"^[ ]{4}The %V Bur of the molecule is:[ ]{4,5}(\d*\.\d*)$")
        return float(m[1])

    def clean_files(self):
        """
        Remove all input and output files associated with py2sambvca.
        """
        shutil.rmtree(self.tmp_dir)

    def parse_output(self):
        """Parse output file for total, quandrant, and octant results.

        Returns:
            total_results (dict): Results for total
            quadrant_results (dict): Quadrant-decomposed results
            octant_results (dict): Octant-decomposed results
        """
        # total results
        m1 = self.get_regex(
            r"^[ ]{5,6}(\d*\.\d*)[ ]{5,6}(\d*\.\d*)[ ]{5,6}(\d*\.\d*)[ ]{5,6}(\d*\.\d*)$"
        )

        m2 = self.get_regex(r"^[ ]{5,6}(\d*\.\d*)[ ]{5,6}(\d*\.\d*)[ ]{5,6}(\d*\.\d*)$")

        if m1 is None:
            print(
                f"The calculation for {self.sphere_radius} could not produce a suitable result. Most likely this is due to the chosen radius beeing too large or too small."
            )
            m1 = [0, 0, 0, 0, 0, 0]
        if m2 is None:
            m2 = [0, 0, 0, 0, 0, 0]

        total_results = {
            "free_volume": float(m1[1]),
            "buried_volume": float(m1[2]),
            "total_volume": float(m1[3]),
            "exact_volume": float(m1[4]),
            "percent_buried_volume": float(m2[2]),
            "percent_free_volume": float(m2[1]),
            "percent_total_volume": float(m2[3]),
        }

        # quadrant and octant results
        quadrant_regions = ["SW", "NW", "NE", "SE"]
        quadrant_results = {
            "free_volume": {},
            "buried_volume": {},
            "total_volume": {},
            "percent_free_volume": {},
            "percent_buried_volume": {},
        }
        octant_regions = [
            r"SW\-z",
            r"NW\-z",
            r"NE\-z",
            r"SE\-z",
            r"SW\+z",
            r"NW\+z",
            r"NE\+z",
            r"SE\+z",
        ]
        octant_results = quadrant_results.copy()

        for region, result_dict in zip(
            [quadrant_regions, octant_regions],
            [quadrant_results, octant_results],
        ):
            for r in region:
                m = self.get_regex(
                    r"^ "
                    + r
                    + r"\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)$"
                )
                result_dict["free_volume"][r.replace("\\", "")] = float(m[1])
                result_dict["buried_volume"][r.replace("\\", "")] = float(m[2])
                result_dict["total_volume"][r.replace("\\", "")] = float(m[3])
                result_dict["percent_free_volume"][r.replace("\\", "")] = float(m[4])
                result_dict["percent_buried_volume"][r.replace("\\", "")] = float(m[5])

        self.total_results = total_results
        self.quadrant_results = quadrant_results
        self.octant_results = octant_results
        return total_results, quadrant_results, octant_results

    def get_quadrant_result(self, key):
        """Get a result for the quadrants.

        Args:
            key (str): type of result

        Returns:
            dict: values requested per quadrant
        """
        return self.get(key, quadrant=True)

    def get_quadrant_free_volume(self):
        """Get the quadrant free volume.

        Returns:
            dict: free volume per quadrant
        """
        return self.get_quadrant_result("free_volume")

    def get_quadrant_buried_volume(self):
        """Get the quadrant buried volume.

        Returns:
            dict: buried volume per quadrant
        """
        return self.get_quadrant_result("buried_volume")

    def get_quadrant_total_volume(self):
        """Get the quadrant total volume.

        Returns:
            dict: total volume per quadrant
        """
        return self.get_quadrant_result("total_volume")

    def get_quadrant_percent_buried_volume(self):
        """Get the quadrant percent buried volume

        Returns:
            dict: total percent buried volume per quadrant
        """
        return self.get_quadrant_result("percent_buried_volume")

    def get_quadrant_percent_free_volume(self):
        """Get the total percent buried volume

        Returns:
            dict: total percent buried volume per quadrant
        """
        return self.get_quadrant_result("percent_buried_volume")

    def get_octant_result(self, key):
        """Get a result for the octant.

        Args:
            key (str): type of result

        Returns:
            dict: values requested per octant
        """
        return self.get(key, octant=True)

    def get_octant_free_volume(self):
        """Get the octant free volume.

        Returns:
            dict: free volume per octant
        """
        return self.get_octant_result("free_volume")

    def get_octant_buried_volume(self):
        """Get the octant buried volume.

        Returns:
            dict: buried volume per octant
        """
        return self.get_octant_result("buried_volume")

    def get_octant_total_volume(self):
        """Get the octant total volume.

        Returns:
            dict: total volume per octant
        """
        return self.get_octant_result("total_volume")

    def get_octant_percent_buried_volume(self):
        """Get the octant percent buried volume

        Returns:
            dict: total percent buried volume per octant
        """
        return self.get_octant_result("percent_buried_volume")

    def get_octant_percent_free_volume(self):
        """Get the octant percent buried volume

        Returns:
            dict: total percent buried volumes per octant
        """
        return self.get_octant_result("percent_buried_volume")

    def get_free_volume(self):
        """Get the free volume.

        Returns:
            float: free volume
        """
        return self.get("free_volume")

    def get_buried_volume(self):
        """Get the buried volume.

        Returns:
            float: buried volume
        """
        return self.get("buried_volume")

    def get_exact_volume(self):
        """Get the exact volume.

        Returns:
            float: exact volume
        """
        return self.get("exact_volume")

    def get_total_volume(self):
        """Get the total volume.

        Returns:
            float: total volume
        """
        return self.get("total_volume")

    def get_percent_buried_volume(self):
        """Get the total percent buried volume

        Returns:
            float: total percent buried volume
        """
        return self.get("percent_buried_volume")

    def get_percent_free_volume(self):
        """Get the total percent buried volume

        Returns:
            float: total percent buried volume
        """
        return self.get("percent_free_volume")

    def get_percent_total_volume(self):
        """Get the percent total volume

        Returns:
            float: percent total volume
        """
        return self.get("percent_total_volume")

    def get_regex(self, regex):
        """Open the output file and search for a line matching a regex pattern.

        Args:
            regex (str): regex to search
        """
        output_file = os.path.join(os.path.splitext(self.input_file)[0] + ".out")
        try:
            with open(output_file, "r") as file:
                file_data = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(
                f"""
                Results not yet retrieved ({output_file} not found).
                Call p2s.run() or p2s.parse_output() before using this function.
                """
            )
        pattern = re.compile(regex)
        for line in file_data:
            m = pattern.search(line)
            if m:
                return m

    def get(self, key, quadrant=False, octant=False):
        """
        Accept a key in the output of parse results, return it.
        """
        if self.total_results is None:
            raise RuntimeError(
                """
                Results not yet retrieved (py2sambvca_input.out not found).
                Call p2s.run() or p2s.parse_output() before using this function.
                """
            )
        if not (octant or quadrant):
            return self.total_results[key]
        elif quadrant:
            return self.quadrant_results[key]
        elif octant:
            return self.octant_results[key]
        else:
            return False

    def run(self):
        self.write_input()
        self.calc()
        self.parse_output()
        self.clean_files()
        return self.total_results, self.quadrant_results, self.octant_results

    def run_range(self, r_min, r_max, nstep=50):
        """
        This function is designed to scan a range of sphere_radii.
        The results are stores in a pandas dataframe.
        Note: If the given radius is too big or too small no output will be stored for this radius.
        Args:
            r_min (number): minimum radius
            r_max (number): maximum radius
            nstep (int, optional): Number of steps. Defaults to 50.

        Returns:
            pandas.DataFrame: A DateFrame object for easy accsess to the results.
        """

        # store original input_file and radius values:
        old_radius = self.sphere_radius
        old_input_filee = self.input_file

        dict_total_results = defaultdict(list)

        for r in np.linspace(r_min, r_max, nstep):
            self.input_file = os.path.join(self.tmp_dir, f"py2sambvca_input_r_{r}.inp")
            self.sphere_radius = r
            self.write_input()
            self.calc()
            total_results, quadrant_results, octant_results = self.parse_output()
            if total_results["total_volume"] != 0:
                dict_total_results["r"].append(r)
                [
                    dict_total_results[key].append(value)
                    for key, value in total_results.items()
                ]

        df_total_results = pd.DataFrame(dict_total_results)
        # return original values
        self.sphere_radius = old_radius
        self.input_file = old_input_filee

        return df_total_results


radii_table = [
    "H       1.28\n",
    "HE      1.64\n",
    "LI      2.13\n",
    "BE      1.79\n",
    "B       2.25\n",
    "C       1.99\n",
    "N       1.81\n",
    "O       1.78\n",
    "F       1.72\n",
    "NE      1.80\n",
    "NA      2.66\n",
    "MG      2.02\n",
    "AL      2.15\n",
    "SI      2.46\n",
    "P       2.11\n",
    "S       2.11\n",
    "CL      2.05\n",
    "AR      2.20\n",
    "K       3.22\n",
    "CA      2.70\n",
    "SC      2.52\n",
    "TI      2.47\n",
    "V       2.42\n",
    "CR      2.41\n",
    "MN      2.40\n",
    "FE      2.39\n",
    "CO      2.34\n",
    "NI      1.91\n",
    "CU      1.64\n",
    "ZN      1.63\n",
    "GA      2.19\n",
    "GE      2.47\n",
    "AS      2.16\n",
    "SE      2.22\n",
    "BR      2.16\n",
    "KR      2.36\n",
    "RB      3.55\n",
    "SR      2.91\n",
    "Y       2.71\n",
    "ZR      2.61\n",
    "NB      2.55\n",
    "MO      2.54\n",
    "TC      2.53\n",
    "RU      2.49\n",
    "RH      2.46\n",
    "PD      1.91\n",
    "AG      2.01\n",
    "CD      1.85\n",
    "IN      2.26\n",
    "SN      2.54\n",
    "SB      2.41\n",
    "TE      2.41\n",
    "I       2.32\n",
    "XE      2.53\n",
    "CS      4.01\n",
    "BA      3.14\n",
    "LA      2.84\n",
    "CE      2.83\n",
    "PR      2.81\n",
    "ND      2.80\n",
    "PM      2.78\n",
    "SM      2.76\n",
    "EU      2.75\n",
    "GD      2.74\n",
    "TB      2.73\n",
    "DY      2.70\n",
    "HO      2.69\n",
    "ER      2.68\n",
    "TM      2.66\n",
    "YB      2.64\n",
    "LU      2.62\n",
    "HF      2.61\n",
    "TA      2.60\n",
    "W       2.55\n",
    "RE      2.53\n",
    "OS      2.53\n",
    "IR      2.49\n",
    "PT      2.01\n",
    "AU      1.94\n",
    "HG      1.81\n",
    "TL      2.29\n",
    "PB      2.36\n",
    "BI      2.42\n",
    "PO      2.30\n",
    "AT      2.36\n",
    "RN      2.57\n",
    "FR      4.07\n",
    "RA      3.31\n",
    "AC      2.89\n",
    "TH      2.87\n",
    "PA      2.84\n",
    "U       2.18\n",
    "NP      2.80\n",
    "PU      2.84\n",
    "AM      2.85\n",
    "CM      2.87\n",
    "BK      2.85\n",
    "CF      2.87\n",
    "ES      2.87\n",
    "FM      2.87\n",
    "M       2.88\n",
    "NO      2.88\n",
    "LR      2.88\n",
]
