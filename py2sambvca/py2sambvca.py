import subprocess
import re
import glob
import os
import uuid
from typing import List, Union

from py2sambvca.radii_tables import table_lookup, format_radii_table


class py2sambvca:
    """
    Wrapper class for py2sambvca functions.

    Call this class to instantiate a py2sambvca object, which has methods to write input, call SambVca,
    and retrieve output.

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
    working_dir (path): Path to the working directory where the output and input files are generated (default py2sambvca_temp_UUID)
    verbose (int): 0 for no output, 1 for some output, 2 for the most output
    radii_table (dict or str): atomic symbols (ie. H, LI) and radii (Angstrom), "default" (bondii radii*1.17), or "vdw" for van Der Waals
    """

    def __init__(
        self,
        xyz_filepath: str,
        sphere_center_atom_ids: List[int],
        z_ax_atom_ids: List[int],
        xz_plane_atoms_ids: List[int],
        atoms_to_delete_ids: List[int] = None,
        sphere_radius: float = 3.5,
        displacement: float = 0.0,
        mesh_size: float = 0.10,
        remove_H: int = 1,
        orient_z: int = 1,
        write_surf_files: int = 1,
        path_to_sambvcax: str = "sambvca.exe",
        working_dir: str = "py2sambvca_temp_UUID",
        verbose: int = 1,
        radii_table: Union[str, dict] = "default",
    ):
        """
        Wrapper class for py2sambvca functions.

        Call this class to instantiate a py2sambvca object, which has methods to write input, call SambVca,
        and retrieve output.

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
        working_dir (path): Path to the working directory where the output and input files are generated (default py2sambvca_temp_UUID)
        verbose (int): 0 for no output, 1 for some output, 2 for the most output
        radii_table (dict or str): atomic symbols (ie. H, LI) and radii (Angstrom), "default" (bondii radii*1.17), or "vdw" for van Der Waals
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
        if xyz_filepath.endswith(".xyz"):
            with open(xyz_filepath, "r") as file:
                self.xyz_data = file.readlines()
        else:
            raise RuntimeError(f"Invalid xyz_filepath ({xyz_filepath})")

        # assign the path to the calculator
        self.path_to_sambvcax = path_to_sambvcax

        # assign working directory path
        if working_dir == "py2sambvca_temp_UUID":
            working_dir = working_dir.replace("UUID", uuid.uuid4().hex)
        self.working_dir = working_dir
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # integer verbosity
        self.verbose = verbose

        # make results accesible from object directly
        self.total_results = None
        self.quadrant_results = None
        self.octant_results = None

        # data table
        self.__radii_table = format_radii_table(
            table_lookup[radii_table] if type(radii_table) is str else radii_table
        )

    def write_input(self):
        """
        Write input for the Sambvca buried-volume Fortran calculator based on the data entered
        when object was initialized.

        """
        # make file in the same cwd, which is where sambvca will look
        with open(os.path.join(self.working_dir, "py2sambvca_input.inp"), "w") as file:
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
            file.writelines(self.__radii_table)
            # write the atom coordinates
            file.writelines(self.xyz_data)

    def calc(self):
        """
        Call SambVca based on the executable path given on initialization of py2sambvca.

        Be sure to write_input() before calling this function.

        """
        if not os.path.exists(self.path_to_sambvcax):
            raise RuntimeError(
                "sambvca executable not found at provided path ({:s})".format(
                    self.path_to_sambvcax
                )
            )
        try:
            result = subprocess.run(
                [
                    self.path_to_sambvcax,
                    os.path.join(self.working_dir, "py2sambvca_input"),
                ],
                stderr=subprocess.DEVNULL,
            )
            result.check_returncode()
            return True
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(e)
            return False

    def clean_files(self):
        """
        Remove all input and output files associated with py2sambvca.

        """
        for f in glob.glob(os.path.join(self.working_dir, "py2sambvca_input*")):
            os.remove(f)

    def parse_output(self):
        """Parse output file for total, quadrant, and octant results.

        Returns:
            total_results (dict): Results for total
            quadrant_results (dict): Quadrant-decomposed results
            octant_results (dict): Octant-decomposed results
        """
        try:
            return self._parse_output()
        except TypeError as te:
            raise RuntimeError(
                "Unable to retrieve output from sambvca21 output file ({:s}). "
                "Did sambvca21 termiante normally?".format(
                    os.path.join(self.working_dir, "py2sambvca_input.out")
                )
            ) from te

    def _parse_output(self):
        """Helper funciton for parse output to better wrap with try/except."""
        # total results
        m1 = self.get_regex(r"^\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)$")

        m2 = self.get_regex(r"^\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)$")

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
        try:
            with open(
                os.path.join(self.working_dir, "py2sambvca_input.out"), "r"
            ) as file:
                file_data = file.readlines()
        except FileNotFoundError:
            raise RuntimeError(
                "Results not yet retrieved ({:s} not found). "
                "Call p2s.run() or p2s.calc() before using this function.".format(
                    os.path.join(self.working_dir, "py2sambvca_input.out")
                )
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
                "Results not yet retrieved ({} not found). "
                "Call p2s.run() or p2s.parse_output() before using this function.".format(
                    os.path.join(self.working_dir, "py2sambvca_input.out")
                )
            )
        try:
            if octant and quadrant:
                raise RuntimeError("Specify either quadrant or octant, not both")
            elif quadrant:
                return self.quadrant_results[key]
            elif octant:
                return self.octant_results[key]
            else:
                return self.total_results[key]
        except KeyError as e:
            if self.verbose:
                print(e)
            raise RuntimeError(f'Invalid parameter name "{key}"')

    def run(self):
        self.write_input()
        self.calc()
        self.parse_output()
        self.clean_files()
        return self.total_results, self.quadrant_results, self.octant_results
