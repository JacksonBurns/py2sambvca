import os
import sys
import unittest

from py2sambvca import p2s


class Testpy2sambvca(unittest.TestCase):
    """
    Test the various functionalities of py2sambvca.
    """

    @classmethod
    def setUpClass(self):
        """Set input attirbutes for py2sambvca."""
        cwd = os.getcwd()
        self.xyz_file = os.path.join(cwd, "test", "data", "nhc.xyz")
        self.sphere_ids = [22]
        self.z_ids = [5]
        self.xz_ids = [1]
        self.delete_ids = [22]
        self.radius = 3.5
        self.displacement = 0.00
        self.mesh_size = 0.10
        self.remove_H = 1
        self.orient_Z = 0
        self.write_surf_files = 0
        if sys.platform == "win32":
            self.exe_path = os.path.abspath(
                os.path.join("executables", "sambvca21.exe")
            )
        else:
            self.exe_path = os.path.abspath(
                os.path.join("..", "executables", "sambvca21.x")
            )

    def test_no_init_error(self):
        """
        Call calc without first writing input, assure result is false.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.clean_files()
        self.assertFalse(test_p2s.calc())

    def test_full_init(self):
        """
        Initialization and attribute check.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            self.delete_ids,
            self.radius,
            self.displacement,
            self.mesh_size,
            self.remove_H,
            self.orient_Z,
            self.write_surf_files,
            self.exe_path,
        )
        self.assertEqual(
            self.sphere_ids,
            test_p2s.sphere_center_atom_ids,
            "sphere_center_atom_ids not set correctly.",
        )
        self.assertEqual(
            self.z_ids,
            test_p2s.z_ax_atom_ids,
            "z_ax_atom_ids not set correctly.",
        )
        self.assertEqual(
            self.xz_ids,
            test_p2s.xz_plane_atoms_ids,
            "xz_plane_atoms_ids not set correctly.",
        )
        self.assertEqual(
            self.delete_ids,
            test_p2s.atoms_to_delete_ids,
            "atoms_to_delete_ids not set correctly.",
        )
        self.assertEqual(
            self.radius,
            test_p2s.sphere_radius,
            "sphere_radius not set correctly.",
        )
        self.assertEqual(
            self.displacement,
            test_p2s.displacement,
            "displacement not set correctly.",
        )
        self.assertEqual(
            self.mesh_size,
            test_p2s.mesh_size,
            "mesh_size not set correctly.",
        )
        self.assertEqual(
            self.remove_H,
            test_p2s.remove_H,
            "remove_H not set correctly.",
        )
        self.assertEqual(
            self.orient_Z,
            test_p2s.orient_z,
            "orient_z not set correctly.",
        )
        self.assertEqual(
            self.write_surf_files,
            test_p2s.write_surf_files,
            "write_surf_files not set correctly.",
        )
        self.assertEqual(
            self.exe_path,
            test_p2s.path_to_sambvcax,
            "path_to_sambvcax not set correctly.",
        )

    def test_partial_init(self):
        """
        Provide only the bare minimum inputs, check all attributes.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
        )
        self.assertEqual(
            self.sphere_ids,
            test_p2s.sphere_center_atom_ids,
            "sphere_center_atom_ids not set correctly.",
        )
        self.assertEqual(
            self.z_ids,
            test_p2s.z_ax_atom_ids,
            "z_ax_atom_ids not set correctly.",
        )
        self.assertEqual(
            self.xz_ids,
            test_p2s.xz_plane_atoms_ids,
            "xz_plane_atoms_ids not set correctly.",
        )
        self.assertEqual(
            None,
            test_p2s.atoms_to_delete_ids,
            "atoms_to_delete_ids not set correctly.",
        )
        self.assertEqual(
            3.5,
            test_p2s.sphere_radius,
            "sphere_radius not set correctly.",
        )
        self.assertEqual(
            0.0,
            test_p2s.displacement,
            "displacement not set correctly.",
        )
        self.assertEqual(
            0.10,
            test_p2s.mesh_size,
            "mesh_size not set correctly.",
        )
        self.assertEqual(
            1,
            test_p2s.remove_H,
            "remove_H not set correctly.",
        )
        self.assertEqual(
            1,
            test_p2s.orient_z,
            "orient_z not set correctly.",
        )
        self.assertEqual(
            1,
            test_p2s.write_surf_files,
            "write_surf_files not set correctly.",
        )
        self.assertEqual(
            self.exe_path,
            test_p2s.path_to_sambvcax,
            "path_to_sambvcax not set correctly.",
        )

    def test_input_writer(self):
        """
        Test ability to write input file separately.
        Skip for the moment
        """
        return
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.write_input()
        with open("py2sambvca_input.inp") as file:
            self.assertIsNotNone(file)

    def test_calc(self):
        """
        Test call to calculator.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.write_input()
        self.assertTrue(test_p2s.calc())

    def test_buried_volume(self):
        """
        Test retrieval of buried volume.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.write_input()
        test_p2s.calc()
        self.assertEqual(test_p2s.get_buried_vol(), 55.7)

    def test_clean_files(self):
        """
        File cleaner should remove all py2sambvca files.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.write_input()
        test_p2s.clean_files()
        with self.assertRaises(FileNotFoundError):
            open("py2sambvca_input.inp")

    def test_run(self):
        """
        Full call to all functions using run.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        print(test_p2s.path_to_sambvcax)
        test_p2s.run()

    def test_getters(self):
        """
        Test getter methods.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            path_to_sambvcax=self.exe_path,
            verbose=0,
        )
        test_p2s.run()
        self.assertEqual(
            test_p2s.get_quadrant_free_volume()["SW"],
            21.0,
        )
        self.assertEqual(
            test_p2s.get_quadrant_buried_volume()["SW"],
            23.9,
        )
        self.assertEqual(
            test_p2s.get_quadrant_total_volume()["NW"],
            44.9,
        )
        self.assertEqual(
            test_p2s.get_quadrant_percent_buried_volume()["SE"],
            58.1,
        )
        self.assertEqual(
            test_p2s.get_quadrant_percent_free_volume()["SW"],
            53.2,
        )
        self.assertEqual(
            test_p2s.get_octant_free_volume()["SW-z"],
            16.7,
        )
        self.assertEqual(
            test_p2s.get_octant_buried_volume()["SW+z"],
            18.1,
        )
        self.assertEqual(
            test_p2s.get_octant_total_volume()["SE-z"],
            22.4,
        )
        self.assertEqual(
            test_p2s.get_octant_percent_buried_volume()["SE+z"],
            90.2,
        )
        self.assertEqual(
            test_p2s.get_octant_percent_free_volume()["NW+z"],
            90.4,
        )
        self.assertEqual(
            test_p2s.get_free_volume(),
            79.4,
        )
        self.assertEqual(
            test_p2s.get_buried_volume(),
            100.0,
        )
        self.assertEqual(
            test_p2s.get_exact_volume(),
            179.6,
        )
        self.assertEqual(
            test_p2s.get_total_volume(),
            179.4,
        )
        self.assertEqual(
            test_p2s.get_percent_buried_volume(),
            55.7,
        )
        self.assertEqual(
            test_p2s.get_percent_free_volume(),
            44.3,
        )
        self.assertEqual(
            test_p2s.get_percent_total_volume(),
            99.9,
        )

    @classmethod
    def tearDownClass(self):
        """Clean up the mess."""
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
        )
        test_p2s.clean_files()


if __name__ == "__main__":
    unittest.main()
