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
        """Set input attirbutes for py2sambvca.
        """
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
            self.exe_path = os.path.join(cwd, "sambvca21.exe")
        else:
            self.exe_path = os.path.join(cwd, "sambvca21.x")

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
        )
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
            self.xyz_file,
            test_p2s.xyz_data,
            "xyz coordinates not set correctly.",
        )
        # repeat above for each attribute

    def test_partial_init(self):
        """
        Provide only the bare minimum inputs, check all attributes.
        """
        pass

    def test_input_writer(self):
        pass

    def test_calc(self):
        pass

    def test_buried_volume(self):
        pass

    def test_parse_output(self):
        pass

    def test_clean_files(self):
        pass

        # test_p2s.write_input()
        # test_p2s.calc()
        # test_p2s.get_buried_vol()
        # test_p2s.parse_output()
        # test_p2s.clean_files()


if __name__ == '__main__':
    unittest.main()
