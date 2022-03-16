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

    def test_py2sambvca(self):
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
        test_p2s.write_input()
        test_p2s.calc()
        test_p2s.get_buried_vol()
        test_p2s.parse_output()
        test_p2s.clean_files()


if __name__ == '__main__':
    unittest.main()
