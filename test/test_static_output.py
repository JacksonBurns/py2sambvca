import os
import unittest
import shutil

from py2sambvca import p2s


class Teststaticoutput(unittest.TestCase):
    """
    Test the various functionalities of py2sambvca on a static output file.
    """

    @classmethod
    def setUpClass(self):
        """Set input attirbutes for py2sambvca."""
        # current working directory
        cwd = os.getcwd()
        # directory for temporary test files
        self.working_dir = os.path.join(os.getcwd(), "test_dir")
        os.makedirs(self.working_dir, exist_ok=True)
        # get the saved data file, copy to working dir
        static_filepath = os.path.join(cwd, "test", "data", "static_output.out")
        self.static_output = os.path.join(self.working_dir, "py2sambvca_input.out")
        shutil.copy(static_filepath, self.static_output)
        # others needed to initialize a dummy py2sambvca
        self.xyz_file = os.path.join(cwd, "test", "data", "nhc.xyz")
        self.sphere_ids = [22]
        self.z_ids = [5]
        self.xz_ids = [1]

    def test_get_free_volume(self):
        """
        Retrieve the percent free volume from the saved output file.
        """
        test_p2s = p2s(
            self.xyz_file,
            self.sphere_ids,
            self.z_ids,
            self.xz_ids,
            working_dir=self.working_dir,
        )
        test_p2s.parse_output()
        result = test_p2s.get_percent_free_volume()
        self.assertEqual(result, 66.2)

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
        shutil.rmtree(self.working_dir)


if __name__ == "__main__":
    unittest.main()
