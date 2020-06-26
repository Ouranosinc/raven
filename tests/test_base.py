from .common import TESTDATA
import raven
from raven.models import Raven, Ostrich
from raven.models.base import get_diff_level
import tempfile
import numpy as np
from pathlib import Path
import pytest


has_singularity = raven.raven_simg.exists()


class TestRaven:
    def test_gr4j(self):
        rvs = TESTDATA["raven-gr4j-cemaneige-nc-rv"]
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]

        model = Raven()
        model.configure(rvs)
        model(ts)

    def test_mohyse(self):
        rvs = TESTDATA["raven-mohyse-rv"]
        ts = list(TESTDATA["raven-mohyse-ts"])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model(ts)

    def test_hmets(self):
        rvs = TESTDATA["raven-hmets-rv"]
        ts = TESTDATA["raven-hmets-ts"]

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model(ts)

    def test_hbvec(self):
        rvs = TESTDATA["raven-hbv-ec-rv"]
        ts = TESTDATA["raven-hbv-ec-ts"]

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model(ts)

    @pytest.mark.skipif(not has_singularity, reason="Singularity is not available.")
    def test_singularity(self):
        rvs = TESTDATA["raven-gr4j-cemaneige-nc-rv"]
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]

        model = Raven()
        model.singularity = True
        model.configure(rvs)
        model.run(ts)


class TestOstrich:
    def test_gr4j_with_no_tags(self):
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        ost = TESTDATA["ostrich-gr4j-cemaneige-rv"]

        model = Ostrich()
        model.configure(ost)
        model(ts)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        assert len(opt_para) == 6
        assert isinstance(opt_func, float)

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(
            opt_para,
            [2.424726, 3.758972, 204.3856, 5.866946, 16.60408, 0.3728098],
            3,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            -0.50717,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-gr4j-with-cema-neige
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -0.5779910, 4,
        #                                err_msg='calibrated NSE is not matching expected value')

        assert Path(model.outputs["calibration"]).exists()

    def test_mohyse_with_no_tags(self):
        ts = TESTDATA["raven-mohyse-nc-ts"]
        ost = TESTDATA["ostrich-mohyse-rv"]

        model = Ostrich()
        model.configure(ost)
        model(ts)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        assert len(opt_para) == 10
        assert isinstance(opt_func, float)

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(
            opt_para,
            [
                7.721801e00,
                8.551484e-01,
                1.774571e01,
                1.627677e00,
                7.702450e-02,
                9.409600e-01,
                6.941596e-01,
                8.207870e-01,
                8.154455e00,
                1.018226e01,
            ],
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            -0.3826810,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-mohyse
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [1.517286E+01, 7.112556E-01, 1.981243E+01, -4.193046E+00,
        #                                           1.791486E-01, 9.774897E-01, 5.353541E-01, 6.686806E-01,
        #                                           1.040908E+01, 1.132304E+01, 8.831552E-02], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -0.3857010, 4,
        #                                err_msg='calibrated NSE is not matching expected value')

        assert Path(model.outputs["calibration"]).exists()

    def test_hmets_with_no_tags(self):
        ts = TESTDATA["raven-hmets-nc-ts"]
        ost = TESTDATA["ostrich-hmets-rv"]

        model = Ostrich()
        model.configure(ost)
        model(ts)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        assert len(opt_para) == 21
        assert isinstance(opt_func, float)

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #

        expected_value = [
            1.777842e01,
            3.317211e00,
            5.727342e00,
            1.419491e00,
            1.382141e01,
            1.637954e01,
            7.166296e-01,
            1.389346e-01,
            2.620464e-02,
            2.245525e-01,
            2.839426e-02,
            -2.003810e00,
            9.479623e-01,
            4.803857e-01,
            2.524914e00,
            4.117232e-01,
            1.950058e-02,
            4.494123e-02,
            1.405815e-03,
            2.815803e-02,
            1.007823e00,
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            1.43474,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-hmets
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.008045E+00, 7.960246E-02, 4.332698E+00, 4.978125E-01,
        #                                           1.997029E+00, 6.269773E-01, 1.516961E+00, 8.180383E-02,
        #                                           6.730663E-02, 2.137822E-02, 2.097163E-02, 1.773348E+00,
        #                                           3.036039E-01, 1.928524E-02, 1.758471E+00, 8.942299E-01,
        #                                           8.741980E-03, 5.036474E-02, 9.465804E-03, 1.851839E-01,
        #                                           1.653934E-01, 2.624006E+00, 8.868485E-02, 9.259195E+01,
        #                                           8.269670E+01], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.350490E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')

        assert Path(model.outputs["calibration"]).exists()

    def test_hbvec_with_no_tags(self):
        ts = TESTDATA["raven-hbv-ec-nc-ts"]
        ost = TESTDATA["ostrich-hbv-ec-rv"]

        model = Ostrich()
        model.configure(ost)
        model(ts)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        assert len(opt_para) == 21
        assert isinstance(opt_func, float)

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(
            opt_para,
            [
                -8.317931e-01,
                4.072232e00,
                2.001574e00,
                5.736299e-03,
                9.985144e-02,
                4.422529e-01,
                3.438486e00,
                8.055843e01,
                4.440133e-01,
                8.451082e-02,
                2.814201e00,
                7.327970e-01,
                1.119773e00,
                1.161223e-03,
                4.597179e-01,
                1.545857e01,
                1.223865e00,
                4.452843e-01,
                9.492006e-01,
                9.948123e-01,
                1.110682e00,
            ],
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            0.242187,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-environment-
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.984519E-02, 4.072232E+00, 2.001574E+00, 3.473693E-02,
        #                                           9.985144E-02, 5.060520E-01, 2.944343E+00, 3.832455E+01,
        #                                           4.606565E-01, 6.303738E-02, 2.277781E+00, 4.873686E+00,
        #                                           5.718813E-01, 4.505643E-02, 8.776511E-01, 1.894145E+01,
        #                                           2.036937E+00, 4.452843E-01, 6.771759E-01, 1.206053E+00,
        #                                           1.024278E+00], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.034670E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')

        assert Path(model.outputs["calibration"]).exists()


def test_get_diff_level():
    fn = Path("/") / "a" / "b" / "c.txt"
    files = [fn, Path("/") / "a" / "b" / "d.txt"]
    assert get_diff_level(files) == 3
    assert fn.relative_to(Path(*fn.parts[:3])) == Path("c.txt")

    files = [fn, Path("/") / "a" / "b1" / "c.txt"]
    assert get_diff_level(files) == 2
    assert fn.relative_to(Path(*fn.parts[:2])) == Path("b/c.txt")

    files = [fn, Path("/") / "a" / "b1" / "b2" / "c.txt"]
    assert get_diff_level(files) == 2
    assert files[0].relative_to(Path(*fn.parts[:2])) == Path("b/c.txt")
    assert files[1].relative_to(Path(*files[1].parts[:2])) == Path("b1/b2/c.txt")
