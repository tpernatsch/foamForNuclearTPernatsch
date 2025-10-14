import numpy as np
import pyvista as pv

from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import NonUniformList, addParameter
from foamForNuclear.openfoamFile import OpenFOAMFile


class QuadratureSet(OpenFOAMFile):
    """
    Quadrature set object used in discrete ordinates solver.

    Parameters
    ----------
    region : str
        Name of the region
    discreteDirectionsOct : int
        Number of directions
    directionUxOct : NonUniformList | list
        List of X component per direction (default `None`).
    directionUyOct : NonUniformList | list
        List of Y component per direction (default `None`).
    directionUzOct : NonUniformList | list
        List of Z component per direction (default `None`).
    directionWeightsOct : NonUniformList | list
        List of weight per direction (default `None`).
    default : int | str {1, 4, 16, 'S4', 'S8'}
        Apply default quadrature set
    """

    def __init__(
            self,
            region: str,
            discreteDirectionsOct: int=None,
            directionUxOct: NonUniformList | list=None,
            directionUyOct: NonUniformList | list=None,
            directionUzOct: NonUniformList | list=None,
            directionWeightsOct: NonUniformList | list=None,
            default: int | str=None
        ):
        super().__init__("quadratureSet", folder="constant", region=region)

        self.discreteDirectionsOct = discreteDirectionsOct
        self.directionUxOct = directionUxOct
        self.directionUyOct = directionUyOct
        self.directionUzOct = directionUzOct
        self.directionWeightsOct = directionWeightsOct

        if (default is not None):
            check_value("default", default, {1, 4, 16, 'S4', 'S8'})

        if (default == 1 or default == '1'):
            self.set_default_1direction()
        elif (default == 4 or default == '4'):
            self.set_default_4direction()
        elif (default == 16 or default == '16'):
            self.set_default_16direction()
        elif (default == 'S4'):
            self.set_default_S4()
        elif (default == 'S8'):
            self.set_default_S8()


    @property
    def discreteDirectionsOct(self):
        return self._discreteDirectionsOct

    @discreteDirectionsOct.setter
    def discreteDirectionsOct(self, discreteDirectionsOct):
        check_type("discreteDirectionsOct", discreteDirectionsOct, int, none_ok=True)
        if (discreteDirectionsOct is not None):
            check_positive("discreteDirectionsOct", discreteDirectionsOct, is_strict=True)
        self._discreteDirectionsOct = discreteDirectionsOct


    @property
    def directionUxOct(self):
        return self._directionUxOct

    @directionUxOct.setter
    def directionUxOct(self, directionUxOct):
        check_type("directionUxOct", directionUxOct, (list, NonUniformList), none_ok=True)
        if (isinstance(directionUxOct, list)):
            self._directionUxOct = NonUniformList(directionUxOct)
        else:
            self._directionUxOct = directionUxOct


    @property
    def directionUyOct(self):
        return self._directionUyOct

    @directionUyOct.setter
    def directionUyOct(self, directionUyOct):
        check_type("directionUyOct", directionUyOct, (list, NonUniformList), none_ok=True)
        if (isinstance(directionUyOct, list)):
            self._directionUyOct = NonUniformList(directionUyOct)
        else:
            self._directionUyOct = directionUyOct


    @property
    def directionUzOct(self):
        return self._directionUzOct

    @directionUzOct.setter
    def directionUzOct(self, directionUzOct):
        check_type("directionUzOct", directionUzOct, (list, NonUniformList), none_ok=True)
        if (isinstance(directionUzOct, list)):
            self._directionUzOct = NonUniformList(directionUzOct)
        else:
            self._directionUzOct = directionUzOct


    @property
    def directionWeightsOct(self):
        return self._directionWeightsOct

    @directionWeightsOct.setter
    def directionWeightsOct(self, directionWeightsOct):
        check_type("directionWeightsOct", directionWeightsOct, (list, NonUniformList), none_ok=True)
        if (isinstance(directionWeightsOct, list)):
            self._directionWeightsOct = NonUniformList(directionWeightsOct)
        else:
            self._directionWeightsOct = directionWeightsOct


    def set_default_1direction(self):
        """
        From Godiva_SN tutorial quadratueSet1
        """
        self.discreteDirectionsOct = 1
        self.directionUxOct = [1]
        self.directionUyOct = [1]
        self.directionUzOct = [1]
        self.directionWeightsOct = [1]


    def set_default_4direction(self):
        """
        From Godiva_SN tutorial quadratueSet4
        """
        self.discreteDirectionsOct = 4
        self.directionUxOct = [
            2.09769106510711E-01,
            9.54983687770317E-01,
            2.09769106510711E-01,
            5.77350269189625E-01
        ]
        self.directionUyOct = [
            2.09769106510711E-01,
            2.09769106510711E-01,
            9.54983687770317E-01,
            5.77350269189625E-01
        ]
        self.directionUzOct = [
            9.54983687770317E-01,
            2.09769106510711E-01,
            2.09769106510711E-01,
            5.77350269189625E-01
        ]
        self.directionWeightsOct = [
            3.39836909454130E-01,
            3.39836909454103E-01,
            3.39836909454103E-01,
            5.51285598432481E-01
        ]


    def set_default_16direction(self):
        """
        From Godiva_SN tutorial quadratueSet16
        """
        self.discreteDirectionsOct = 16
        self.directionUxOct = [
            1.22785291159728E-01,
            5.31419325021509E-01,
            1.03129501638901E-01,
            2.35702260395515E-01,
            8.40807829938193E-01,
            9.84808379609774E-01,
            8.40807829938210E-01,
            9.42809041582063E-01,
            1.03129501638887E-01,
            5.31419325021508E-01,
            1.22785291159756E-01,
            2.35702260395515E-01,
            6.86947072008541E-01,
            2.37081084268253E-01,
            6.86947072008531E-01,
            5.77350269189625E-01
        ]
        self.directionUyOct = [
            1.22785291159728E-01,
            1.03129501638898E-01,
            5.31419325021504E-01,
            2.35702260395515E-01,
            1.03129501638887E-01,
            1.22785291159753E-01,
            5.31419325021502E-01,
            2.35702260395515E-01,
            8.40807829938193E-01,
            8.40807829938206E-01,
            9.84808379609773E-01,
            9.42809041582063E-01,
            6.86947072008541E-01,
            6.86947072008531E-01,
            2.37081084268254E-01,
            5.77350269189625E-01
        ]
        self.directionUzOct = [
            9.84808379609780E-01,
            8.40807829938206E-01,
            8.40807829938209E-01,
            9.42809041582063E-01,
            5.31419325021532E-01,
            1.22785291159753E-01,
            1.03129501638902E-01,
            2.35702260395515E-01,
            5.31419325021531E-01,
            1.03129501638899E-01,
            1.22785291159756E-01,
            2.35702260395515E-01,
            2.37081084268193E-01,
            6.86947072008531E-01,
            6.86947072008531E-01,
            5.77350269189625E-01
        ]
        self.directionWeightsOct = [
            5.26559082615718E-02,
            9.95720041972116E-02,
            9.95720041972139E-02,
            8.80369927980968E-02,
            9.95720041972128E-02,
            5.26559082615709E-02,
            9.95720041972133E-02,
            8.80369927981024E-02,
            9.95720041972128E-02,
            9.95720041972108E-02,
            5.26559082615713E-02,
            8.80369927981055E-02,
            1.32024927825519E-01,
            1.32024927825525E-01,
            1.32024927825526E-01,
            1.55210814955923E-01
        ]


    def set_default_S4(self):
        """
        Quadrature set based on Chebichev Legendre S4
        """
        self.discreteDirectionsOct = 3
        self.directionUxOct = [
            3.5989E-01,
            8.6885E-01,
            3.5947E-01
        ]
        self.directionUyOct = [
            8.6885E-01,
            3.5989E-01,
            3.5947E-01
        ]
        self.directionUzOct = [
            3.3998E-01,
            3.3998E-01,
            8.6114E-01
        ]
        self.directionWeightsOct = [
            5.1219E-01,
            5.1219E-01,
            5.4641E-01
        ]


    def set_default_S8(self):
        """
        Quadrature set based on Chebichev Legendre S8
        """
        self.discreteDirectionsOct = 10
        self.directionUxOct = [
            1.9178E-01,
            5.4614E-01,
            8.1736E-01,
            9.6414E-01,
            2.2020E-01,
            6.0159E-01,
            8.2178E-01,
            2.3130E-01,
            5.5841E-01,
            1.9729E-01
        ]
        self.directionUyOct = [
            9.6414E-01,
            8.1736E-01,
            5.4614E-01,
            1.9178E-01,
            8.2178E-01,
            6.0159E-01,
            2.2020E-01,
            5.5841E-01,
            2.3130E-01,
            1.9729E-01
        ]
        self.directionUzOct = [
            1.8343E-01,
            1.8343E-01,
            1.8343E-01,
            1.8343E-01,
            5.2553E-01,
            5.2553E-01,
            5.2553E-01,
            7.9667E-01,
            7.9667E-01,
            9.6029E-01
        ]
        self.directionWeightsOct = [
            1.4243E-01,
            1.4243E-01,
            1.4243E-01,
            1.4243E-01,
            1.6426E-01,
            1.6426E-01,
            1.6426E-01,
            1.7466E-01,
            1.7466E-01,
            1.5901E-01
        ]


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter('discreteDirectionsOct', self.discreteDirectionsOct, isAddExtraLine=True)
        text += addParameter('directionUxOct', self.directionUxOct, isAddExtraLine=True)
        text += addParameter('directionUyOct', self.directionUyOct, isAddExtraLine=True)
        text += addParameter('directionUzOct', self.directionUzOct, isAddExtraLine=True)
        text += addParameter('directionWeightsOct', self.directionWeightsOct, isAddExtraLine=True)

        return(text)


    def plot_directions(self, point_size: float=20) -> None:
        """
        Plot the quadrature set on a sphere octant.

        Parameters
        ----------
        point_size : float
            Point size used for the plotting (default `20`).
        """
        plotter = pv.Plotter(off_screen=True)

        sphere = pv.Sphere(radius=1, end_theta=90, end_phi=90)
        plotter.add_mesh(sphere, opacity=0.3, use_transparency=True)

        points = np.transpose([
            self.directionUxOct,
            self.directionUyOct,
            self.directionUzOct
        ])

        # Normalize
        for i, point in enumerate(points):
            norm = np.sqrt(sum([e**2 for e in point]))
            points[i] = [p/norm for p in point]

        point_cloud = pv.PolyData(points)
        point_cloud["weight"] = self.directionWeightsOct
        plotter.add_mesh(
            point_cloud,
            point_size=point_size,
            render_points_as_spheres=True
        )

        # plotter.enable_parallel_projection()
        plotter.screenshot(f"fig_quadratureSet_{self.region}.png")
