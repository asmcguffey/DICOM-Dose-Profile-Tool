import os

import pydicom

import numpy as np

from scipy.interpolate import RegularGridInterpolator


class RTDoseReader:

    def __init__(self, fp: str = '', ds: pydicom.Dataset = None):
        self._fp: str = fp
        self._ds: object = ds
        self._ipp: tuple = ()
        self._iop: tuple = ()
        self._dose_grid_scaling: float = None
        self._slice_thickness: float = None
        self._pixel_spacing: tuple = ()
        self._pixel_array: np.ndarray = None
        self._pixel_to_coords_matrix: np.ndarray = None
        self._interpolation_method: str = 'linear'
        self._interpolation_method_options = ('nearest', 'linear')
        self._interpolator = None
        self._dose_units: str = None

        # This code will prioritize the dataset `self._ds` over the file `self._fp` if the dataset is provided.
        if os.path.exists(self._fp) and not isinstance(self._ds, pydicom.Dataset):
            self._read_dataset_from_file()

        elif isinstance(self._ds, pydicom.Dataset):
            self._configure_dose_array()
    

    def _configure_dose_array(self):
        self._ipp = self._ds.ImagePositionPatient
        self._iop = self._ds.ImageOrientationPatient
        self._pixel_spacing = self._ds.PixelSpacing
        self._slice_thickness = self._ds.SliceThickness
        self._dose_grid_scaling = self._ds.DoseGridScaling
        self._dose_units = self._ds.DoseUnits
        self._pixel_array = self._ds.pixel_array
        self._pixel_array = self._pixel_array.transpose((2, 1, 0))
        self._pixel_to_coords_matrix = self._dicom_pixel_to_patient_slice(
            self._ipp,
            self._iop,
            self._pixel_spacing[0],
            self._pixel_spacing[1],
            self._slice_thickness
        )

    def _read_dataset_from_file(self):
        self._ds = pydicom.dcmread(self._fp)
        self._configure_dose_array()

    @staticmethod
    def _dicom_pixel_to_patient_slice(
            image_position_patient,
            image_orientation_patient,
            row_pixel_spacing,
            column_pixel_spacing,
            slice_thickness
    ):
        ipp_x, ipp_y, ipp_z = image_position_patient
        iop_x = image_orientation_patient[:3]
        iop_y = image_orientation_patient[3:]

        z_orientation = np.cross(iop_x, iop_y)

        # DICOM Standard Equation C.7.6.2.1-1
        # Maps pixel location to RCS
        transformation = np.array([
            [iop_x[0] * column_pixel_spacing, iop_y[0] * row_pixel_spacing, z_orientation[0] * slice_thickness, ipp_x],
            [iop_x[1] * column_pixel_spacing, iop_y[1] * row_pixel_spacing, z_orientation[1] * slice_thickness, ipp_y],
            [iop_x[2] * column_pixel_spacing, iop_y[2] * row_pixel_spacing, z_orientation[2] * slice_thickness, ipp_z],
            [0, 0, 0, 1]
        ])

        return transformation

    def _configure_interpolator(self):
        dose = self.dose_array
        dim = dose.shape
        self._interpolator = RegularGridInterpolator(
            [range(dim[0]), range(dim[1]), range(dim[2])],
            dose,
            method=self._interpolation_method,
            bounds_error=False,
            fill_value=0
        )

    def set_interpolation_method(self, method):
        if method not in self._interpolation_method_options:
            raise ValueError(
                'Interpolation method {method} not supported. Supported are "nearest" and "linear"')
        else:
            self._interpolation_method = method
            self._configure_interpolator()

    @property
    def dose_units(self):
        return self._dose_units

    @property
    def dose_array(self):
        return self._pixel_array * self._dose_grid_scaling
    
    @property
    def ijk_to_xyz(self):
        return self._pixel_to_coords_matrix
    
    @property
    def xyz_to_ijk(self):
        return np.linalg.inv(self._pixel_to_coords_matrix)

    def interpolate_points(self, dicom_pts):
        if self._interpolator is None:
            self._configure_interpolator()
        dicom_pts = np.hstack([dicom_pts, np.ones_like(dicom_pts[:, 0])[:, np.newaxis]])
        pts_ijk = (np.linalg.inv(self._pixel_to_coords_matrix) @ dicom_pts.T)
        return self._interpolator(pts_ijk[:3, :].T)

    def get_profile(self, p0=None, p1=None, spacing=1.0, return_coords=False):

        if p0 is not None and p1 is not None:
            if len(p0) != len(p1):
                raise AssertionError('p0 must be the same length as p1')

            else:
                p0 = np.asarray(p0)
                p1 = np.asarray(p1)

                x0, y0, z0 = p0
                x1, y1, z1 = p1

                dist = np.linalg.norm(p1 - p0)
                samples = int(dist // spacing) + 1

                xs = np.linspace(x0, x1, samples)
                ys = np.linspace(y0, y1, samples)
                zs = np.linspace(z0, z1, samples)

                line_xyz = np.vstack([xs, ys, zs, np.ones_like(xs)])
                line_ijk = (np.linalg.inv(self._pixel_to_coords_matrix) @ line_xyz)[:3, :]

                if self._interpolator is None:
                    self._configure_interpolator()

                if not isinstance(self._interpolator, RegularGridInterpolator):
                    raise ValueError('Interpolator must be instance of RegularGridInterpolator.')

                else:
                    dose_ijk = self._interpolator(line_ijk.T)
                    if not return_coords:
                        return dose_ijk
                    else:
                        return np.vstack([xs, ys, zs, dose_ijk]).T
