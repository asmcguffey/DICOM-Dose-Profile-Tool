import os

import numpy as np

from parse_dicom import parse_folder
from read_rtdose import RTDoseReader
from read_rtplan import RTPlanReader


class ProfileTool:

    def __init__(self, folder=None, interpolation_method='linear'):

        self._rtdose = None
        self._rtplan = None
        self._interpolation_method = interpolation_method

        if isinstance(folder, str):
            if os.path.isdir(folder):
                self._folder = folder
            else:
                raise FileNotFoundError('Provided folder does not exist.')

            self._get_datasets_from_folder()
    
    def load_folder(self, folder):

        self._folder = folder
        self._get_datasets_from_folder()
    
    
    def _get_datasets_from_folder(self):

        required_modalities = ('RTDOSE', 'RTPLAN')
        dict_of_datasets = parse_folder(self._folder)

        if dict_of_datasets:

            if all([(modality in dict_of_datasets) for modality in required_modalities]):
                
                self._rtdose = dict_of_datasets['RTDOSE']
                self._rtplan = dict_of_datasets['RTPLAN']
                self._configure_readers()

    def _configure_readers(self):
        
        self._plan_reader = RTPlanReader(ds=self._rtplan)
        self._dose_reader = RTDoseReader(ds=self._rtdose)
        self._dose_reader.set_interpolation_method(self._interpolation_method)

    def get_profile(
        self, p0, p1, spacing=0.1, return_coords=True,
        profile_dose_units='Gy', normalize_by_monitor_units = False,
        zero_point='default'):

        """
        p0: array-like
            DICOM x, y, z coordinates of the starting point of the profile
        p1: array-like
            DICOM x, y, z coordinates of the ending point of the profile
        return_coords: bool default True
            Whether to return coordinates in the output array of the function.
            If set to False, the coordinates can be made externally using np.linspace() 
            over each of the dimensions.
        profile_dose_units: str default 'Gy'
            The units to use for the dose profile. Options are 'Gy' and
        zero_point: array-like or str or None
            DICOM x, y, z coordinates of the zero point of the profile. 
        """

        spacing = float(spacing)
        spacing *= 10

        if zero_point is None:
            zero_point = np.asarray([0, 0, 0])

        elif zero_point == 'default':
            zero_point = np.asarray(self.get_zero_point(mode='INTERCEPT_CAX_SURFACE'))

        else:
            zero_point = np.asarray(zero_point)
        
        if isinstance(self._dose_reader, RTDoseReader):

            profile = self._dose_reader.get_profile(
                p0=p0 + zero_point,
                p1=p1 + zero_point,
                spacing=spacing,
                return_coords=return_coords
            )
            profile[:, :3] = profile[:, :3] - zero_point[None, :]

            if self._dose_reader.dose_units == 'GY':
                if profile_dose_units.upper() == 'CGY':
                    profile[:, 3] = profile[:, 3] * 100
            
            elif self._dose_reader.dose_units == 'RELATIVE':
                raise ValueError('DICOM dose units == "RELATIVE". Relative dose not supported.')

            else:
                raise ValueError(f'Invalid dose units {self._dose_reader.dose_units}')
                
            if normalize_by_monitor_units:
                if isinstance(self._plan_reader, RTPlanReader):
                    if 'Monitor Units' in self._plan_reader.beam_parameters:
                        monitor_units = float(self._plan_reader.beam_parameters['Monitor Units'])
                        profile[:, 3] = profile[:, 3] / monitor_units
            
            return profile

        else:
            return None

    def get_zero_point(self, mode='INTERCEPT_CAX_SURFACE') -> list:
        """
        Parameters
        ----------
        mode: str default 'INTERCEPT_CAX_SURFACE'
            Mode by which to get zero point.
                'INTERCEPT_CAX_SURFACE'
                Calculates the coordinates of the point at which the beam's central axis crosses the surface.

        """
        if mode == 'INTERCEPT_CAX_SURFACE':
            # The gantry rotates in the XY Plane. 
            beam_params_dict = self._plan_reader.beam_parameters
            if not all([key in beam_params_dict for key in ['Isocenter', 'Gantry Angle', 'SSD']]):
                return None
            else:
                iso = [float(v) for v in beam_params_dict['Isocenter']]
                gantry_angle = float(beam_params_dict['Gantry Angle'])
                ssd = float(beam_params_dict['SSD'])
                sad = float(beam_params_dict['SAD'])

                x_zero_pt = iso[0] - (ssd - sad) * np.sin(np.radians(gantry_angle))
                y_zero_pt = iso[1] + (ssd - sad) * np.cos(np.radians(gantry_angle))
                z_zero_pt = iso[2]

                return [x_zero_pt, y_zero_pt, z_zero_pt]



    @property
    def dicom_datasets(self) -> dict:
        return {'RTDOSE': self._rtdose, 'RTPLAN': self._rtplan}

    def profiles_from_file(
        self, fp, save_as, delimiter=',', skip_rows=2,
        spacing=0.1, profile_dose_units='Gy', 
        normalize_by_monitor_units = False, zero_point='default'):
        """
        Parameters
        ----------
        fp: str
            The file path (must be csv).
        save_as: str
            The path to save the data to. File d
        delimiter: str default ','
            The delimiter of the data file.
        spacing: float default 0.1
            The point spacing in cm of the interpolation.
        Returns
        -------
        None
            Writes csv file with profile data.
        """

        profiles = []


        if isinstance(self._dose_reader, RTDoseReader):

            with open(fp, 'r') as csv_file:
                for i, line in enumerate(csv_file.readlines()):
                    if i > skip_rows - 1:
                        points = [float(v) for v in line.strip().split(delimiter) if v != '']
                        z0, z1, x0, x1, y0, y1 = points
                        profile = self.get_profile(
                            p0=np.asarray([x0, y0, z0]) * 10,
                            p1=np.asarray([x1, y1, z1]) * 10,
                            spacing=spacing,
                            profile_dose_units=profile_dose_units,
                            normalize_by_monitor_units=normalize_by_monitor_units,
                            zero_point=zero_point)

                        profile[:, :3] = profile[:, :3] / 10

                        profiles.append(profile)
            
            if os.path.exists(save_as):
                os.remove(save_as)

            if normalize_by_monitor_units:
                profile_dose_units = profile_dose_units + '/MU'

            write_header = ['Crossline (X)', 'Inline (Z)', 'Depth (Y)', f'Dose ({profile_dose_units})']
            with open(save_as, 'w') as write_file:
                ext = os.path.splitext(save_as)[-1].lower()
                if ext == '.csv':
                    write_dlm = ','
                else:
                    raise ValueError('Only writing to CSV files is supported.')

                for i, profile in enumerate(profiles):
                    write_file.write(str(i + 1) + '\n')
                    write_file.write(write_dlm.join(write_header) + '\n')
                    for i in range(profile.shape[0]):
                        write_file.write(
                            write_dlm.join([f'{v:.3f}' for v in profile[i, :]]) + '\n')
                    write_file.write('\n\n')
        

                    


            

            

        


                    
        






