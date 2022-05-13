import pydicom
import os


def parse_folder(folder: str) -> dict:

    out = dict()
    
    if os.path.isdir(folder):
        for root, _, files in os.walk(folder):
            for file in files:
                path = os.path.join(root, file)
                if file.endswith('.dcm'):
                    tmp_ds = pydicom.dcmread(path)
                    modality = tmp_ds.Modality
                    if modality != 'CT':
                        out[modality] = tmp_ds
    
    return out


