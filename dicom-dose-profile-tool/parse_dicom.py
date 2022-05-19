import pydicom
import os


def parse_folder(folder: str) -> dict:
    """
    Parameters
    ----------
    folder: str
        The folder which contains dicom files.

    Returns
    -------
    dict
        Dictionary of pydicom.Dataset objects.
    """
    out = dict()

    ignore_modalities = ["CT"]

    if os.path.isdir(folder):
        for root, _, files in os.walk(folder):
            for file in files:
                path = os.path.join(root, file)
                if file.endswith(".dcm"):
                    tmp_ds = pydicom.dcmread(path)
                    modality = tmp_ds.Modality
                    # Ignore CT images
                    if modality not in ignore_modalities:
                        out[modality] = tmp_ds

    return out
