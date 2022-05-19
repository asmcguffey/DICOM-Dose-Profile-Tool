import os
import pydicom


class RTPlanReader:
    def __init__(self, fp=None, ds=None):

        self._fp = fp
        self._ds = ds

        self._beam_parameters: dict = {
            "Isocenter": {
                "tag": [(0x300A, 0x00B0), (0x300A, 0x0111), (0x300A, 0x012C)],
            },
            "Monitor Units": {"tag": [(0x300A, 0x00B0), (0x300A, 0x010E)]},
            "SSD": {"tag": [(0x300A, 0x00B0), (0x300A, 0x0111), (0x300A, 0x0130)]},
            "SAD": {"tag": [(0x300A, 0x00B0), (0x300A, 0x00B4)]},
            "Gantry Angle": {
                "tag": [(0x300A, 0x00B0), (0x300A, 0x0111), (0x300A, 0x011E)]
            },
        }

        self._iso = None
        self._monitor_units = None
        self._source_to_surface_distance = None

        self.load_dataset()

    def load_file(self, fp):
        self._fp = fp
        self.load_dataset()

    def load_dataset(self):
        if isinstance(self._ds, pydicom.Dataset):
            self.read_beam_parameters_from_dataset()
        else:
            if self._fp is not None:
                if isinstance(self._fp, str):
                    if os.path.exists(self._fp):
                        self._ds = pydicom.dcmread(self._fp)
                        self.read_beam_parameters_from_dataset()

    def read_beam_parameters_from_dataset(self):
        """
        Reads the beam parameters from an RT Plan file. Tag is provided as a list of tuples.
        The DICOM metadata uses a tree-like structure (e.g., A block is inside a beam, which is inside a plan.)
        The for-loop is used to access data nested in the tree using the list of tags.
        Tags are arranged in the list ordered from highest-level to lowest-level.
        """
        for key in self._beam_parameters:
            contents = self._ds
            tag_string = self._beam_parameters[key]["tag"]
            for i, tag in enumerate(tag_string):
                if i < len(tag_string) - 1:
                    contents = contents[tag][0]
                else:
                    contents = contents[tag]

            self._beam_parameters[key]["value"] = contents.value

    @property
    def iso(self):
        return self._iso

    @property
    def ds(self):
        return self._ds

    @property
    def monitor_units(self):
        return self._monitor_units

    @property
    def beam_parameters(self) -> dict:
        return {
            key: items["value"]
            for key, items in self._beam_parameters.items()
            if "value" in items
        }
