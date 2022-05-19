from parse_dicom import parse_folder
from profiles_from_dicom import ProfileTool
import PySimpleGUI as sg


def popup_window(layout, title, window_kw):
    if "icon" not in window_kw:
        window_kw["icon"] = "icon.ico"
    window = sg.Window(title, layout, **window_kw)
    while True:
        event, values = window.read()
        if event in ["Exit", sg.WIN_CLOSED]:
            break


def popup_view_dicom_info(dict_of_datasets):
    multiline_kw = {
        "size": (100, 20),
        "expand_x": True,
        "expand_y": True,
        "auto_refresh": True,
    }
    tabs = sg.TabGroup(
        [
            [
                sg.Tab(key, [[sg.Multiline(str(dataset), **multiline_kw)]])
                for key, dataset in dict_of_datasets.items()
            ]
        ]
    )
    layout = [[tabs]]
    popup_window(layout, title="View DICOM Info", window_kw={"resizable": True})


def display():

    sg.theme("DarkGrey13")

    config_frame_layout = [
        [
            sg.Text("DICOM Folder"),
            sg.Push(),
            sg.Input(
                key="-DICOM FOLDER PATH-",
                tooltip="Path to DICOM folder. Must contain \nexactly 1 RTDOSE and 1 RTPLAN file.",
            ),
            sg.FolderBrowse("Browse..."),
        ],
        [
            sg.Text("Points File (.csv)"),
            sg.Push(),
            sg.Input(
                key="-POINTS FILE PATH-", tooltip="Path to csv file defining profiles."
            ),
            sg.FileBrowse("Browse...", file_types=(("CSV Files", "*.csv"))),
        ],
    ]

    interpolation_options = {"Nearest Neighbor": "nearest", "Linear": "linear"}

    interpolation_frame_layout = [
        [
            sg.Text("Spacing (cm)"),
            sg.Push(),
            sg.Input(default_text=0.1, size=(5, 1), key="-SPACING-"),
            sg.OptionMenu(
                list(interpolation_options.keys()),
                default_value="Linear",
                key="-INTERPOLATION METHOD-",
                tooltip="The method of interpolation to use on the DICOM RTDOSE 3D grid.",
            ),
        ]
    ]

    coordinates_frame_layout = [
        [
            sg.Checkbox(
                "Set Zero Point at Surface / CAX Intersection",
                default=True,
                key="-CENTER AT SURFACE CAX-",
                tooltip="Whether to set the zero point (0, 0, 0) of the coordinates \n at the intersection of beam central axis and water surface.",
            )
        ]
    ]

    dose_unit_options = ("cGy", "Gy")

    dose_output_settings_layout = [
        [
            sg.Text("Dose units"),
            sg.OptionMenu(dose_unit_options, default_value="Gy", key="-DOSE UNIT-"),
            sg.Checkbox(
                "Normalize by MU",
                default=False,
                key="-NORMALIZE BY MU-",
                tooltip="Whether to divide the dose values by the MU setting of the beam.",
            ),
        ]
    ]

    layout = [
        [sg.Frame("Configuration", layout=config_frame_layout)],
        [sg.Frame("Interpolation", layout=interpolation_frame_layout)],
        [sg.Frame("Coordinates", layout=coordinates_frame_layout)],
        [sg.Frame("Output", layout=dose_output_settings_layout)],
        [
            sg.In(visible=False, key="-SAVE AS-", enable_events=True),
            sg.FileSaveAs(
                "Write Profiles",
                target="-SAVE AS-",
                file_types=(("CSV Files", ".csv"),),
            ),
            sg.Cancel(),
            sg.Push(),
            sg.B("View DICOM Info"),
        ],
    ]

    window = sg.Window("DICOM Dose Profile Tool", layout, icon="icon.ico")

    pt = None

    while True:
        event, values = window.read()
        if event == "-SAVE AS-":
            if all(
                [
                    values[key]
                    for key in [
                        "-DICOM FOLDER PATH-",
                        "-POINTS FILE PATH-",
                        "-SAVE AS-",
                    ]
                ]
            ):
                pt = ProfileTool(
                    folder=values["-DICOM FOLDER PATH-"],
                    interpolation_method=interpolation_options[
                        values["-INTERPOLATION METHOD-"]
                    ],
                )

                pt.profiles_from_file(
                    fp=values["-POINTS FILE PATH-"],
                    save_as=values["-SAVE AS-"],
                    skip_rows=2,
                    spacing=values["-SPACING-"],
                    profile_dose_units=values["-DOSE UNIT-"],
                    normalize_by_monitor_units=values["-NORMALIZE BY MU-"],
                    zero_point="default" if values["-CENTER AT SURFACE CAX-"] else None,
                )

        if event in (sg.WIN_CLOSED, "Cancel"):
            break

        if event == "View DICOM Info":
            if values["-DICOM FOLDER PATH-"]:
                dicom_datasets = parse_folder(values["-DICOM FOLDER PATH-"])
                popup_view_dicom_info(dicom_datasets)

    window.close()
