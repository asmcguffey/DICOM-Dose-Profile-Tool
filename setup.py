from setuptools import setup

setup(
    name='dicom-dose-profile-tool',
    version='0.1.0',
    description='Extracts dose profiles from a DICOM RT dataset.',
    author='Andrew McGuffey',
    author_email='amcguf3@lsu.edu',
    url='https://github.com/asmcguffey/DICOM-Dose-Profile-Tool'
    packages=['dicom-dose-profile-tool'],
    python_requires=">=3.6"
)