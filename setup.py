import setuptools
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

with open(os.path.join(here, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]


setuptools.setup(
    name="extrabol", # Replace with your own username
    version="0.0.30",
    py_modules=['extrabol', 'snana2extrabol'],
    author="Ian M. Thornton",
    author_email="iot5037@psu.edu",
    description="Estimate SN bolometric light curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/villrv/extrabol",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
            'console_scripts': [
                    'extrabol = extrabol.extrabol:main',
                    'snana2extrabol = extrabol.snana2extrabol:convert'
                ]},
    python_requires='>=3.6'
)