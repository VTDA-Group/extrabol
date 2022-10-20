import setuptools
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'docs/README.md'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

setuptools.setup(
    name="extrabol", # Replace with your own username
    version="0.1.14",
    py_modules=['extrabol'],
    author="Ian M. Thornton",
    author_email="iot5037@psu.edu",
    description="Estimate SN bolometric light curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/villrv/extrabol",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'extrabol': ['template_bank/*.npz', 'example/PSc000174.dat']},
    install_requires=[
        "numpy",
        "astropy",
        "astroquery",
        "matplotlib",
        "george",
        "extinction",
        "emcee",
        "mkdocs >= 1.2.2",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
            'console_scripts': [
                    'extrabol = extrabol.extrabol:main'
                ]},
    python_requires='>=3.6'
)