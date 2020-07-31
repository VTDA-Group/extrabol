import setuptools
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()


setuptools.setup(
    name="extrabol", # Replace with your own username
    version="0.0.10",
    author="V. Ashley Villar",
    author_email="vav2110@columbia.edu",
    description="Estimate SN bolometric light curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/villrv/extrabol",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
            'console_scripts': [
                    'extrabol = extrabol.extrabol:main',
                ]},
    python_requires='>=3.6'
)