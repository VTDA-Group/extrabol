import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="extrabol", # Replace with your own username
    version="0.0.3",
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
    python_requires='>=3.6',
)