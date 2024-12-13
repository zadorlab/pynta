from setuptools import setup, find_packages

VERSION = '0.0.0'
DESCRIPTION = 'Connection between ASE + PWDFT'
LONG_DESCRIPTION = 'Module that runs PWDFT using ASE (calculator), made in ANL'

# Setting up
setup(
        name="ase_pwdft",
        version=VERSION,
        author="RayHEsparza",
        author_email="rayhe88@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[],
        keywords=['ASE','PWDFT','Electronic Structure'],
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: Linux"
            ] 
    )
