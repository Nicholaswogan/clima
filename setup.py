from skbuild import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open("CMakeLists.txt",'r') as fil:
    lines = fil.readlines()
    for line in lines:
        if line.startswith("project(Clima"):
            version = line.split('"')[1]
            break
            
setup(
    name="clima",
    packages=['clima'],
    python_requires='>=3.6',
    version=version,
    license="GNU General Public License v3.0",
    install_requires=['numpy','numba'], 
    author='Nicholas Wogan',
    author_email = 'nicholaswogan@gmail.com',
    description = "Climate model of planet's atmospheres.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = "https://github.com/Nicholaswogan/clima",
    cmake_args=['-DCMAKE_BUILD_TYPE=Release',\
                '-DBUILD_PYTHON_CLIMA=ON',\
                '-DBUILD_EXECUTABLES=OFF',\
                '-DBUILD_WITH_OPENMP=ON']
)


