# Docker environment for ubuntu, python3.6
#
# Usage:
#  * build the image:
#    $ docker build -t pynta .
#  * start the image:
#    docker run -it pynta

# Latest version of ubuntu
FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

# Install system packages
RUN apt-get update && \
    apt-get install -y git openmpi-bin libopenmpi-dev && \
    apt-get install -y wget cmake meson libopenblas-dev

# Install python

RUN apt-get install -y python3.6 python3-pip

# Install Python packages

RUN pip3 install --upgrade pip && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install flake8 pytest && \
    python3 -m pip install numpy scipy cffi

WORKDIR /build

# Install xtb-python

RUN git clone https://github.com/grimme-lab/xtb-python.git && \
    cd xtb-python %% \
    git submodule update --init %% \
    meson setup build --prefix=$PWD --libdir=xtb/xtb --buildtype release --optimization 2 -Dla_backend=openblas %% \
    ninja -C build install && \
    pip install --user -e .

# Install PostgreSQL
RUN wget https://sbp.enterprisedb.com/getfile.jsp?fileid=1257417&_ga=2.119848136.1712599154.1611463162-349055021.1611463162 && \
    wait && \
    mv 'getfile.jsp?fileid=1257417' postgresql.tar.gz && \
    tar -xvzf postgresql.tar.gz && \
    rm postgresql.tar.gz && \
    echo 'export PATH=./pgsql/bin:$PATH' >> ~/.bashrc

# Install Balsam

RUN git clone https://github.com/balsam-alcf/balsam.git && \
    cd balsam && \
    python3 setup.py install

# RUN git clone https://github.com/zadorlab/sella.git &&\
#     cd /build/sella && \
#     git reset --hard 463e0556089f6af3385bed4aa7b49ee040cd19d1 && \
#     python3 setup.py install

# Set pynta code path
ENV CODE_DIR /build/pynta

# Copy pynta source files

COPY . $CODE_DIR

# Install pynta

RUN cd $CODE_DIR && \
    python3 setup.py install

# RUN mkdir -p $CODE_DIR
