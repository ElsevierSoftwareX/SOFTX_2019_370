FROM ubuntu:trusty 
MAINTAINER Bernie Pope
ENV BASE_DIR /src
ADD ./score $BASE_DIR/score
ADD ./maxima $BASE_DIR/maxima
RUN apt-get -y update && apt-get install -y \
   autoconf \
   automake \
   cmake \
   git \
   g++ \
   libtool \
   make \
   libboost-all-dev \
   libbz2-dev \
   libglpk-dev \
   libqt4-dev \
   libsvm-dev \
   libxerces-c-dev \
   libzip-dev \
   patch \
   zlib1g-dev \
   && rm -rf /var/lib/apt/lists/*
WORKDIR $BASE_DIR 
RUN mkdir openms
WORKDIR $BASE_DIR/openms
RUN git clone https://github.com/OpenMS/contrib.git
RUN mkdir contrib-build
WORKDIR $BASE_DIR/openms/contrib-build
RUN cmake -DBUILD_TYPE=SEQAN ../contrib
RUN cmake -DBUILD_TYPE=WILDMAGIC ../contrib
RUN cmake -DBUILD_TYPE=EIGEN ../contrib
WORKDIR $BASE_DIR/openms
RUN git clone https://github.com/OpenMS/OpenMS.git
RUN mkdir openms_build
WORKDIR $BASE_DIR/openms/openms_build
RUN cmake -DCMAKE_PREFIX_PATH="$BASE_DIR/openms/contrib-build;/usr;/usr/local" -DBOOST_USE_STATIC=OFF -DWITH_GUI=OFF -DHAS_XSERVER=OFF ../OpenMS
RUN make
WORKDIR $BASE_DIR/score
RUN cmake -D OpenMS_DIR=$BASE_DIR/openms/openms_build/ -D CMAKE_PREFIX_PATH="$BASE_DIR/openms/contrib-build;/usr;/usr/local" .
RUN make
WORKDIR $BASE_DIR/maxima
RUN cmake -D OpenMS_DIR=$BASEDIR/openms/openms_build/ -D CMAKE_PREFIX_PATH="$BASE_DIR/openms/contrib-build;/usr;/usr/local" .
RUN make
COPY ./hitime /
ENTRYPOINT ["/hitime"]
