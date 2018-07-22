FROM hroest/openms-lib-2.2 
MAINTAINER Bernie Pope
ENV BASE_DIR /src
RUN add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get -y update && \
    apt-get install -y gcc-4.9 g++-4.9 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 && \
    rm -rf /var/lib/apt/lists/*
ADD ./score $BASE_DIR/score
WORKDIR $BASE_DIR/score
RUN cmake -D OpenMS_DIR=/openms_build/ -D CMAKE_PREFIX_PATH="/contrib-build/;/usr/;/usr/local" .
RUN make
WORKDIR $BASE_DIR
COPY ./hitime /
ENTRYPOINT ["/hitime"]
