FROM hroest/openms-lib-2.2 
MAINTAINER Bernie Pope
ENV BASE_DIR /src
RUN apt-get -y update && apt-get install -y \
    libboost-program-options-dev \ 
    && rm -rf /var/lib/apt/lists/*
ADD ./score $BASE_DIR/score
ADD ./maxima $BASE_DIR/maxima
WORKDIR $BASE_DIR/score
RUN cmake -D OpenMS_DIR=/openms_build/ -D CMAKE_PREFIX_PATH="/contrib-build/;/usr/;/usr/local" .
RUN make
#WORKDIR $BASE_DIR/maxima
#RUN cmake -D OpenMS_DIR=/openms_build/ -D CMAKE_PREFIX_PATH="/contrib-build/;/usr/;/usr/local" .
#RUN make
WORKDIR $BASE_DIR
COPY ./hitime /
ENTRYPOINT ["/hitime"]
