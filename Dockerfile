FROM hroest/openms-lib-2.2 
MAINTAINER Bernie Pope
ENV BASE_DIR /src
RUN apt-get -y update && rm -rf /var/lib/apt/lists/*
ADD ./score $BASE_DIR/score
WORKDIR $BASE_DIR/score
RUN cmake -D OpenMS_DIR=/openms_build/ -D CMAKE_PREFIX_PATH="/contrib-build/;/usr/;/usr/local" .
RUN make
WORKDIR $BASE_DIR
COPY ./hitime /
ENTRYPOINT ["/hitime"]
