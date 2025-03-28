#! /bin/bash -e
# script to run NERVE standalone trough docker

PSORT_VERSION="v0.0.1"
NERVE_VERSION="v0.0.8"

# create network
[ ! "$(docker network ls | grep nerve-network)" ] && docker network create nerve-network --attachable

# build psortb image
#[ ! "$(docker images -q psortb_http_api:v0.0.1)" ] && docker build -t psortb_http_api:v0.0.1 -f $(pwd)/docker/psortb/Dockerfile .

# run psortb container if not alredy running
[ ! "$(docker ps | grep francecosta/psortb_http_api:$PSORT_VERSION)" ] &&  docker run --rm -p 8080:8080 --network nerve-network --name psortb -d francecosta/psortb_http_api:$PSORT_VERSION

# build nerve if not already built
[ ! "$(docker images -q nerve:$NERVE_VERSION)" ] && docker build -t nerve:$NERVE_VERSION .

# run nerve container
docker run --network nerve-network -p 8880:8880 -it -v $(pwd):/workdir nerve:$NERVE_VERSION "$@"
