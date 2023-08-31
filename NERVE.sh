#! /bin/bash -e
# script to run NERVE standalone trough docker

# create network
[ ! "$(docker network ls | grep nerve-network)" ] && docker network create nerve-network --attachable

# build psortb image
#[ ! "$(docker images -q psortb_http_api:v0.0.1)" ] && docker build -t psortb_http_api:v0.0.1 -f $(pwd)/docker/psortb/Dockerfile .

# run psortb container if not alredy running
[ ! "$(docker ps | grep francecosta/psortb_http_api:v0.0.1)" ] &&  docker run --rm -p 8080:8080 --network nerve-network --name psortb -d francecosta/psortb_http_api:v0.0.1

# build nerve if not already built
[ ! "$(docker images -q nerve:v0.0.3)" ] && docker build -t nerve:v0.0.3 .

# run nerve container
docker run --network nerve-network -p 8880:8880 -it -v $(pwd)/:/workdir nerve:v0.0.3 "$@"
