#docker-compose up -d --build --no-cache
#docker-compose up --force-recreate
#docker-compose down && docker-compose build --no-cache && docker-compose up
docker-compose down && docker-compose build --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) && docker-compose up -d

