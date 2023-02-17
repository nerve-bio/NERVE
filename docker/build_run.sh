#docker-compose up -d --build --no-cache
#docker-compose up --force-recreate
#docker-compose down && docker-compose build --no-cache && docker-compose up -d
docker-compose down && docker-compose build --no-cache --build-arg USER_ID="$(id -u)" --build-arg GROUP_ID="$(id -g)" --build-arg VERSION="interactive" && docker-compose up -d

