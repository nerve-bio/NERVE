### Build the image, passing current User and Group (numeric)
### so that the built image will use those values for managing
### files in mounted folders
docker build --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) -t --no-cache nerve .
