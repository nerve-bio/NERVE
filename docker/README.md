NERVE developer version
This setting is based on docker-compose. It permits to run multiple docker containers parallely.
- ./web: contains NERVE and its dependencies;
- ./psortb: contains psortb, its dependencies and an API to allow calling the program from within the web container.
- ./workdir is a shared "volume" between the computer and the web container that conists on a folder seen both by the container itself and the machine in 
which it is run.
./web and ./psortb contains folder contains a "Dockerfile" that define the respective environments.

Usage:
- download and install docker-compose;
- $tar –xvzf NERVE_compose.tar.gz
- $cd NERVE_compose
- $./build_run.sh
- navigate to http://localhost:8880

Notebooks:
- executable_tests.ipynb: contains some example code to run NERVE
