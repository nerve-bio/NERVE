###
### NERVE dependencies installation
###
### base image found here https://hub.docker.com/_/python
### first stage (base, with Python 3.8)
FROM python:3.8 AS base
 # set where NERVE is contained
ENV PYTHONPATH "${PYTHONPATH}:/NERVE"
# set where spaan is contained
ENV PYTHONPATH "${PYTHONPATH}:/"
ENV LANG C

FROM base AS intermediate
# install python dependencies
COPY requirements.txt .
RUN pip install --upgrade setuptools
RUN pip install -r ./requirements.txt && \
    python -m pip install git+https://github.com/nicolagulmini/tmhmm.py

FROM intermediate AS dependencies
# clone repositories
RUN git clone https://github.com/flatironinstitute/DeepFRI.git
# download model weights for DeepFri
RUN wget https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz
RUN tar xvzf trained_models.tar.gz -C ./DeepFRI
RUN mv trained_models ./DeepFri/trained_models
RUN git clone https://github.com/Superzchen/iFeature.git
#COPY trained_models ./trained_models
#COPY DeepFRI ./DeepFri
#COPY iFeature ./iFeature

RUN git clone https://github.com/nicolagulmini/spaan && \
    git clone https://github.com/FranceCosta/NERVE && \
    python3 ./DeepFri/setup.py install

RUN apt-get update && \
    apt-get install -y apt-utils ncbi-blast+ &&\
    apt-get install nano

FROM dependencies AS final
USER user
RUN mkdir /workdir
RUN chown user:user /workdir

ENTRYPOINT ["/NERVE/NERVE/NERVE_main.py"]

#ARG USER_ID
#ARG GROUP_ID

#EXPOSE 8888

# user and group id are defined in build.sh command

#RUN addgroup --gid $GROUP_ID user # name user as user
#RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user # name group user as user

#RUN chown user:user /NERVE



#CMD jupyter notebook --notebook-dir=/my_data --ip='0.0.0.0' --port='8888' --allow-root --NotebookApp.token='' --NotebookApp.password=''

