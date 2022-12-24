###
### NERVE dependencies installation
###
### base image found here https://hub.docker.com/_/python
### first stage (base, with Python 3.8)
ARG VERSION
FROM python:3.8 AS base
 # set where NERVE is contained
ENV PYTHONPATH "${PYTHONPATH}:/NERVE"
# set where spaan is contained
ENV PYTHONPATH "${PYTHONPATH}:/"
ENV LANG C

FROM base AS intermediate

RUN pip install setuptools==65.4.1
RUN git clone https://github.com/flatironinstitute/DeepFRI.git
RUN mv ./DeepFRI ./DeepFri 
# install DeepFri dependencies
WORKDIR ./DeepFri
RUN pip install .
WORKDIR /
# download model weights for DeepFri
RUN wget https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz
RUN mv newest_trained_models.tar.gz trained_models.tar.gz
RUN tar xvzf trained_models.tar.gz -C ./DeepFri

# install python dependencies
COPY requirements.txt .
RUN pip install -r ./requirements.txt && \
    python -m pip install git+https://github.com/nicolagulmini/tmhmm.py

FROM intermediate AS dependencies
# clone repositories
RUN git clone https://github.com/Superzchen/iFeature.git && \
    git clone https://github.com/FranceCosta/NERVE && \
    git clone https://github.com/nicolagulmini/spaan

RUN apt-get update && \
    apt-get install -y apt-utils ncbi-blast+ &&\
    apt-get install nano

FROM dependencies AS setting
ARG USER_ID
ARG GROUP_ID
EXPOSE 8880
RUN addgroup --gid $GROUP_ID user # name user as user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user # name group user as user
RUN mkdir /workdir
RUN chown -R user:user /workdir
# give permissions to allow testing
RUN chown -R user:user /NERVE
RUN chmod +x /NERVE/code/NERVE.py
USER user

FROM setting as version-interactive
ENV cmd="jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''"

FROM setting as version-commandline
ENV cmd="['/NERVE/NERVE/NERVE_main.py']"

FROM version-${VERSION} AS final
#RUN echo "command is ${cmd}"
#CMD cmd
CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''

#CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''
#ENTRYPOINT ["/NERVE/NERVE/NERVE_main.py"]

