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
#RUN pip install nbconvert==7.2.0 --force-reinstall

FROM intermediate AS dependencies
# clone repositories
RUN git clone https://github.com/flatironinstitute/DeepFRI.git
RUN mv ./DeepFRI ./DeepFri 
RUN python3 ./DeepFri/setup.py install
# download model weights for DeepFri
RUN wget https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz
RUN mv newest_trained_models.tar.gz trained_models.tar.gz
RUN tar xvzf trained_models.tar.gz -C ./DeepFri
RUN git clone https://github.com/Superzchen/iFeature.git && \
    git clone https://github.com/FranceCosta/NERVE && \
    git clone https://github.com/nicolagulmini/spaan

RUN apt-get update && \
    apt-get install -y apt-utils ncbi-blast+ &&\
    apt-get install nano

FROM dependencies AS final
ARG USER_ID
ARG GROUP_ID
EXPOSE 8880
RUN addgroup --gid $GROUP_ID user # name user as user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user # name group user as user
RUN mkdir /workdir
RUN chown -R user:user /workdir
#VOLUME ["/workdir"]
RUN chmod +x /NERVE/NERVE/NERVE_main.py
USER user

CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''
#ENTRYPOINT ["/NERVE/NERVE/NERVE_main.py"]

