FROM python:3.8 AS base

LABEL maintainer="Francesco Costa <fcosta@ebi.ac.uk>"
ENV PYTHONPATH "${PYTHONPATH}:/usr/nerve_python"
ENV PYTHONPATH "${PYTHONPATH}:/usr/nerve_python/NERVE"
ENV LANG C

FROM base AS dependencies
WORKDIR /usr/nerve_python
RUN pip install setuptools==65.4.1
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install git-lfs && \
    git lfs install && \
    git lfs clone https://github.com/FranceCosta/nerve_data.git && \
    git clone https://github.com/nicolagulmini/spaan && \
    git clone https://github.com/FranceCosta/NERVE

WORKDIR /usr/nerve_python/assets

RUN mv /nerve_data/DeepFri.tar.gz DeepFri.tar.gz && \
    tar xvzf ./DeepFri.tar.gz && \
    mv /nerve_data/newest_trained_models.tar.gz trained_models.tar.gz && \
    tar xvzf trained_models.tar.gz -C ./DeepFri && \
    mv /nerve_data/iFeature.tar.gz iFeature.tar.gz && \
    tar xvzf iFeature.tar.gz

WORKDIR /usr/nerve_python/assets/DeepFri
RUN pip install .

WORKDIR /usr/nerve_python

RUN pip install -r ./requirements.txt && \
    python -m pip install git+https://github.com/nicolagulmini/tmhmm.py && \
    apt-get update && \
    apt-get install -y apt-utils ncbi-blast+ nano

FROM dependencies AS final
