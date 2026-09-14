FROM francecosta/nerve:v0.0.7 AS final

ARG NERVE_SOURCE_COMMIT=unknown
ARG NERVE_SOURCE_DIRTY=unknown
ARG NERVE_SOURCE_TREE=unknown
LABEL org.opencontainers.image.source="https://github.com/nerve-bio/NERVE" \
      org.opencontainers.image.revision="${NERVE_SOURCE_COMMIT}" \
      org.opencontainers.image.source-tree="${NERVE_SOURCE_TREE}" \
      org.opencontainers.image.version="v0.0.8" \
      org.opencontainers.image.vendor="NERVE"
ENV NERVE_SOURCE_COMMIT="${NERVE_SOURCE_COMMIT}" \
    NERVE_SOURCE_DIRTY="${NERVE_SOURCE_DIRTY}" \
    NERVE_SOURCE_TREE="${NERVE_SOURCE_TREE}" \
    PYTHONPATH="/usr/nerve_python/NERVE:/usr/nerve_python"

COPY ./code /usr/nerve_python/NERVE/code
WORKDIR /workdir
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]
