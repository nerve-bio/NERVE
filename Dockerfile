# Pinned by digest: the v0.0.7 tag is mutable. Re-pin with
# `docker image inspect francecosta/nerve:v0.0.7 --format '{{.RepoDigests}}'`.
FROM francecosta/nerve:v0.0.7@sha256:f49feae627370f12bde8a4ec12a117b53e5840caf10a21531e7fc1e9c46fd9f9 AS final

ARG NERVE_SOURCE_COMMIT=unknown
ARG NERVE_SOURCE_DIRTY=unknown
ARG NERVE_SOURCE_TREE=unknown
LABEL org.opencontainers.image.source="https://github.com/nerve-bio/NERVE" \
      org.opencontainers.image.revision="${NERVE_SOURCE_COMMIT}" \
      org.opencontainers.image.source-tree="${NERVE_SOURCE_TREE}" \
      org.opencontainers.image.version="v0.1.0" \
      org.opencontainers.image.vendor="NERVE"
ENV NERVE_SOURCE_COMMIT="${NERVE_SOURCE_COMMIT}" \
    NERVE_SOURCE_DIRTY="${NERVE_SOURCE_DIRTY}" \
    NERVE_SOURCE_TREE="${NERVE_SOURCE_TREE}" \
    PYTHONPATH="/usr/nerve_python/NERVE:/usr/nerve_python"

COPY ./code /usr/nerve_python/NERVE/code
RUN find /usr/nerve_python/NERVE/code -type d -name __pycache__ -prune -exec rm -rf {} +
WORKDIR /workdir
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]
