FROM francecosta/nerve:v0.0.5 AS final

COPY ./ /usr/nerve_python/NERVE
WORKDIR /workdir
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]
