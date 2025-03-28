FROM francecosta/nerve:v0.0.7 AS final

COPY ./code /usr/nerve_python/NERVE/code
WORKDIR /workdir
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]