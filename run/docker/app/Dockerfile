FROM python:3.6-stretch
ENV PYTHONUNBUFFERED 1

RUN apt-get update
RUN apt-get install -y curl
RUN apt-get install -y shellcheck
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ stretch-pgdg main" >> /etc/apt/sources.list
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ xenial-pgdg main 12" | tee /etc/apt/sources.list.d/pgdg.list
RUN curl -L https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add -
RUN apt-get update
RUN apt-get install -y gcc postgresql-client-12 libpq-dev musl-dev bash postgresql-12
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y ncbi-blast+

RUN apt-get autoremove -y
RUN mkdir /pip_requirements
ADD run/requirements.txt /pip_requirements.txt
RUN pip install -r /pip_requirements.txt
COPY config.json /code/config.json

RUN mkdir -p /code
WORKDIR /code

RUN mkdir -p /data/media/db

ADD run/docker/wait_for_it.sh /wait_for_it.sh
RUN chmod +x /wait_for_it.sh

ADD . /code/
RUN ls /code

EXPOSE 8000
