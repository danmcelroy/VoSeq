FROM postgres:12.2

ENV PGDATA /srv/postgresql/data/
ADD run/docker/postgres/init-user-db.sh /docker-entrypoint-initdb.d/init-user-db.sh
