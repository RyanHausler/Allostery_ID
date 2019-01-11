FROM ubuntu:latest
MAINTAINER Ryan Hausler (rmh1995@gmail.com)

ENV PACKAGES clustalw python3-tk


RUN apt-get update \
&& apt-get install -y python3-pip python3-dev \
&& cd /usr/local/bin \
&& ln -s /usr/bin/python3 python \
&& pip3 install --upgrade pip \
&& apt-get update \ 
&& apt-get upgrade -y \ 
&& useradd -ms /bin/bash luser \ 
&& printf "#!/bin/bash\nexport DEBIAN_FRONTEND=noninteractive\napt-get install -y tzdata\nln -fs /usr/share/zoneinfo/America/New_York /etc/localtime\ndpkg-reconfigure --frontend noninteractive tzdata\n" >tzinst \
&& chmod 700 tzinst \
&& ./tzinst \
&& rm -f tzinst \
&& apt-get install -y ${PACKAGES} \
&& apt-get clean

COPY . /

RUN pip install --no-cache-dir -r requirements.txt

RUN chmod -R 777 /usr/
RUN chmod -R 777 /usr/bin/
RUN chmod 777 /usr/bin/clustalw

COPY . /home/luser

USER luser

WORKDIR /home/luser

ENTRYPOINT ["python", "AllosteryID.py"]
