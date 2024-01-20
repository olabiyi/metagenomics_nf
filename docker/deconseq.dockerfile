FROM quay.io/grbot/deconseq:latest
RUN apt-get update -y 
RUN sed -E 's/db\//database\/deconseq\//g' /deconseq-standalone-0.4.3/DeconSeqConfig.pm | \
    sed -E 's/bact =>/plant =>/g' | sed -E 's/Bacterial genomes/Maize genome/g' | \
    sed -E 's/bactDB/host_prinseqDB/g' > temp && mv temp /deconseq-standalone-0.4.3/DeconSeqConfig.pm
RUN pip install awscli --upgrade
RUN export PATH=$PATH:/home/ubuntu/.local/bin
