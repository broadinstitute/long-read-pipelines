FROM us.gcr.io/broad-dsp-lrma/lr-pb:0.1.40

RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.4.0/trgt-v0.4.0-linux_x86_64.gz" -O "trgt.gz" && \
        gunzip trgt.gz && \
        chmod 777 trgt && \
        mv trgt /usr/bin/

RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.4.0/trvz-v0.4.0-linux_x86_64.gz" -O "trvz.gz" && \
        gunzip trvz.gz && \
        chmod 777 trvz && \
        mv trvz /usr/bin/
