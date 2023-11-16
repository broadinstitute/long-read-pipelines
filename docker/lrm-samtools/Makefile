VERSION = 0.1.1
TAG1 = us.gcr.io/broad-dsp-lrma/lrm-samtools:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lrm-samtools:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
