VERSION = 2.4.0 
TAG1 = us.gcr.io/broad-dsp-lrma/lr-seqkit:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-seqkit:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
